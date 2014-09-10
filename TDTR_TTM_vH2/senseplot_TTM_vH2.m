function [S_Mcell,S_sys] = senseplot_TTM_vH2(datparams,sysparams,calparams,...
                                             matparams,Tparams,fignum)
%senseplot_TTM_vH2 - Calculates sensitivity plots dlogR/dlogX for TTM
%thermal model. The sens_consider booleans can be edited to change which
%sensitivities are calculated. For twofit functionality, this script is
%unchanged; it is simply called twice in MANUALFIT or MAIN to produce
%two plots, hence the fignum input.
%
% Inputs:
%    datparams - {tdelay ratio_data datadir}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind rorfit twofit intscheme nnodes consider_error Mcell_err T0_err}
%    matparams - {Mcell aniso BI n_toplayer TCR}
%      Mcell     - {LCTE LCTEG gg}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%    fignum    - number for sensitivity figure
%[See INITIALIZE_CELLPARAMS_TTM_vH2 for details.]
%
% Outputs
%    S_Mcell - {S_LCTE,S_LCTEG,S_gg}
%       S_LCTE - sensitivities to one-channel overlayer parameters
%       S_LCTEG - sensitivities to TTM substrate thermal parameters
%       S_gg - sensitivities to coupling parameters
%    S_sys - {S_f, S_r_pump, S_r_probe} - system parameter sensitivities
%
% Other m-files required: TDTR_REFL_TTM_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: errorbars_TTM_vH2.m, senseplot_vH2.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 8-Apr-2014 - comments updated, harmonized with
%                                senseplot_vH1.
%                   14-July-2014 - vH2: twofit functionality
%------------- BEGIN CODE --------------
%% 
INITIALIZE_CELLPARAMS_TTM_vH2; % unpacks/cleans cellparams (the five inputs)
fprintf('Generating Sensitivity Plot for all Variables...please be patient\n')

if nargin < 6
    fignum = 202; % default fignum
end

%% Generates a reference model based on anticipated parameters
sysparams{2} = f(1);
[deltaR_model,ratio_model]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
if rorfit || twofit
    sysparams{2} = f(2); % second frequency
    [deltaR_model2,ratio_model2]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
    
    dr_model = ratio_model ./ ratio_model2; % the reference data
    sysparams{2} = f(1); % return to first frequency
end

%% initialize matrices for sensitivities
LCTEtemp = LCTE;
LCTEGtemp = LCTEG;
ggtemp = gg;

%nt = length(tdelay); 
nl = length(LCTE(1,:));
nc = length(gg(:,1)); % number of channels in substrate

S_LCTE = cell(4,nl);
S_LCTEG = cell(5,nc);
S_gg = cell(nc,nc);

% these will be written with LCTE arrays, then overwritten with LCTEG arrays.
deltaR_temp  = cell(5,max(nl,nc));
ratio_temp   = cell(5,max(nl,nc));
deltaR_temp2 = cell(5,max(nl,nc));
ratio_temp2  = cell(5,max(nl,nc));

%% which sensitivities to consider? (saves time)
LCTE_sens_consider = ones(4,nl);
LCTEG_sens_consider = ones(5,nc);
gg_sens_consider = ones(nc,nc);

% Skip absorption layer; instead, couple its perturbations into that of
% the transducer layer. If there is no transducer layer, just a "1 nm"
% absorption layer, set jabs = 0 in your analyze script.
if jabs ~= 0, LCTE_sens_consider(:,jabs) = 0; end 

for j = 1:nl
    % Skip thicknesses and heat capacities of all interfaces in LCTE
    if LCTE(2,j) == 1e5 && LCTE(3,j) == 1e-9 % reliable interface indicators
        LCTE_sens_consider(2,j) = 0;
        LCTE_sens_consider(3,j) = 0;
        LCTE_sens_consider(4,j) = 0; % interfaces are not anisotropic.
    end
    
    % skip eta for all isotropic layers
    if ~aniso(j), LCTE_sens_consider(4,j) = 0; end
end

% LCTEG, gg %
LCTEG_sens_consider(3,:) = 0; % semi-infinite substrate
for n = 1:nc
    % skip 1D and isotropic channels
    if LCTEG(4,n) == 0 || LCTEG(4,n) == 1, LCTEG_sens_consider(4,n) = 0; end
    
    % skip adiabatic interfaces
    if LCTEG(5,n) == 0, LCTEG_sens_consider(5,n) = 0; end 
    
    gg_sens_consider(n,1:n) = 0;  % skip non-unique couplings
end

%% Make sensitivity plots
%-----------------LCTE--------------
for i = 1:4
    for j=1:nl
        % skip conditions; includes the aniso variable.
        if LCTE_sens_consider(i,j) == 0, continue; end

        LCTEtemp(i,j) = LCTE(i,j)*1.01;

        % couplings between model parameters, e.g. absorption layers and
        % anisotropies.
        if jabs ~= 0 && j == jabs, LCTEtemp(i,jtrans) = LCTE(i,jtrans)*1.01; end
        if jabs ~= 0 && j == jtrans, LCTEtemp(i,jabs) = LCTE(i,jabs)*1.01; end
        if i == 1 && aniso(j), LCTEtemp(4,j) = LCTE(i,j)/1.01; end % eta and Lz are coupled.

        % Perturbing eta is assumed to be perturbing Lx, not Lz. %

        % Perform sensitivity calculation with LCTEtemp
        matparams{1} = {LCTEtemp LCTEG gg};
        sysparams = {tau_rep f(1) r_pump r_probe};
        [deltaR_temp{i,j},ratio_temp{i,j}] = TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
        if rorfit
            sysparams = {tau_rep f(2) r_pump r_probe};
            [deltaR_temp2{i,j},ratio_temp2{i,j}]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);

            dr_temp = ratio_temp{i,j} ./ ratio_temp2{i,j};
            Num=log(dr_temp)-log(dr_model);
        else
            Num=log(ratio_temp{i,j})-log(ratio_model);
        end
        Denom=log(LCTEtemp(i,j))-log(LCTE(i,j)); % log(1.01)
        S_LCTE{i,j}=Num/Denom;

        LCTEtemp = LCTE;
        matparams{1} = Mcell;
        fprintf('Calculated S_LCTE{%i,%i}...\n',i,j);
    end
end
clear LCTEtemp;

%% LCTEG sensitivities
%-----------------LCTEG--------------
for i = 1:5
    for j=1:nc
        if LCTEG_sens_consider(i,j) == 0, continue; end
        
        LCTEGtemp(i,j) = LCTEG(i,j)*1.01;
        
        % eta and Lz are coupled
        if i == 1 && aniso(j), LCTEGtemp(4,j) = LCTEGtemp(4,j)/1.01; end
        
        % Perturbing eta is assumed to be perturbing Lx, not Lz. %

        % Perform sensitivity calculation with LCTEGtemp
        matparams{1} = {LCTE LCTEGtemp gg};
        sysparams = {tau_rep f(1) r_pump r_probe};
        [deltaR_temp{i,j},ratio_temp{i,j}] = TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
        if rorfit
            sysparams = {tau_rep f(2) r_pump r_probe};
            [deltaR_temp2{i,j},ratio_temp2{i,j}]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);

            dr_temp = ratio_temp{i,j} ./ ratio_temp2{i,j};
            Num=log(dr_temp)-log(dr_model);
        else
            Num=log(ratio_temp{i,j})-log(ratio_model);
        end
        Denom=log(LCTEGtemp(i,j))-log(LCTEG(i,j)); % log(1.01)
        S_LCTEG{i,j}=Num/Denom;
        
        LCTEGtemp = LCTEG;
        matparams{1} = Mcell;
        sysparams{2} = f(1);
        fprintf('Calculated S_LCTEG{%i,%i}...\n',i,j);
    end
end
clear LCTEGtemp;

%% gg sensitivities
for i = 1:nc
    for j = i+1:nc % only consider unique couplings
        if gg_sens_consider(i,j) == 0, continue; end
        ggtemp(i,j) = gg(i,j)*1.01;
        
        % Perform sensitivity calculation with ggtemp
        matparams{1} = {LCTE,LCTEG,ggtemp};
        sysparams = {tau_rep f(1) r_pump r_probe};
        [deltaR_temp{i,j},ratio_temp{i,j}] = TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
        if rorfit
            sysparams = {tau_rep f(2) r_pump r_probe};
            [deltaR_temp2{i,j},ratio_temp2{i,j}]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);

            dr_temp = ratio_temp{i,j} ./ ratio_temp2{i,j};
            Num=log(dr_temp)-log(dr_model);
        else
            Num=log(ratio_temp{i,j})-log(ratio_model);
        end
        Denom=log(ggtemp(i,j))-log(gg(i,j)); % log(1.01)
        S_gg{i,j}=Num/Denom;
        
        ggtemp = gg;
        matparams{1} = Mcell;
        sysparams{2} = f(1);
        fprintf('Calculated S_gg{%i,%i}...\n',i,j);
    end
end
S_g = S_gg{1,2}; % restrict to 2-channel case, one unique coupling "g".

%% System sensitivities
sys_consider = [1 2 3];
S_f = []; S_r_pump = []; S_r_probe = [];
for i = sys_consider;
    sysbump = 1.01;
    
    switch i
        case 1, sysparams{2} = f(1)*sysbump;
        case 2, sysparams{3} = r_pump*sysbump;
        case 3, sysparams{4} = r_probe*sysbump;
    end
    [temp_deltaR,temp_ratio]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
    if rorfit
        sysparams{2} = f(2)*sysbump;
        [temp_deltaR2,temp_ratio2]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);

        dr_temp = temp_ratio ./ temp_ratio2;
        Num=log(dr_temp)-log(dr_model);
    else
        Num=log(temp_ratio)-log(ratio_model);
    end
    Denom=log(sysbump);
    
    switch i
        case 1, S_f=Num/Denom;
        case 2, S_r_pump=Num/Denom;
        case 3, S_r_probe=Num/Denom;
    end
    sysparams = {tau_rep, f(1), r_pump, r_probe}; % return to reference value
    fprintf('Calculated S_sys #%i...\n',i);
end

%% Plot sensitivities
figure(fignum)
clf
axes('XScale','log');
set(gca,'Box','on');
hold on;
ColorOrder = get(gcf,'DefaultAxesColorOrder');

%% label LCTE sensitivities
LCTElegend = [];
LCTElab = cell(size(LCTE));
LCTEmarker = {'o','*','x','+'};

for i = 1:4
    for j = 1:nl
        switch i
            case 1, LCTElab{i,j} = sprintf('L-%i',j);
            case 2, LCTElab{i,j} = sprintf('C-%i',j);
            case 3, LCTElab{i,j} = sprintf('t-%i',j);
            case 4, LCTElab{i,j} = sprintf('e-%i',j);
        end
        
        if isempty(S_LCTE{i,j}), continue;
        else
            semilogx(tdelay,S_LCTE{i,j},LCTEmarker{i},'Color',ColorOrder(j,:));
            LCTElegend = [LCTElegend;LCTElab{i,j}];
        end
    end
end

%% label and plot LCTEG sensitivities
LCTEGlegend = [];
LCTEGlab = cell(size(LCTEG));
LCTEGmarker = {'s','*','x','+','^'};

for i = 1:5
    for j = 1:nc
        switch i
            case 1, LCTEGlab{i,j} = sprintf('L%i*',j);
            case 2, LCTEGlab{i,j} = sprintf('C%i*',j);
            case 3, LCTEGlab{i,j} = sprintf('t%i*',j);
            case 4, LCTEGlab{i,j} = sprintf('e%i*',j);
            case 5, LCTEGlab{i,j} = sprintf('G%i*',j);
        end
        
        if isempty(S_LCTEG{i,j}), continue;
        else
            if i == 5
                semilogx(tdelay,S_LCTEG{i,j},LCTEGmarker{i},'Color',ColorOrder(7-j,:),...
                                                            'MarkerFaceColor',ColorOrder(7-j,:));
            else
                semilogx(tdelay,S_LCTEG{i,j},LCTEGmarker{i},'Color',ColorOrder(7-j,:));
            end
            LCTEGlegend = [LCTEGlegend;LCTEGlab{i,j}];
        end
    end
end

%% plot gg sensitivities
if ~isempty(S_g)
    gglegend = 'g12'; % for TTM
    semilogx(tdelay,S_g,'k--','LineWidth',2); % big dashed line
else
    gglegend = '';
end

%% plot and label system sensitivities
syslegend = [];
for i = sys_consider
    switch i
        case 1, semilogx(tdelay,S_f,'-','Color',ColorOrder(i,:));
            syslegend = [syslegend;'frq'];
        case 2, semilogx(tdelay,S_r_pump,'-','Color',ColorOrder(i,:));
            syslegend = [syslegend;'Rpm'];
        case 3, semilogx(tdelay,S_r_probe,'-','Color',ColorOrder(i,:));
            syslegend = [syslegend;'Rpb'];
    end
end

%% Other plot details
legend([LCTElegend;LCTEGlegend;gglegend;syslegend])
xlabel('time delay (ps)','FontSize',16)
if rorfit
    ylabel(sprintf('Sensitivity:  dlog[R(%0.2f MHz)/R(%0.2f MHz)/dlogX',f(1)/1e6,f(2)/1e6), 'FontSize',16)
else
    ylabel(sprintf('Sensitivity:  dlogR/dlogX @ %0.2f MHz', f/1e6), 'FontSize',16)
end
set(gca, 'XTick', [1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
set(gca, 'XTickLabel', [100, 200, 500, 1000, 2000, 5000, 1e4]);
set(gca, 'XMinorTick', 'off');
set(gca,'FontSize',16);
axis([100e-12 10e-9 -2 2]);

%% export sensitivities
S_Mcell = {S_LCTE,S_LCTEG,S_gg};
S_sys = {S_f,S_r_pump,S_r_probe};
end
