function [S_LCTE,S_sys] = senseplot_vH2(datparams,sysparams, calparams, matparams, Tparams)
%senseplot_vH2 - Calculates sensitivity plots dlogR/dlogX for thermal
%model. The sens_consider booleans can be edited to change which
%sensitivities are calculated.
%
% Inputs:
%    datparams - {tdelay ratio_data datadir}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind rorfit intscheme nnodes consider_error Mcell_err T0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_vH2 for details.]
%
% Outputs
%    S_LCTE - sensitivities to one-channel overlayer parameters
%    S_sys  - {S_f,S_r_pump,S_r_probe};
%       S_f, S_r_pump, S_r_probe - system parameter sensitivities
%
% Other m-files required: TDTR_REFL_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: errorbars_vH2.m, senseplot_TTM_vH2.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 8-Apr-2014 - made into function, comments updated,
%                                harmonized with TTM version.
%                   14-July-2014 - vH2. No changes since June 11th.
%------------- BEGIN CODE --------------
%%
INITIALIZE_CELLPARAMS_vH2; % unpacks/cleans cellparams (the five inputs)
fprintf('Generating Sensitivity Plot for all Variables...please be patient\n')

%% Generates a reference model based on anticipated parameters
[deltaR_model,ratio_model]=TDTR_REFL_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);

switch sigfit
    case 1 % V(in) fit
        % Construct normalized V(in) model and data,
        % relative to its value at Zdelay picoseconds.
        Vin_model = real(deltaR_model);
        Vin_model_Zdelay = Vin_model(Zind);
        nVin_model = Vin_model / Vin_model_Zdelay;

        Vin_data_Zdelay = Vin_data(Zind);
        nVin_data = Vin_data / Vin_data_Zdelay;
    case 2 % V(out) fit
        % Construct normalized V(out) model and data,
        % relative to its mean value near Zdelay picoseconds.
        Vout_model = imag(deltaR_model);
        Vout_model_Zdelay = mean(Vout_model(Zind:Zind+3));
        nVout_model = Vout_model / Vout_model_Zdelay;

        Vout_data_Zdelay = mean(Vout_data(Zind:Zind+3));
        nVout_data = Vout_data / Vout_data_Zdelay;
    otherwise % ratio fit
        % do nothing here
end

%% initialize matrices for sensitivities
LCTEtemp = LCTE;

%nt = length(tdelay); 
nl = length(LCTE(1,:));

S_LCTE = cell(4,nl);

deltaR_temp  = cell(4,nl);
ratio_temp   = cell(4,nl);

%% which sensitivities to consider? (saves time)
LCTE_sens_consider = ones(4,nl);

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
jabs
LCTE_sens_consider

%% -----------------Compute sensitivities for LCTE--------------
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
        matparams{1} = LCTEtemp;
        [deltaR_temp{i,j},ratio_temp{i,j}] = TDTR_REFL_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
        switch sigfit
            case 1
                Vin_temp = real(deltaR_temp{i,j}); 
                norm = Vin_temp(Zind);
                nVin_temp = Vin_temp / norm;
                Num=log(nVin_temp)-log(nVin_model);
            case 2
                Vout_temp = imag(deltaR_temp{i,j}); 
                norm = Vout_temp(Zind);
                nVout_temp = Vout_temp / norm;
                Num=log(nVout_temp)-log(nVout_model);
            otherwise
                Num=log(ratio_temp{i,j})-log(ratio_model);
        end
        Denom=log(LCTEtemp(i,j))-log(LCTE(i,j));
        S_LCTE{i,j}=Num/Denom;

        LCTEtemp = LCTE;
        matparams{1} = LCTE;
        fprintf('Calculated S_LCTE{%i,%i}...\n',i,j);
    end
end
clear LCTEtemp;

%% -----------------Compute sensitivities for system parameters-----------
sys_consider = [1 2 3];
S_f = []; S_r_pump = []; S_r_probe = [];
for i = sys_consider
    sysbump = 1.01;
    switch i
        case 1, sysparams{2} = f*sysbump;
        case 2, sysparams{3} = r_pump*sysbump;
        case 3, sysparams{4} = r_probe*sysbump;
    end
    
    [temp_deltaR,temp_ratio]=TDTR_REFL_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
    switch sigfit
        case 1
            Vin_temp = real(temp_deltaR); 
            norm = Vin_temp(Zind);
            nVin_temp = Vin_temp / norm;
            Num=log(nVin_temp)-log(nVin_model);
        case 2
            Vout_temp = imag(temp_deltaR); 
            norm = Vout_temp(Zind);
            nVout_temp = Vout_temp / norm;
            Num=log(nVout_temp)-log(nVout_model);
        otherwise
            Num=log(temp_ratio)-log(ratio_model);
    end
    Denom=log(sysbump);
    
    switch i
        case 1, S_f=Num/Denom;
        case 2, S_r_pump=Num/Denom;
        case 3, S_r_probe=Num/Denom;
    end
    sysparams = {tau_rep, f, r_pump, r_probe}; % return to reference value
    fprintf('Calculated S_sys #%i...\n',i);
end

%% Plot sensitivities
figure(202)
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
            case 1, LCTElab{i,j} = sprintf('L%i',j);
            case 2, LCTElab{i,j} = sprintf('C%i',j);
            case 3, LCTElab{i,j} = sprintf('t%i',j);
            case 4, LCTElab{i,j} = sprintf('e%i',j);
        end
        
        if isempty(S_LCTE{i,j}), continue;
        else
            semilogx(tdelay,S_LCTE{i,j},LCTEmarker{i},'Color',ColorOrder(j,:));
            LCTElegend = [LCTElegend;LCTElab{i,j}];
        end
    end
end

%% plot and label system sensitivities
syslegend = [];
for i = sys_consider
    switch i
        case 1, semilogx(tdelay,S_f,'-','Color',ColorOrder(i,:));
            syslegend = [syslegend;'ff'];
        case 2, semilogx(tdelay,S_r_pump,'-','Color',ColorOrder(i,:));
            syslegend = [syslegend;'Rp'];
        case 3, semilogx(tdelay,S_r_probe,'-','Color',ColorOrder(i,:));
            syslegend = [syslegend;'Rb'];
    end
end

%% Other plot details
figure(202);
set(gca,'FontSize',16);
set(gca, 'TickLength' , [.02 .02]);

legend([LCTElegend;syslegend])
xlabel('time delay (ps)','FontSize',16)
switch sigfit
    case 1, ylabel(sprintf('Sensitivity:  dlog[nV(in) @ %0.0f ps]/dlogX',tdelay(Zind)*1e12), 'FontSize',16);
    case 2, ylabel(sprintf('Sensitivity:  dlog[nV(out) @ %0.0f ps]/dlogX',tdelay(Zind)*1e12), 'FontSize',16);
    otherwise, ylabel('Sensitivity:  dlogR/dlogX', 'FontSize',16);

set(gca, 'XTick', [1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
set(gca, 'XTickLabel', [100, 200, 500, 1000, 2000, 5000, 1e4]);
set(gca, 'XMinorTick', 'off');
axis([100e-12 10e-9 -2 2]);

%% export sensitivities
%S_LCTE = S_LCTE;
S_sys = {S_f,S_r_pump,S_r_probe};
end
