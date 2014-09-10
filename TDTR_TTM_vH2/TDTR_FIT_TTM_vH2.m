function [Z,deltaR_model,ratio_model,Mcell,T_adj,dr_model]= ...
TDTR_FIT_TTM_vH2(X, Xijc, datparams,sysparams,calparams,matparams,Tparams)
%TDTR_FIT_TTM_vH2 - automatically fit TTM thermal model to TDTR ratio data.
%The main program tries to minimize "Z" by optimizing the variable(s) X.
% This version is capable of self-consistent temperature adjustment using
% the Tparams. If you do not wish to do self-consistent T-corrections,
% simply enter the one-element cell {-1} as Tparams.
%
% Syntax:  [Z,deltaR_model,ratio_model,Mcell,T_adj,dr_model]= ...
%TDTR_FIT_TTM_vH2(X, Xijc, datparams,sysparams,calparams,matparams,Tparams)
%
% Inputs:
%    X        - Fit parameters. Has to be separate from everything else
%             - for fminsearch functionality.
%    Xijc     - Mx3 matrix: mth row is for X(m), representing either
%               LCTE(i,j), LCTEG(i,j), or gg(i,j), depending on c.
%    datparams - {tdelay ratio_data datadir}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind rorfit twofit intscheme nnodes consider_error Mcell_err T0_err}
%    matparams - {Mcell aniso BI n_toplayer TCR}
%      Mcell     - {LCTE LCTEG gg}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_TTM_vH2.m for details on the params inputs.]
%
% Outputs:
%    Z           - Goodness-of-fit. In Joe's words: "Typically, this is 
%                  the sum of the squares of the residuals, but you might 
%                  want to weight these by the sensitivity, particularly 
%                  if you don't intend to calculate the errorbars!"
%    deltaR      - Complex number array. Real part is the model V(in),
%                  imaginary part is the model V(out). Represents change
%                  in reflectance from pump heating.
%                  If twofit is true, this is a 2-element cell array.
%    ratio_model - Ratio -V(in)/V(out) from TDTR thermal model.
%                  If twofit is true, this is a 2-element cell array.
%    Mcell       - Self-consistent laser heating can change the thermal
%                  parameters depending on the fit parameters, so the
%                  final Mcell is also an output of this function.
%    T_adj       - T_adj = T0 + dTss + dTpp, the "actual" temperature
%                  in Kelvin adjusted for steady-state (and per-pulse
%                  if perpulse is TRUE) heating.
%    dr_model    - model data for ratio-over-ratio fitting (rorfit).
%
% Example:
%    Xsol = fminsearch(@(X),TDTR_FIT_TTM_vH2(X, Xijc, datparams,...
%               sysparams, calparams, matparams, Tparams),Xguess,options);
%
% Other m-files required: TDTR_REFL_TTM_vH2.m, TDTR_TEMP_TTM_vH2.m, 
%                         TDTR_TDEP_TTM_vH2.m, SS_Heating_TTM_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_MANUALFIT_TTM_vH2.m, TDTR_MANUALFIT_vH2.m, TDTR_FIT_vH2.m

% Author: Gregory Hohensee
% Acknowledgement: built from TDTR_FIT_V4B, my bi-directional tweak to the
% original TDTR_FIT_V4.m, by Joseph P. Feser.
% University of Illinois at Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Apr 2014; Last revision: 7-Apr-2014 - vH1
%                          14-Jul-2014 - vH2: twofit functionality.

%------------- BEGIN CODE --------------
%% Check input parameters, assign defaults, errors, warnings as necessary
INITIALIZE_CELLPARAMS_TTM_vH2; % unpacks and re-packs cellparams

%% Assign fit variables X according to their Xijc index.
% Xijc = [i,j,c] indexes X across (i,j) of the Mcell{c}.
nrows = length(Xijc(:,1));
for x = 1:nrows % iterate over fit parameters
    c = Xijc(x,3);
    Xij = Xijc(:,1:2);
    switch c % switch over Mcell structures
        case 1, LCTE(Xij(x,1),Xij(x,2)) = X(x); 
        case 2, LCTEG(Xij(x,1),Xij(x,2)) = X(x); 
        case 3, gg(Xij(x,1),Xij(x,2)) = X(x);
    end
end

%% Perturbations for error bar calculation.
if consider_error
    % Note that only one element of LCTE_err should be nonzero at any
    % one call of this function.
    LCTE = LCTE .* (1 + LCTE_err);
    LCTEG = LCTEG .* (1 + LCTEG_err);
    gg = gg .* (1 + gg_err);
    
    if jabs ~= 0 % if there exists an absorption layer...
        % its uncertainty is coupled to that of the transducer layer.
        if sum(LCTE(1:2,jabs) == 0) == 0
            LCTE(1:2,jabs) = LCTE(1:2,jabs) .* (1 + LCTE_err(1:2,jtrans));
        end
    end
    
    % update eta for a perturbation in cross-plane conductivity,
    % if the layer has anisotropy that can vary.
    LCTE(4,:) = LCTE(4,:) ./ (1 + aniso.*LCTE_err(1,:));
    LCTEG(4,:) = LCTEG(4,:) ./ (1 + aniso.*LCTEG_err(1,:));
    
    % measured temperature or pressure uncertainty
    T0 = T0 + T0_err; % perturbation will take effect in TDEP
    %P0 = P0 + P0_err; % perturbation... needs writeLCTE(P).
    
end
Mcell = {LCTE,LCTEG,gg}; % update Mcell after error propagation
matparams{1} = Mcell; % update matparams before TDEP

%% Self-consistent steady-state (and optionally per-pulse) heating 
if T0 ~= -1
    % re-use Tparams for the TDEP function, although now Tparams
    % may have some default values in it that were unspecified before.
    Tparams{1} = T0; % update in case of consider_error
    [dTss, dTpp, Mcell] = TDTR_TDEP_TTM_vH2(matparams,sysparams,...
                                           Tparams,intscheme,nnodes);
    fprintf('T0 = %0.2f K, dTss = %0.2f, dTpp = %0.2f\n',T0,dTss,dTpp)
end
matparams{1} = Mcell; % update matparams after TDEP
%% Compute model; plot model vs. data, and return goodness of fit.
sysparams{2} = f(1); % first frequency
[deltaR_model,ratio_model]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,...
                                            A_pump,intscheme,nnodes);
dr_model = 0; % default
if rorfit || twofit
    sysparams{2} = f(2); % second frequency
    [~,ratio_model2]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,...
                                        A_pump,intscheme,nnodes);
    sysparams{2} = f; % return to reference
    dr_model = ratio_model ./ ratio_model2;
end
    
%% Uncomment next few lines to see the non-linear optimization in action!
figure(10)
if rorfit
    semilogx(tdelay,dr_data,'ok',tdelay,dr_model,'r'); 
    
    ylabel('R_1 / R_2','FontSize',16);
    axis([1e-10 10e-9 0 max(1.2*max(dr_data),2)]);
elseif twofit
    semilogx(tdelay,ratio_data, 'ok', tdelay,ratio_model,'r',...
             tdelay,ratio_data2,'or',tdelay,ratio_model2,'b');  

    ylabel('R_1 and R_2','FontSize',16);
    axis([1e-10 10e-9 min([ratio_data;ratio_model;ratio_data2;ratio_model2])*0.8 ...
                      max([ratio_data;ratio_model;ratio_data2;ratio_model2])*1.2]);
    legend('D1','M1','D2','M2');    
else % standard ratio fit
    loglog(tdelay,ratio_data,'ob',tdelay,ratio_model,'r');
    ylabel('Ratio','FontSize',16);
    axis([1e-10 10e-9 min([ratio_data;ratio_model])/2 max([ratio_data;ratio_model])*2])    
    set(gca, 'YTick', [0.1 0.2 0.5, 1, 2, 5, 10, 20, 50, 100, 200]);
    set(gca, 'YTickLabel', [0.1 0.2 0.5, 1, 2, 5, 10, 20, 50, 100, 200]);
end
xlabel('Time delay (ps)','FontSize',16);
set(gca, 'XTick', [1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
set(gca, 'XTickLabel', [100, 200, 500, 1000, 2000, 5000, 1e4]);
set(gca, 'XMinorTick', 'off');
set(gca, 'YMinorTick', 'off')
set(gca, 'FontSize', 16);
pause(0.1) % may be necessary so MATLAB can refresh the figure between iterations

% Define symbols and units for all possible fit parameters
tag = {'L','C','t','eta','L*','C*','t*','eta*','G*','g'};
units = {'W/m-K', 'J/cm^3-K', 'nm', '(Lx/Lz)',...
         'W/m-K', 'J/cm^3-K', 'mm', '(Lx/Lz)','MW/m^3-K',...
         'pW/nm^3-K'};
scale = [1 1e-6 1e9 1 1 1e-6 1e3 1 1e-3 1];

% Report current fit parameters
Mat = sprintf('Current material fit parameters (*: %i-channel substrate)',length(LCTEG(1,:)));
for x = 1:nrows % iterate through fit parameters
    c = Xijc(x,3);
    switch c
        case 3 % gg
            itag = 10;
            Mat = char(Mat,sprintf('%s(%i,%i) = %0.3f %s',...
                  tag{itag},Xijc(x,1),Xijc(x,2),...
                  gg(Xijc(x,1),Xijc(x,2))/scale(itag),units{itag}));
        otherwise % LCTE or LCTEG
            itag = Xijc(x,1) + 4*(c-1);
            if itag <= 4 % LCTE parameters
                Mat = char(Mat,sprintf('%s(%i) = %0.4f %s',...
                  tag{itag},Xijc(x,1),...
                  LCTE(Xijc(x,1),Xijc(x,2))/scale(itag),units{itag}));
            else % itag between 5 and 9, LCTEG
                Mat = char(Mat,sprintf('%s(%i) = %0.4f %s',...
                  tag{itag},Xijc(x,1),...
                  LCTEG(Xijc(x,1),Xijc(x,2))/scale(itag),units{itag}));
            end
    end
end
Mat

%% Assign final temperature adjusted for SS and PP heating.
if T0 ~= -1, T_adj = T0 + dTss + dTpp; else T_adj = T0; end

%% Goodness-of-fit by fractional residuals, starting from Zdelay ps
if rorfit
    Lr = length(dr_model);
    res=(1-(dr_model(Zind:Lr) ...
            ./dr_data(Zind:Lr)) ).^2;
elseif twofit % res = (1 - (R1_model / R1_data))^2 + (1 - (R2_model / R2_data))^2
    L1 = length(ratio_model);
    L2 = length(ratio_model2);
    res=(1-(ratio_model(Zind:L1)./ratio_data(Zind:L1))).^2 ...
            + (1-(ratio_model2(Zind:L2) ./ratio_data2(Zind:L2))).^2;
        
    % combine models into cell arrays for convenient output.
    ratio_model = {ratio_model,ratio_model2};
    deltaR_model = {deltaR_model,deltaR_model2};
else % ordinary ratio fitting
    L1 = length(ratio_model);
    res=(1-(ratio_model(Zind:L1) ...
            ./ratio_data(Zind:L1) )).^2;
end
Z=sum(res)
end
