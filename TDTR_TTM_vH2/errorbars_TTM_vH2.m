function [kErr_perc_TTM, kErr_abs_TTM, ErrSummary_perc_TTM, ErrSummary_abs_TTM] = ...
    errorbars_TTM_vH2(XguessIJC,datparams,sysparams,calparams,matparams,Tparams,options)
%errorbars_TTM_vH2 - calculates error bars for the Xijc fit parameters
%based on uncertainties in other thermal/system parameters. Sequentially
%calls TDTR_FIT_TTM_vH2.m to calculate how much the fitted Xsol changes
%with a perturbation in another parameter equal to its uncertainty. Sums up
%changes in Xsol in quadrature for all considered perturbations, as
%specified by errorbar_conditions_TTM_vH2.m script.
%
% Inputs:
%    Xijc      - XguessIJC(:,2:4). For each row, [i j c] indicates the
%                (i,j)th element of Mcell{c}.
%    datparams - {tdelay ratio_data datadir}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind rorfit twofit intscheme nnodes consider_error Mcell_err T0_err}
%    matparams - {Mcell aniso BI n_toplayer TCR}
%      Mcell     - {LCTE LCTEG gg}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%    options   - fminsearch options, precision tolerances on the fit.
%
% Outputs:
%    kErr_perc_TTM - total fractional errorbars from quadrature sum over all
%                    perturbations from uncertainties. Each row corresponds
%                    to the same row in Xguess / the fit parameters.
%    kErr_abs_TTM  - Absolute error bars.
%    ErrSummary_perc_TTM - individual fractional errors from individual
%                          perturbations in other material / system
%                          parameters.
%    ErrSummary_abs_TTM - individual absolute errors in Xsol due to
%                         uncertainty in other parameters.
%
% Other m-files required: TDTR_FIT_TTM_vH2.m, TDTR_REFL_TTM_vH2.m,
%                         TDTR_TDEP_TTM_vH2.m, SS_Heating_TTM_vH2.m,
%                         TDTR_TEMP_TTM_vH2.m.
% Subfunctions: none
% MAT-files required: none
%
% See also: errorbars_vH2.m, senseplot_TTM_vH2.m

% Author: Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 7-Apr-2014 - converted to function
%                   8-Apr-2014 - harmonized with errorbars_vH1.
%                   14-July-2014 - vH2, now with twofit.
%------------- BEGIN CODE --------------

%% unpack inputs
INITIALIZE_CELLPARAMS_TTM_vH2; % unpacks/cleans the five cellparams

%% Initialize variables and uncertainties
Xijc = XguessIJC(:,2:4);
Xguess = XguessIJC(:,1);
errorbar_conditions_TTM_vH2; % default uncertainties
                
%% Compute reference fit parameters Xsol
% TolX = 1e-2 reduces strictness of fminsearch
if nargin < 7, options = optimset('TolX',1e-2); end 

fprintf('YErr(n,:) = uncertainty (absolute) in X due to uncertainty in parameter Y(n)\n')
Xsol=fminsearch(@(X) TDTR_FIT_TTM_vH2(X,Xijc,datparams,sysparams,calparams,...
                                  matparams,Tparams),Xguess,options);
Xguess = Xsol; % these are interchangeable in this script from here onward.
%% Compute change in fit parameters due to independent perturbations
% Errors from LCTE
for i = 1:4
    for j=1:nl
        % skip conditions
        if LCTE_consider(i,j) == 0, continue; end
        if j == jabs, continue; end % skip error propagation from absorption layer
                                    % to avoid double-counting; the LCTE_err
                                    % elements for the absorption layer are
                                    % applied simultaneously with the
                                    % LCTE_err elements for the transducer.
                                    % If you have a "1nm transducer", where
                                    % the entire transducer is modeled as
                                    % an absorption layer, you should set
                                    % jabs = 0 in your analyze script.
        
        LCTE_err_temp(i,j) = LCTE_err(i,j);
            
        % couplings between model parameters, e.g. absorption layer and
        % anisotropies.
        %if jabs ~= 0 && j == jabs, LCTE_err_temp(i,jtrans) = LCTE_err(i,jtrans); end
        if jabs ~= 0 && j == jtrans, LCTE_err_temp(i,jabs) = LCTE_err(i,jabs); end
        if i == 1, LCTE_err_temp(4,j) = LCTE_err(i,j); end % eta and Lz are coupled.
        
        % Compute change in X due to uncertainties in LCTE
        matparams{1} = LCTE_err_temp;
        Xsoltemp=fminsearch(@(X) TDTR_FIT_TTM_vH2(X,Xijc,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
        
        LCTE_abs_err{i,j}=abs(Xsoltemp-Xguess); %Errors in X(:) due to variable LCTE(i,j)
        
        LCTE_err_temp = zeros(4,nl); % re-initialize LCTE_err_temp
        sprintf('Done with LCTE(%i,%i) error...',i,j)
    end
end
matparams{1} = Mcell; % re-initialize matparams.

%% Errors from LCTEG
for i = 1:5
    for j=1:nc
        % skip conditions
        if LCTEG_consider(i,j) == 0, continue; end
        
        LCTEG_err_temp(i,j) = LCTEG_err(i,j);
            
        % couplings between model parameters, e.g. anisotropies.
        % Isotropic or low-dimensional channels are already filtered out
        % by LCTEG_consider.
        if i == 1, LCTEG_err_temp(4,j) = LCTEG_err(i,j); end % eta and Lz are coupled.
        
        % Compute change in X due to uncertainties in LCTEG
        matparams{1} = {LCTE, LCTEG_err_temp, gg};
        Xsoltemp=fminsearch(@(X) TDTR_FIT_TTM_vH2(X,Xijc,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
        
        LCTEG_abs_err{i,j}=abs(Xsoltemp-Xguess); %Errors in X(:) due to variable LCTEG(i,j)
        
        LCTEG_err_temp = zeros(5,nc); % re-initialize LCTEG_err_temp
        sprintf('Done with LCTEG(%i,%i) error...',i,j)
    end
end
matparams{1} = Mcell; % re-initialize matparams.

%% Errors from gg
for i = 1:nc
    for j=1+1:nc
        % skip conditions
        if gg_consider(i,j) == 0, continue; end
        
        gg_err_temp(i,j) = gg_err(i,j);
        
        % Compute change in X due to uncertainties in gg
        matparams{1} = {LCTE, LCTEG, gg_err_temp};
        Xsoltemp=fminsearch(@(X) TDTR_FIT_TTM_vH2(X,Xijc,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
        
        gg_abs_err{i,j}=abs(Xsoltemp-Xguess); %Errors in X(:) due to variable gg(i,j)
        
        gg_err_temp = zeros(size(gg)); % re-initialize gg_err_temp
        sprintf('Done with gg(%i,%i) error...',i,j)
    end
end
matparams{1} = Mcell; % re-initialize matparams.

%% Errors from sys parameters
%-------Probe Radius--------------
if r_probe_consider==1
    sysparams{4} =r_probe*(1+r_err);
    Xsoltemp=fminsearch(@(X) TDTR_FIT_TTM_vH2(X,Xijc,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);

    for n=1:length(Xguess)
        r_probe_err(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_probe(ii)
    end
    sprintf('Done with r_probe error...')
    r_probe_err
    
    sysparams{4} = r_probe; % return to reference value
end

%-------Pump Radius--------------
if r_pump_consider==1
    sysparams{3} =r_pump*(1+r_err);
    Xsoltemp=fminsearch(@(X) TDTR_FIT_TTM_vH2(X,Xijc,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
    for n=1:length(Xguess)
        r_pump_err(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_pump(ii)
    end
    sprintf('Done with r_pump error...')
    r_pump_err
    
    sysparams{3} = r_pump; % return to reference value
end
%--------Phase Error-------------
if phase_consider==1
    % for rorfit, assumes equal degphase for each data set.
    % I expect r-over-r analysis greatly reduces the influence
    % of phase error on any fit parameters.
    radphase=pi/180*degphase;
    Vtemp=(Vin_data+sqrt(-1)*Vout_data)*exp(sqrt(-1)*radphase);
    Vin_phaseshifted=real(Vtemp);
    Vout_phaseshifted=imag(Vtemp);
    ratio_phaseshifted=-Vin_phaseshifted./Vout_phaseshifted;
    
    if rorfit || twofit
        Vtemp2=(Vin_data2+sqrt(-1)*Vout_data2)*exp(sqrt(-1)*radphase);
        Vin_phaseshifted2=real(Vtemp2);
        Vout_phaseshifted2=imag(Vtemp2);
        ratio_phaseshifted2=-Vin_phaseshifted2./Vout_phaseshifted2;
        
        dr_phaseshifted = ratio_phaseshifted ./ ratio_phaseshifted2;
        if rorfit
            datparams{2} = dr_phaseshifted;
        else % twofit
            datparams{2} = {ratio_phaseshifted,ratio_phaseshifted2};
        end
    else % ordinary fit
        datparams{2} = ratio_phaseshifted;
    end
    
    Xsoltemp=fminsearch(@(X) TDTR_FIT_TTM_vH2(X,Xijc,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
    for n=1:length(Xguess)
        phase_err(1,n)=abs(Xsoltemp(n)-Xguess(n));
    end
    sprintf('Done with phase error...')
    phase_err
    
    % return to reference value for signal data "datparams{2}"
    if rorfit
        datparams{2} = dr_data; 
    elseif twofit
        datparams{2} = {ratio_data,ratio_data2};
    else
        datparams{2} = ratio_data;
    end
end

%% Assemble error summaries
ErrSummary_abs_LCTE = [];
ErrSummary_abs_LCTEG = [];
ErrSummary_abs_gg = [];
for i = 1:4
    for j = 1:nl
        % empty cells of LCTE_err{i,j} will have no effect here
        ErrSummary_abs_LCTE = vertcat(ErrSummary_abs_LCTE,LCTE_err{i,j});
    end
end
for i = 1:5
    for j = 1:nc
        ErrSummary_abs_LCTEG = vertcat(ErrSummary_abs_LCTEG,LCTEG_err{i,j});
    end
end
for i = 1:nc
    for j = i+1:nc
        ErrSummary_abs_gg = vertcat(ErrSummary_abs_gg,gg_err{i,j});
    end
end
ErrSummary_abs_sys = vertcat(r_probe_err,r_pump_err,phase_err);

%% Assemble reports for error summaries
ErrSummary_abs_TTM = vertcat(ErrSummary_abs_LCTE,...
                             ErrSummary_LCTEG,...
                             ErrSummary_abs_gg,...
                             ErrSummary_abs_sys);

repeat_Xsol = ones(length(ErrSummary_abs_TTM(:,1)),1)*Xsol; % rows for parameters and 3 spacers
ErrSummary_perc_TTM = ErrSummary_abs_TTM ./ repeat_Xsol     %percent error broken by variable

kErr_perc_TTM=sqrt(sum(ErrSummary_perc_TTM.^2,1)) %total percent error in each fitted parameter
kErr_abs_TTM=kErr_perc_TTM.*Xsol                  %total absolute error in each fitted parameter