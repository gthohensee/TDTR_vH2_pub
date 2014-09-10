function [dTss, dTpp, Mcell] = TDTR_TDEP_TTM_vH2(matparams,sysparams,...
                                           Tparams,intscheme,nnodes)
% TDTR_TDEP_TTM_vH2 - handles steady state heating, looks at pump transient 
% heating, and adjusts thermal parameters of the layers accordingly.
%
% Syntax:
% [dTss, dTpp, Mcell] = TDTR_TDEP_TTM_vH2(matparams,sysparams,...
%                                           Tparams,intscheme,nnodes)
%
% Inputs:
%    sysparams - {tau_rep f r_pump r_probe}
%    matparams - {LCTE aniso BI n_toplayer TCR}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_TTM_vH2.m for details on the params inputs.]
%
%    intscheme  - 0 = Legendre-Gauss, 1 = rombint, 2 = Simpson. For the
%                 thermal model numerical integration.
%    nnodes     - For Legendre-Gauss integration, number of nodes. Default
%                 should be at least 35.
%
% Outputs:
%    dTss - steady-state heating (T = T0 + dTss + dTpp)
%    dTpp - per-pulse heating
%    Mcell - {LCTE,LCTEG,gg}
%
% Example: 
%    --
%
% Other m-files required: TDTR_TEMP_TTM_vH2.m, SS_Heating_TTM_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_TDEP_vH2.m

% Author: Gregory Hohensee
% University of Illinois at Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Apr 2014; Last revision: 1-Apr-2014
%                          8-Apr-2014 - harmonized with non-TTM version
%                          14-July-2014 - vH2. Fixed bug in per-pulse
%                                         calculation re: absC.
%------------- BEGIN CODE --------------
%% initialize cell params
PL = [0 length(sysparams) 0 length(matparams) length(Tparams)];

% unpack matparams {Mcell aniso BI n_toplayer TCR};
if PL(4) < 1
    error('Mcell not specified in FIT_TTM.'); 
else
    Mcell = matparams{1};
    if length(Mcell) < 3
        error('Insufficient elements in Mcell for FIT_TTM.');
    else
        LCTE = Mcell{1};
        LCTEG = Mcell{2};
        gg = Mcell{3};
    end
end
%if PL(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
%if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end

% unpack pump-probe system parameters
%sysparams = {tau_rep f r_pump r_probe};
if PL(2) < 1, tau_rep = 1/80e6; 
    warning('Defaulting to 80 MHz repetition rate'); else tau_rep = sysparams{1}; end
% if PL(2) < 2
%     warning('Defaulting to 9.8 MHz modulation frequency');
%     if rorfit
%         f = [9.8e6 9.8e6]; % two modulation frequencies
%     else % simple ratio fit, one frequency
%         f = 9.8e6;
%     end
% else f = sysparams{2}; % 1x1 if not rorfit, 1x2 if rorfit.
% end 
if PL(2) < 3, error('pump spot size not specified.'); else r_pump = sysparams{3}; end
if PL(2) < 4, r_probe = r_pump;
    warning('Defaulting to r_probe = r_pump'); else r_probe = sysparams{4}; end

% unpack Tparams (TDEP check different from INITIALIZE check).
if PL(5) < 1, return; else T0 = Tparams{1}; end
if PL(5) < 2, return; else T_Mcell = Tparams{2}; end
if PL(5) < 3, return; else A_pump = Tparams{3}; end
if PL(5) < 4, A_probe = 0; warning('Defaulting to A_probe = 0.'); 
   else A_probe = Tparams{4}; end
if PL(5) < 5, absC = [295 0.13]; else absC = Tparams{5}; end
if PL(5) < 6, perpulse = 0; else perpulse = Tparams{6}; end
if PL(5) < 7, jabs = 0; else jabs = Tparams{7}; end
if PL(5) < 8, jtrans = n_toplayer+2; else jtrans = Tparams{8}; end

%% initialize remaining variables
% time-averaged laser power at sample surface, before absorption.
A_tot = A_pump + A_probe;
dTss = 0;
dTpp = 0;
if A_tot == 0, return; end

nc = length(gg(:,1)); % number of channels

T_LCTE = T_Mcell{1};
T_LCTEG = T_Mcell{2};
T_gg = T_Mcell{3};
[~,iabs] = min(abs(absC(:,1) - T0));
absC_at_T = absC(iabs,2); % should rename absC as TabsC, then use absC for this variable.


%% What is the transducer thickness, given j(abs) and j(transducer)?
if jabs ~= 0 % if there exists an absorption layer...
    tabs = LCTE(2,jabs) / LCTE(2,jtrans);
    ttrans = tabs*1e-9 + LCTE(3,jtrans); % this won't update with temperature.
else % jabs = 0: assume no absorption layer
    ttrans = LCTE(3,jtrans);
end

%% Initialize dT temperature rise(s)
dTss = SS_Heating_TTM_vH2(matparams,r_pump,r_probe,absC_at_T,A_tot,intscheme,nnodes);
if perpulse
    % factor of 2 for the 50% duty cycle of square wave pump modulation at
    % EOM. See David's RSI 2004 paper for per-pulse heating equation?
    dTpp = (2*absC_at_T*A_pump*tau_rep) / ...
           (pi * r_pump^2 * ttrans * LCTE(2,jtrans));
end

% Given: T_LCTE(i,j,m,n) = [T(m) V(n)], Vij => LCTE(i,j)
% Steps:
% 1) get individual V(n,i,j) for T(n) = T0 + dTss + dTpp; n can vary
% 2) assign V(nn,i,j) to all elements LCTE(i,j) for the various nn(i,j).
% 3) Same for LCTEG and gg.
% 4) Recompute dT until self-consistent with thermal parameters.
dTerr = dTss+dTpp; % initialize as deviation from T0
V = LCTE; % temporary values
W = LCTEG;
U = gg;
while (dTerr > 0.01) % Loop until dT converges to within 1%. Doesn't take many iterations.
    T = T0 + dTss + dTpp;
    [~,iabs] = min(abs(absC(:,1) - T));
    absC_at_T = absC(iabs,2);
    
    % Update V(i,j) from T_LCTE
    for i = 1:4
        for j = 1:length(LCTE(1,:))
            TV = T_LCTE{i,j}; % [T, V(T)] for LCTE(i,j)
            [~,iV] = min(abs(TV(:,1) - T));
            
            % No fit parameter should have any temperature
            % dependence. As such, the writeLCTE(G) databases should be
            % given the (Xijc) indicator, and return a null
            % value (negative) to prevent TDEP from adjusting the fit.
            if TV(iV,2) < 0
                V(i,j) = LCTE(i,j);
            else
                V(i,j) = TV(iV,2);
            end
            
            % % QUESTION: what about anisotropy depending on Lz(T)? % %
            % T-dep of the eta part of LCTE is assumed to be due to changes
            % in Lx/Lz, NOT just the in-plane component. For isotropic layers,
            % T_LCTE(4,j,:,2) should just be 1.
            
            % So, T_L and T_E should be constructed in the writeLCTE database
            % such that they are self-consistent: T_E = T(Lx/Lz), given T_Lx.
            
            % Therefore, I need make no compensation in eta for changes in 
            % Lz with temperature.
            
            % This applies equally to the N-channel substrate; writeLCTEG
            % should contain the necessary information in T_E.
            
            % errorbar perturbations will still have to track coupled
            % variables: eta and absorption layer.
        end
    end
    
    % (TTM) Update W(i,j) from T_LCTEG
    for i = 1:5
        for j = 1:length(LCTEG(1,:))
            TW = T_LCTEG{i,j}; % [T, W(T)] for LCTEG(i,j)
            [~,iW] = min(abs(TW(:,1) - T));
            if TW(iW,2) < 0
                W(i,j) = LCTEG(i,j);
            else
                W(i,j) = TW(iW,2);
            end
        end
    end
    
    % (TTM) Update gg(i,j) from T_gg (in case you actually know g(T))
    % Two-temperature model: the only valid g is gg(1,2).
    for i = 1:nc
        for j = i+1:nc 
            TU = T_gg{i,j}; % [T, U(T)] for gg(i,j), two-temperature model
            [~,ig] = min(abs(TU(:,1) - T));
            if TU(ig,2) < 0
                U(1,2) = gg(1,2);
            else
                U(1,2) = TU(ig,2);
            end
        end
    end
    
    Mcelltemp = {V,W,U};
    matparamstemp = matparams;
    matparamstemp{1} = Mcelltemp;
    
    % recompute dT_SS
    dTss_new = SS_Heating_TTM_vH2(matparamstemp,r_pump,r_probe,absC_at_T,A_tot,intscheme,nnodes);
    if perpulse
        dTpp_new = (2*absC_at_T*A_pump*tau_rep) / ...
                   (pi * r_pump^2 * ttrans * LCTE(2,jtrans));
    else
        dTpp_new = 0;
    end
    dT_new = dTss_new + dTpp_new;
    dT = dTss + dTpp;
    dTerr = (dT - dT_new) / (dT_new); % fractional inconsistency in dT
    
    dTss = dTss_new; % update dTss
    dTpp = dTpp_new; % update dTpp
end

% Assign final results for thermal properties at T = T0 + dTss + dTpp.
LCTE = V;
LCTEG = W;
gg = U;
Mcell = {LCTE, LCTEG, gg};
end
%------------- END CODE --------------