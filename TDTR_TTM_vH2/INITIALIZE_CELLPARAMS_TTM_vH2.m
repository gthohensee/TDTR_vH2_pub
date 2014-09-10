%INITIALIZE_CELLPARAMS_TTM_vH2 - Unpacks cell arrays into variables for the
% % FIT, MANUALFIT, errorbar, and senseplot functions. If the cell
% arrays are incomplete (not as many elements as expected), then this
% script assigns defaults and issues warnings or errors as necessary.
%
% Inputs:
%    datparams - {tdelay ratio_data datadir}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind rorfit twofit intscheme nnodes consider_error Mcell_err T0_err}
%    matparams - {Mcell aniso BI n_toplayer TCR}
%      Mcell     - {LCTE LCTEG gg}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%
%
% Outputs
%    tdelay  - COLUMN VECTOR of time delay data
%    ratio_data - COLUMN VECTOR of ratio data; can also be "dr_data",
%                 representing the ratio of the ratios of two data sets at
%                 different TDTR system conditions.
%    datadir - STRING directory path in which (processed) data, fit
%              results, and sensitivity plots are kept and saved.
%
%
%    tau_rep - laser repetition time (sec) = 1/f_rep
%    f       - modulation frequency 
%    r_pump  - pump 1/e2 radius (m)
%    r_probe - probe 1/e2 radius (m)
%
%
%    Zind    - Closest index of tdelay such that Zdelay = tdelay(Zind).
%              Zdelay is user-specified time delay (ps) for shortest time
%              delay from which to evaluate goodness-of-fit.
%    rorfit  - TRUE if fitting dr_data, ratio-of-ratio data from two TDTR
%              data sets on the same sample.
%    twofit  - TRUE if simultaneously fitting two TDTR data sets 
%              for two different modulation frequencies, same sample.
%    intscheme - Numerical integration algorithm.
%                0 = Legendre-Gauss, 1 = rombint, 2 = Simpson.
%    nnodes   - For Legendre-Gauss and Simpson integration,
%               number of nodes or integration points. Default >= 35.
%    consider_error - TRUE if FIT function is called from an errorbars
%                     script. Supports errorbar calculation functionality.
%    Mcell_err - {LCTE_err LCTEG_err gg_err}; a cell array of uncertainties
%                assigned to each element in Mcell.
%    T0_err - assigned uncertainty (absolute) in the measured temperature.
%
%
%    LCTE    - vertcat(lambda,C,t,eta) for the one-channel overlayers
%    LCTEG   - vertcat(lambda,C,t,eta,G) for N-channels of substrate
%     lambda   - VECTOR of thermal conductivities (W/mK) (layer 1 is the topmost layer)
%     C        - VECTOR of specific heats (J/m3-K)
%     t        - VECTOR of layer thicknesses (last layer is alway treated semi-inf, but
%                you still have to enter something).  (m)
%     eta      - VECTOR of anisotropic ratio (kx/ky), use ones(length(lambda)) for
%                isotropic
%     G        - conductance between N-channel substrate and overlayers
%    gg      - NxN matrix representing the heat coupling rates (W/m^3-K)
%              between the N heat channels in the substrate. Only the i < j
%              indices of gg(i,j) are used.
%    aniso   - VECTOR of booleans for each overlayer: aniso(j) is TRUE if
%              jth layer is allowed to have a variable eta between 0 and 1.
%    BI      - TRUE if bidirectional heat flow is present.
%    n_toplayer - Specifies number of layers above the plane where
%                 heat is deposited in the bidirectional heat flow model.
%    TCR     - temperature coefficient of reflectivity  (1/K)...doesn't affect ratio
%
%
%    T0      - nominal measured temperature (K)
%    T_Mcell - {T_LCTE,T_LCTEG,T_gg}
%       T_LCTE - Cell array, where T_LCTE{i,j} is a two-column array
%                of [T, V(T)] reference information to determine LCTE(i,j)
%                at the self-consistent temperature T0+dT indicated by
%                steady-state and per-pulse laser heating.
%       T_LCTEG - same, for LCTEG
%       T_gg    - same, for gg
%    A_pump  - pump intensity (W)...doesn't effect ratio
%    A_probe - probe intensity (W)...doesn't effect ratio
%    absC - [T, absC(T)] two-column array for temperature dependence of
%           the absorption coefficient of the transducer at the laser
%           wavelength at the specified temperature T.
%    perpulse - TRUE if considering per-pulse heating.
%    jabs - index j for the absorption layer in LCTE.
%    jtrans - index j for the transducer layer in LCTE. jabs and jtrans
%             are necessary for per-pulse heating calculation.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: The TDTR_FIT_TTM_vH2.m, TDTR_MANUALFIT_TTM_vH2.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 7-Apr-2014 - comments updated.
%                   30-Apr-2014 - Jun Liu adds twofit to vH1.
%                   14-Jul-2014 - vH2, comments for twofit.
%------------- BEGIN CODE --------------

%% Check input parameters, assign defaults, errors, warnings as necessary
PL = [length(datparams) length(sysparams) length(calparams) ...
      length(matparams) length(Tparams)];

% calculation parameters
%calparams = {Zind rorfit twofit intscheme nnodes consider_error Mcell_err T0_err}
if PL(3) < 1, [~,Zind] = min(abs(tdelay - 100e-12)); else Zind = calparams{1}; end
if PL(3) < 2, rorfit = 0; else rorfit = calparams{2}; end
if PL(3) < 3, twofit = 0; else twofit = calparams{3}; end
if PL(3) < 4, intscheme = 0; else intscheme = calparams{4}; end
if PL(3) < 5, nnodes = 35; else nnodes = calparams{5}; end
if PL(3) < 6, consider_error = 0; else consider_error = calparams{6}; end
if PL(3) < 7
    Mcell_err = {zeros(size(LCTE)),zeros(size(LCTEG)),zeros(size(gg))};
    LCTE_err = Mcell_err{1};
    LCTEG_err = Mcell_err{2};
    gg_err = Mcell_err{3};
else
    Mcell_err = calparams{7};
    if length(Mcell_err) < 3
        error('Insufficient elements in Mcell_err for FIT_TTM.');
    else
        LCTE_err = Mcell_err{1};
        LCTEG_err = Mcell_err{2};
        gg_err = Mcell_err{3};
    end
end
if PL(3) < 8, T0_err = 0; else T0_err = calparams{8}; end
%if PL(3) < 9, P0_err = 0; else P0_err = calparams{9}; end

% data parameters
%datparams = {tdelay ratio_data datadir}, or dr_data, or two_data
if PL(1) < 1, error('no time vector'); else tdelay = datparams{1}; end
if PL(1) < 2, error('no data vector'); 
else
    if rorfit % ratio-over-ratio fitting
        dr_data = datparams{2};
    elseif twofit % simultaneous fitting to two data sets
        two_data = datparams{2};
        ratio_data = two_data{1};
        ratio_data2 = two_data{2};
    else
        ratio_data = datparams{2}; 
    end
end
% NB: now that I've set the two ratio data sets to one cell in datparams,
% I don't need to decide whether datparams has 3 or 4 elements.
if PL(1) < 3, datadir = pwd; else datadir=datparams{3}; end

% pump-probe system parameters
%sysparams = {tau_rep f r_pump r_probe};
if PL(2) < 1, tau_rep = 1/80e6; 
    warning('Defaulting to 80 MHz repetition rate'); else tau_rep = sysparams{1}; end
if PL(2) < 2
    warning('Defaulting to 9.8 MHz modulation frequency');
    if rorfit || twofit
        f = [9.8e6 9.8e6]; % two modulation frequencies
    else % simple ratio fit, one frequency
        f = 9.8e6;
    end
else f = sysparams{2}; % 1x1 if not rorfit, 1x2 if rorfit.
end 
if PL(2) < 3, error('pump spot size not specified.'); else r_pump = sysparams{3}; end
if PL(2) < 4, r_probe = r_pump;
    warning('Defaulting to r_probe = r_pump'); else r_probe = sysparams{4}; end

% material parameters
%matparams = {Mcell aniso BI n_toplayer TCR};
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
if PL(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end

% Parameters for self-consistent temperature w/laser heating.
%Tparams = {T0, T_Mcell, A_pump, A_probe, absC, perpulse, jabs, jtrans};
if PL(5) < 1, T0 = -1;     else T0 = Tparams{1}; end
if PL(5) < 2, T_Mcell = -1; else T_Mcell = Tparams{2}; end
if PL(5) < 3, A_pump = 10e-3;  warning('Defaulting to A_pump = 10mW.'); 
   else A_pump = Tparams{3}; end
if PL(5) < 4, A_probe = 10e-3; warning('Defaulting to A_probe = 10mW.'); 
   else A_probe = Tparams{4}; end
if PL(5) < 5, absC = [295 0.13]; else absC = Tparams{5}; end
if PL(5) < 6, perpulse = 0; else perpulse = Tparams{6}; end
if PL(5) < 7, jabs = 0; warning('Defaulting to jabs = 0.'); 
    else jabs = Tparams{7}; end
if PL(5) < 8, jtrans = n_toplayer+2; warning('Defaulting to jtrans = n_toplayer+2.');
    else jtrans = Tparams{8}; end
    
%% Reconstruct all cellparams with any new default values
% if rorfit == TRUE, f holds 2 values, and is not directly readable by
% REFL, which expects a scalar.
% Similarly, if twofit == TRUE, datparams contains two ratio data arrays.
sysparams = {tau_rep f r_pump r_probe};

if rorfit
    datparams = {tdelay dr_data datadir};
elseif twofit
    two_data = {ratio_data, ratio_data2};
    datparams = {tdelay two_data datadir};
else
    datparams = {tdelay ratio_data datadir};
end

calparams = {Zind rorfit twofit intscheme nnodes consider_error Mcell_err T0_err};
matparams = {Mcell aniso BI n_toplayer TCR};
Tparams = {T0, T_Mcell, A_pump, A_probe, absC, perpulse, jabs, jtrans};