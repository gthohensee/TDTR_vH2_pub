function [LCTEG,gg,T_LCTEG,T_gg] = writeLCTEGg_vH2(subID, T0, Xijc, refdir,hsubs)
%writeT_LCTEGg - output thermal parameters in the two-temperature model
% for a given substrate at T = T0. LCTEG just represents this one layer;
% LCTE still contains the parameters for the transducer and other
% one-channel layers.
%
% Syntax:  [LCTEG,gg,T_LCTEG,T_gg] = writeLCTEGg_vH2(subID, T0, Xijc, refdir,hsubs)
%
% Inputs:
%    subID - string identifier for n-channel substrate
%    T0    - Nominal (measured) temperature; no laser heating yet.
%    Xijc  - (i,j) indices of fit parameters across c: LCTE, LCTEG, gg.
%    refdir - directory path for T-dependent material reference data
%    hsubs - substrate thickness; check TEMP_TTM to confirm whether it
%            accounts for the possibility of a thermally finite substrate!!
%
% Outputs:
%    LCTEG   - 5xN thermal matrix of N-channel substrate, analogous to
%              LCTE with extra row for conductance G (SI units) into
%              the upper one-channel layers.
%    gg      - NxN matrix of coupling parameters. Only the top half of the
%              matrix contains unique coupling parameters. gg(i,j), i < j.
%    T_LCTEG - Cell matrix of two-column matrices for each element of
%              LCTEG, indicating temperature dependence [T V(T)] of each
%              element of LCTEG.
%    T_gg    - Cell matrix of two-column matrices for each element of gg,
%              indicating temperature dependence [T g(T)] of each
%              element of gg.
%
% Example:
%     ---
%
% Other m-files required: TDTR_TDEP_general.m
% Subfunctions: none
% MAT-files required: none
% DAT-files required: several, depending on your material.
%
% See also: writeLCTE_vH2.m.

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% April 2014; Last revision: 3-April-2014
%                            14-July-2014: vH2
%------------- BEGIN CODE --------------
%% Go through stack, assign thermal parameters.
RT = 295; % default value for room temperature (K)
clear lower; % in case the workspace contains a "lower" variable.

%% check inputs
if nargin < 2, T0 = -1; end % default no temperature dependence
if nargin < 3, Xijc = 0; end % default no fit parameter
if nargin < 4
    refdir = '/Users/gregoryhohensee/Documents/MATLAB/projectbacon/exportbacon/TTM/refdir'; 
end
if nargin < 5, hsubs = 500e3; % some thermally infinite length

switch lower(subID) % lowercase
    case 'sp95' % Ca9La5Cu24O41 spin ladder
        
        n = 2; % 2-channel model for spinladder - phonons, magnons
        LCTEG = zeros(5,n);
        gg = zeros(n,n); % 2x2 matrix for two-channel model; one unique g(1,2).
        T_LCTEG = cell(size(LCTEG));
        T_gg = cell(size(gg));

        % TDTR_TDEP_general expects SI or cm-SI units, [T,V(T)] format.
        TC2_file = 'interp_C_Ca9La5Cu24O41_lattice.txt';
        TL2_file = 'K_SSa_Hess2001_0p94scale_10xinterp.txt';
        TC3_file = 'CmagINSdispersion_magnonstau.txt'; 
        TL3_file = 'K_SSmagnon_NaruseHess.txt'; 
        
        if T0 ~= -1
            [Cp_T0 Lp_T0 TCp TLp] = TDTR_TDEP_general(TC2_file, TL2_file, T0, refdir);
            [Cm_T0 Lm_T0 TCm TLm] = TDTR_TDEP_general(TC3_file, TL3_file, T0, refdir);
            
            T_LCTEG{1,1} = TLp;
            T_LCTEG{2,1} = TCp;
            T_LCTEG{1,2} = TLm;
            T_LCTEG{2,2} = TCm;
        else
            [Cp_T0 Lp_T0] = TDTR_TDEP_general(TC2_file, TL2_file, RT, refdir);
            [Cm_T0 Lm_T0] = TDTR_TDEP_general(TC3_file, TL3_file, RT, refdir);
        end
        Gp = 0.1; % 100 MW/m^2-K, same units as for LCTE interfaces.
        etap = 1; % isotropic phonons
        LCTEGp = [Lp_T0,Cp_T0,hsubs,etap,Gp]';

        Gm = 0; % adiabatic
        etam = 0; % cross-plane magnons
        LCTEGm = [Lm_T0,Cm_T0,hsubs,etam,Gm]';

        LCTEG = horzcat(LCTEGp, LCTEGm);
        
        gg(1,2) = 1; % pW/nm^3-K units are used throughout TTM_vH2 package.
        
    case 'srcuo' % Sr14Cu24O41
        
        n = 2; % 2-channel model for spinladder - phonons, magnons
        LCTEG = zeros(5,n);
        gg = zeros(n,n); % 2x2 matrix for two-channel model; one unique g(1,2).
        T_LCTEG = cell(size(LCTEG));
        T_gg = cell(size(gg));

        % TDTR_TDEP_general expects SI or cm-SI units, [T,V(T)] format.
        TC2_file = 'interp_C_Ca9La5Cu24O41_lattice.txt';
        TL2_file = 'Interp_scale_ 1.00_ka_Hess data.txt';
        TC3_file = 'CmagINSdispersion_magnonstau.txt'; 
        TL3_file = 'Interp_scale_ 0.70_Kcmagnon.txt'; 
        TGp_file = 'Interp_shifted_ 0.000000_G_conductance_Greg.txt';
        
        pathG = strcat(refdir, '\', TGp_file);
        TGp = dlmread(pathG);
        if T0 ~= -1
            [Cp_T0 Lp_T0 TCp TLp] = TDTR_TDEP_general(TC2_file, TL2_file, T0, refdir);
            [Cm_T0 Lm_T0 TCm TLm] = TDTR_TDEP_general(TC3_file, TL3_file, T0, refdir);
            % define TGp = [T, G(T)], where G = 0.1 ==> 100 MW
            T_LCTEG{1,1} = TLp;
            T_LCTEG{2,1} = TCp;
            T_LCTEG{1,2} = TLm;
            T_LCTEG{2,2} = TCm;
            T_LCTEG{5,1} = TGp;
        else
            [Cp_T0 Lp_T0] = TDTR_TDEP_general(TC2_file, TL2_file, RT, refdir);
            [Cm_T0 Lm_T0] = TDTR_TDEP_general(TC3_file, TL3_file, RT, refdir);
        end
       [~,iG] = min(abs(TGp(:,1) - T0));
       Gp = TGp(iG,2);
        
        % Gp = 0.1; % 100 MW/m^2-K, same units as for LCTE interfaces.
        etap = 1; % isotropic phonons
        LCTEGp = [Lp_T0,Cp_T0,hsubs,etap,Gp]';

        Gm = 0; % adiabatic
        etam = 0; % cross-plane magnons
        LCTEGm = [Lm_T0,Cm_T0,hsubs,etam,Gm]';

        LCTEG = horzcat(LCTEGp, LCTEGm);
        
        gg(1,2) = 1; % pW/nm^3-K units are used throughout TTM_vH2 package.
    otherwise
        fprintf('ERROR: unable to parse stack of layer IDs. \nPlease check analyze_yymmdd.m, writeLCT.m\n')
end

% Handles T-independent and fit parameters in T_LCTEG
for i = 1:5
    for j = 1:length(LCTEG(1,:))
        % Assign null values to fit parameters in T_LCTEG, to prevent TDEP from
        % "updating" a fitting parameter when handling self-consistent heating.
        if Xijc ~= 0
            for x = 1:length(Xijc(:,1))
                if Xijc(x,1) == i && Xijc(x,2) == j && Xijc(x,3) == 2
                    T_LCTEG{i,j} = [300 -1];
                end
            end
        end
        
        % Assign T-independent default values to all remaining voids in
        % T_LCTEG
        if isempty(T_LCTEG{i,j})
            T_LCTEG{i,j} = [300 LCTEG(i,j)];
        end
    end
end

% Handles T-independent and fit parameters in T_gg
n = length(gg(:,1));
for i = 1:n
    for j = i+1:n
        % Assign null values to fit parameters in T_gg, to prevent TDEP from
        % "updating" a fitting parameter when handling self-consistent heating.
        if Xijc ~= 0
            for x = 1:length(Xijc(:,1))
                if Xijc(x,1) == i && Xijc(x,2) == j && Xijc(x,3) == 3
                    T_gg{i,j} = [300 -1];
                end
            end
        end
        
        % Assign T-independent default values to all remaining voids in
        % T_gg
        if isempty(T_gg{i,j})
            T_gg{i,j} = [RT gg(i,j)];
        end
    end
end

end
%------------- END CODE --------------