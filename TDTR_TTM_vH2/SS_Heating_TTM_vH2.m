function dTss = SS_Heating_TTM_vH2(matparams,r_pump,r_probe,absC,A_tot,intscheme,nnodes)
%SS_Heating_TTM_vH2 - Calculates the steady state heating of an anisotropic 
%multilayer stack. Can handle bidirectional heat flow.
%
% Syntax:  dTss = SS_Heating_TTM_vH2(matparams,r_pump,r_probe,absC,A_tot,intscheme,nnodes)
%
% Inputs:
%    matparams - {Mcell aniso BI n_toplayer TCR}
%    Mcell   - {LCTE LCTEG gg}
%    LCTE    - vertcat(lambda,C,t,eta);
%       lambda  - VECTOR of thermal conductivities (W/mK) (layer 1 is the topmost layer)
%       C       - VECTOR of specific heats (J/m3-K)
%       t       - VECTOR of layer thicknesses (last layer is alway treated semi-inf, but
%                 you still have to enter something).  (m)
%       eta     - VECTOR of anisotropic ratio (kx/ky), use ones(length(lambda)) for
%                 isotropic
%[See INITIALIZE_CELLPARAMS_vH2.m for more details.]
%    r_pump  - pump 1/e2 radius (m)
%    r_probe - probe 1/e2 radius (m)
%    A_tot   - Total laser power (mW) hitting the sample inside of a duty
%              cycle of the ~200Hz probe modulation. This is AFTER
%              correcting the powermeter reading for power meter
%              calibration (1.08x) and the objective's transmission
%              coefficient. See the analyze template for details.
%    absC    - SCALAR fraction of laser power absorbed by the top layer
%              [Poor design; absC in higher level functions refers
%               to T_absC, not a scalar absC as it does here]
%
% Outputs:
%    dTss - steady state temperature rise due to laser heating.
%
% Example: 
%    --
%
% Other m-files required: TDTR_TEMP_TTM_vH2.m, TDTR_TEMP_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: analyze_yymmdd_template_TTM.m

% Author: Gregory Hohensee / built from Joseph Feser's SS_Heating.m
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 12-Sep-2012 - J. Feser's TDTR_V4 published
%                   26-Mar-2014 - SS_Heating_vH.m
%                   7-Apr-2014 - vH1.m
%                   8-Apr-2014 - harmonized with non-TTM version
%                   14-July-2014 - vH2. No change from vH1.
%------------- BEGIN CODE --------------

%% unpack matparams
PL = [0 0 0 length(matparams) 0];

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

%% assign remaining variables
dTss = 0;

%%
%cap1 = 10; % for fmax, M; irrelevant for SS heating at f = 0 Hz.
cap2 = 4*pi; % for kmax
f=0; %laser Modulation frequency, Hz
A_abs = absC*A_tot; %absorbed power
r = sqrt(r_pump^2+r_probe^2);

kmin=1/(10000*r); %smallest wavevector (can't use exactly zero, 
                  %because of a numerical issue -> 0/0 = NaN
                  %when calculating SS heating with f = 0.
                  
kmax = cap2/sqrt(r_pump^2+r_probe^2); % cap2 is RW convention for TTM
% RW uses cap2 = 2 instead of 1.5 for his REFL/TEMP scripts,
% which use Simpson integration exclusively. 2x catches higher k-modes,
% is more accurate for fixed node density over k-space.
% For TTM, k(RW) = 4*pi^2*k(JPF)? Where k(JPF) is in use in TEMP_V4.

switch intscheme
    case 0 % Legendre-Gauss
        [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...
    case 1 % rombint
        % do nothing %
    case 2 % Simpson (uses nnodes)
        % do nothing %
    otherwise
        error('integration scheme not properly specified. See analyze template.');
end

if BI
    totlayers = length(LCTE(1,:));
    % Split the material parameters into "up" and "down" sections,
    % relative to the heater and thermometer at interface of n_toplayer
    % and n_toplayer+1.
    for j = n_toplayer:-1:1
        LCTEu(:,n_toplayer+1-j) = LCTE(:,j);
    end
    LCTEd = LCTE(:,n_toplayer+1:totlayers);
    % ASSUMPTION: 
    % The N-channel substrate is assumed to be "below" everything else.
    % That's why the "up" half is modeled in regular TEMP_vH, whereas
    % the "down" half is handled in TEMP_TTM_vH.
    % TTM has NOT been tested with bidirectional heat flow! -3/27/2014.
    
    switch intscheme
        case 0 % Legendre-Gauss
            % calculate the heat flow up, layer n_toplayer up to layer 1
            Iu = TDTR_TEMP_vH2(kvect,f,LCTEu,r_pump,r_probe,A_abs);
            dTu = weights'*Iu;
            
            % calculate the heat flow down, layer n_toplayer+1 down to bottom
            Id = TDTR_TEMP_TTM_vH2(kvect,f,LCTEd,LCTEG,gg,r_pump,r_probe,A_abs);
            dTd = weights'*Id;
        case 1 % rombint
            dTu=rombint_VV3(@(kvect) TDTR_TEMP_vH2(kvect,f,LCTEu,r_pump,r_probe,A_abs),kmin,kmax,1);
            dTd=rombint_VV3(@(kvect) TDTR_TEMP_TTM_vH2(kvect,f,LCTEd,LCTEG,gg,r_pump,r_probe,A_abs),kmin,kmax,1);
        case 2 % Simpson integration
            kvect = linspace(kmin,kmax,nnodes)';
            Iu = TDTR_TEMP_vH2(kvect,f,LCTEu,r_pump,r_probe,A_abs);
            Id = TDTR_TEMP_TTM_vH2(kvect,f,LCTEd,LCTEG,gg,r_pump,r_probe,A_abs);
            dTu=SimpsonInt(kvect,Iu);
            dTd=SimpsonInt(kvect,Id);
    end
    % make the parallel sum of the temperatures for up and down
    dTss = dTu .* dTd ./ (dTu + dTd);
else % BI = 0    
    switch intscheme
        case 0 % Legendre-Gauss
            I = TDTR_TEMP_TTM_vH2(kvect,f,LCTE,LCTEG,gg,r_pump,r_probe,A_abs);
            dTss = weights'*I;
        case 1 % rombint
            dTss=rombint_VV3(@(kvect) TDTR_TEMP_TTM_vH2(kvect,f,LCTE,LCTEG,gg,r_pump,r_probe,A_abs),kmin,kmax,1);
        case 2 % Simpson
            kvect = linspace(kmin,kmax,nnodes)';
            I = TDTR_TEMP_TTM_vH2(kvect,f,LCTE,LCTEG,gg,r_pump,r_probe,A_abs);
            dTss=SimpsonInt(kvect,I);
    end
end
end
%----------------- END CODE --------------------