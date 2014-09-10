function [deltaR,ratio]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,...
                                         intscheme,nnodes)
%TDTR_REFL_TTM_vH2 - Calculates the Reflectivity Signal and Ratio for
%samples using a two-temperature model for the substrate.
%
%This function can handle uni-directional or bi-directional heat flow, as
%well as simple anisotropy (in-plane vs. cross-plane). It can perform
%integration of the TDTR_TEMP integrands with either Legendre-Gauss,
%rombint, or Simpson integration, as desired.
%
%In order to speed up the code, it is parallelized...the convention is...
%  tdelay is COLUMN vector of desired delay times
%  Mvect (the fourier components) are ROW vectors
%  Matrices have size, length(tdelay) x length(Mvect)
%
% Syntax:  [deltaR,ratio]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,...
%                                         intscheme,nnodes)
%
% Inputs (scalars unless told otherwise):  
%    tdelay  - COLUMN VECTOR of time delay data
%    sysparams - {tau_rep f r_pump r_probe}
%    matparams - {Mcell aniso BI n_toplayer TCR}; Mcell = {LCTE LCTEG gg};
%    intscheme - Numerical integration algorithm.
%                0 = Legendre-Gauss, 1 = rombint, 2 = Simpson.
%    nnodes   - For Legendre-Gauss and Simpson integration,
%               number of nodes or integration points. Default >= 35.
%See INITIALIZE_CELLPARAMS_TTM_vH2.m for details on sysparams, matparams.
%
% Outputs:
%    deltaR - Complex number VECTOR. Real part is model V(in), imaginary
%             part is model V(out).
%    ratio  - model -V(in)/V(out).
%
% Other m-files required: TDTR_TEMP_TTM_vH2.m, TDTR_TEMP_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: The TDTR_V4 package, TDTR_REFL_vH2.m

% Author: Gregory Hohensee / built from Joseph Feser's TDTR_REFL_V4.
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 12-Sep-2012 - J. Feser's REFL_V4 published
%                   2013 - added bidirectional heat flow, REFL_V4B.
%                   25-Mar-2014 - standardized comments into TDTR_vH package.
%                   26-Mar-2014 - use LCTE instead of lambda,C,t,eta. 
%                                 Allow various integration schemes.
%                   27-Mar-2014 - TDTR_REFL_TTM_vH.m, permits modeling
%                                 of two-temperature substrates.
%                   7-Apr-2014  - vH1.
%                   14-July-2014  - vH2. No change from vH1.
%------------- BEGIN CODE --------------
%% check input
PL = [0 length(sysparams) 0 length(matparams) 0];

if nargin < 6, nnodes = 35; end
if nargin < 5, intscheme = 0; end
if nargin < 4, A_pump = 10e-3; end % arbitrary positive value
if nargin < 3, error('Insufficient parameters for TDTR_REFL_TTM_vH2.'); end

% pump-probe system parameters
%sysparams = {tau_rep f r_pump r_probe};
if PL(2) < 1, tau_rep = 1/80e6; 
    warning('Defaulting to 80 MHz rep. rate'); else tau_rep = sysparams{1}; end
if PL(2) < 2, f = 9.8e6;
    warning('Defaulting to 9.8 MHz modulation'); else f = sysparams{2}; end
if PL(2) < 3, error('pump spot size not specified.'); else r_pump = sysparams{3}; end
if PL(2) < 4, r_probe = r_pump;
    warning('Defaulting to r_probe = r_pump'); else r_probe = sysparams{4}; end

% material parameters
%matparams = {Mcell aniso BI n_toplayer TCR}; Mcell = {LCTE, LCTEG, gg};
if PL(4) < 1
    error('Mcell not specified in REFL_TTM.'); 
else
    Mcell = matparams{1};
    if length(Mcell) < 3
        error('Insufficient elements in Mcell for REFL_TTM.');
    else
        LCTE = Mcell{1};
        LCTEG = Mcell{2};
        gg = Mcell{3};
    end
end
%if pl(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end

%% reminders for users of previous versions of TTM thermal modeling
% C2 = LCTEG(2,1);
% C3 = LCTEG(2,2);
% lambda2 = LCTEG(1,1);
% lambda3 = LCTEG(1,2);
% G2 = LCTEG(5,1);
% G3 = LCTEG(5,2);
% g = gg(1,2);
% TTM = true;
%% FORTRAN thermal model - there is no Fortran code for TTM.
%% MATLAB port, bi-directional anisotropic TTM substrate thermal model.
% That is, I see no reason why this can't also handle anisotropic cases.

% constants
cap1 = 10;   % for fmax, M
cap2 = 4*pi; % for kmax
% History lesson: % % % % % % %
% RW uses cap1 = 15 in his REFL/TEMP scripts with Simpson
% integration. More accuracy in integration at the cost of speed.
%
% cap2 = 4*pi gives the kmax that Rich uses in TDTR_TEMPvect_mod2.m,
% and I've carried over his TEMP code for TEMP_TTM_vH.
% For TDTR_TEMP_vH, I keep the same kmax as Joe had, and as Rich had
% for Simpson integration of regular samples.
%
% RW uses cap2 = 2 instead of 1.5 for his REFL/TEMP scripts,
% which use Simpson integration exclusively. 2x catches higher k-modes,
% is more accurate for fixed node density over k-space.
% For TTM, k(RW) != k(JPF), where k(JPF) is in use in TEMP_V4.
% % % % % % % % % % % % % % % %

fmax=cap1/min(abs(tdelay)); %maximum frequency considered (see Cahill 2004 RSI paper)
ii=sqrt(-1);

M=cap1*ceil(tau_rep/min(abs(tdelay))); %Highest Fourier component considered
mvect=-M:M; %Range of Fourier components to consider (Vector)
fudge1=exp(-pi*((mvect/tau_rep+f)/fmax).^2);%artificial decay (see RSI paper)
fudge2=exp(-pi*((mvect/tau_rep-f)/fmax).^2);

kmin = 0;
kmax = cap2/sqrt(r_pump^2+r_probe^2);

dT1=zeros(1,length(mvect))';
dT2=zeros(1,length(mvect))';

if BI
    dT1u = zeros(1,length(mvect))';
    dT2u = zeros(1,length(mvect))';

    dT1d = zeros(1,length(mvect))';
    dT2d = zeros(1,length(mvect))';
    
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
            [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...
            % calculate the heat flow up, layer n_toplayer up to layer 1
            I1u = TDTR_TEMP_vH2(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump);
            I2u = TDTR_TEMP_vH2(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump);
            dT1u = weights'*I1u;
            dT2u = weights'*I2u;
            % calculate the heat flow down, layer n_toplayer+1 down to bottom
            I1d = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep+f,LCTEd,LCTEG,gg,r_pump,r_probe,A_pump);
            I2d = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep-f,LCTEd,LCTEG,gg,r_pump,r_probe,A_pump);
            dT1d = weights'*I1d;
            dT2d = weights'*I2d;
        case 1 % rombint
            % The Rombint solution doesn't converge quickly for TTM.
            dT1u=rombint_VV3(@(kvect) TDTR_TEMP_vH2(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump),0,kmax,length(mvect));
            dT2u=rombint_VV3(@(kvect) TDTR_TEMP_vH2(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump),0,kmax,length(mvect));
            dT1d=rombint_VV3(@(kvect) TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep+f,LCTEd,LCTEG,gg,r_pump,r_probe,A_pump),0,kmax,length(mvect));
            dT2d=rombint_VV3(@(kvect) TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep+f,LCTEd,LCTEG,gg,r_pump,r_probe,A_pump),0,kmax,length(mvect));
        case 2 % Simpson integration
            kvect = linspace(0,kmax,nnodes)';
            I1u = TDTR_TEMP_vH2(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump);
            I2u = TDTR_TEMP_vH2(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump);
            I1d = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep+f,LCTEd,LCTEG,gg,r_pump,r_probe,A_pump);
            I2d = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep-f,LCTEd,LCTEG,gg,r_pump,r_probe,A_pump);
            dT1u=SimpsonInt(kvect,I1u);
            dT2u=SimpsonInt(kvect,I2u);
            dT1d=SimpsonInt(kvect,I1d);
            dT2d=SimpsonInt(kvect,I2d);
        otherwise
            error('integration scheme not properly specified. See analyze template.');
    end
    % make the parallel sum of the temperatures for up and down
    dT1 = dT1u .* dT1d ./ (dT1u + dT1d);
    dT2 = dT2u .* dT2d ./ (dT2u + dT2d);
else % BI = 0
    switch intscheme
        case 0 % Legendre-Gauss
            [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...

            I1 = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep+f,LCTE,LCTEG,gg,r_pump,r_probe,A_pump);
            I2 = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep-f,LCTE,LCTEG,gg,r_pump,r_probe,A_pump);
            dT1 = weights'*I1;
            dT2 = weights'*I2;
        case 1 % rombint
            % the TTM TDTR calculation requires the solution to 2x2
            % eigenvalue problems across (k,f) space. The Rombint solution
            % doesn't converge quickly.
            dT1=rombint_VV3(@(kvect) TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep+f,LCTE,LCTEG,gg,r_pump,r_probe,A_pump),kmin,kmax,length(mvect));
            dT2=rombint_VV3(@(kvect) TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep-f,LCTE,LCTEG,gg,r_pump,r_probe,A_pump),kmin,kmax,length(mvect));
        case 2 % Simpson
            kvect = linspace(kmin,kmax,nnodes)';
            I1 = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep+f,LCTE,LCTEG,gg,r_pump,r_probe,A_pump);
            I2 = TDTR_TEMP_TTM_vH2(kvect,mvect/tau_rep-f,LCTE,LCTEG,gg,r_pump,r_probe,A_pump);
            dT1=SimpsonInt(kvect,I1);
            dT2=SimpsonInt(kvect,I2);
        otherwise
            error('integration scheme not properly specified. See analyze template.');
    end
end


dT1 = dT1';
dT2 = dT2';

expterm=exp(ii*2*pi/tau_rep*(tdelay*mvect));
%Retemp=(ones(length(tdelay),1)*(dT1.*fudge1+dT2.*fudge2)).*expterm;
Retemp=(ones(length(tdelay),1)*(dT1'.*fudge1+dT2'.*fudge2)).*expterm;
%Imtemp=-ii*(ones(length(tdelay),1)*(dT1-dT2)).*expterm;
Imtemp=-ii*(ones(length(tdelay),1)*(dT1-dT2)').*expterm;
Resum=sum(Retemp,2); %Sum over all Fourier series components
Imsum=sum(Imtemp,2);

deltaRm=TCR*(Resum+ii*Imsum); %
deltaR=deltaRm.*exp(ii*2*pi*f*tdelay); %Reflectance Fluxation (Complex)

ratio=-real(deltaR)./imag(deltaR);
end