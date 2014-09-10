function [Integrand, G]=TDTR_TEMP_TTM_vH2(kvectin,freq,LCTE,LCTEG,gg,r_pump,r_probe,A_pump)
%TDTR_TEMP_TTM_vH2 - Computes frequency domain average temperature response
%to periodic gaussian pump beam, probed by another gaussian beam, for a
%sample whose substrate's heat flow follows the two-temperature model.
%
% This program is vectorized/optimized to handle all frequencies 
% simultaneously (f is a ROW vector)
%
% Syntax: [Integrand, G]=TDTR_TEMP_TTM_vH1(kvectin,freq,LCTE,LCTEG,gg,r_pump,r_probe,A_pump)
%
% Inputs:
%    kvectin - vector of wavenumber (m^-1)
%    freq    - excitation frequency (Hz), ROW vector
%    LCTE    - vertcat(lambda,C,t,eta)
%    LCTEG   - vertcat(lambda,C,t,eta,G) for N-channels of substrate
%    (lambda)  - vector of thermal conductivities, 
%                lambda(1)=top surface,(W/m-K)
%    (C)       - vector of volumetric specific heat (J/m3-K)
%    (t)       - thicknesses of each layer (layer N will NOT be used, semiinfinite)
%    (eta)     - anisotropy of each layer.
%    (G)       - conductance between N-channel substrate and overlayers
%    gg      - NxN matrix representing the heat coupling rates (W/m^3-K)
%              between the N heat channels in the substrate. Only the i < j
%              indices of gg(i,j) are used.
%    r_pump  - Pump spot size (m)
%    r_probe - Probe spot size (m)
%    A_pump  - Pump power (W) hitting the sample, used to ESTIMATE 
%              amplitude (not used for fitting)
%
% Outputs:
%    Integrand - G.*Kernal; see David Cahill's 2004 paper.
%    G         - The layer G(k), where k is spatial wavenumber
%
%
% Other m-files required: AnalyticEig_2x2xkxf.m, the mtimesx package
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_REFL_TTM_vH1.m
% TO DO: remove hardcoded adiabatic condition on second substrate channel

% Author: Rich Wilson / Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 18-Mar-2013 - Rich Wilson's 2D, two-channel substrate
%                                 (TDTR_TEMPvect_mod2.m).
%                   later, 2013 - Greg Hohensee, parallelized the (k,f)
%                                 eigenvalue problem (TDTR_TEMP_vTTM.m).
%                   27-Mar-2014 - header comments reformatted,
%                                 TDTR_TEMP_TTM_vH.m
%                   7-Apr-2014  - Introduced LCTE, LCTEG, gg data
%                                 structures, made substrate channel
%                                 anisotropies depend on eta's, not
%                                 hard-coded for spinladder. vH1.m.
%                   14-July-2014 - vH2. No change.
%------------- BEGIN CODE --------------
%% Assign variables from previous versions of this code.
% They are unpacked from my other data structures so I don't need to edit 
% most of the TEMP_TTM code.
lambda = LCTE(1,:);
C = LCTE(2,:);
t = LCTE(3,:);
eta = LCTE(4,:);

lambda2 = LCTEG(1,1);
lambda3 = LCTEG(1,2);
C2 = LCTEG(2,1);
C3 = LCTEG(2,2);
eta2 = LCTEG(4,1);
eta3 = LCTEG(4,2);
G2 = LCTEG(5,1)*1e9; % interfaces in LCTEG are given in same W/m-K units
                     % as 1nm interface layers in LCTE.
                     % So, convert to W/m^2-K.
%G3 = LCTEG(5,2)*1e9; % defaulted to 0 in TDTR_TEMP_vTTM.m

% coupling parameter is given in units of pW/nm^3-K = 1e15 W/m^3-K
g = gg(1,2)*1e15; % 2-channel assumption

%% Rich's code, plus my parallelization of the eigenvalue problem
Nfreq=length(freq);
kvect=kvectin(:)*ones(1,Nfreq);
Nlayers=length(lambda); %# of layers, not including two-channel substrate
Nint=length(kvectin); %# of different frequencies to calculate for

%k is a COLUMN vector (actually a matrix that changes down the rows)
%f is a ROW vector

% definitions
ii=sqrt(-1);
omega=2*pi*freq;
kvect2 = kvect.^2;


% Preallocate for speed
M = zeros(Nlayers,Nint,Nfreq,2,2);

for n = 1:Nlayers;
    q = sqrt(eta(n)*kvect2+ii*C(n)*ones(Nint,1)*omega/lambda(n));
    phi = q*t(n);
    phi(real(phi)>50)=50+ii*imag(phi(real(phi)>50)); % caps the real part of phi at 50?
    M(n,:,:,1,1) = cosh(phi) ;
    M(n,:,:,1,2) =  -1/lambda(n)*1./q.*sinh(phi);  
    M(n,:,:,2,1) =  -lambda(n)*q.*sinh(phi) ;
    M(n,:,:,2,2) =  cosh(phi);
end
% n = 1, calculate M(1,...)
    aa1 = M(1,:,:,1,1);
    bb1 = M(1,:,:,1,2);
    cc1 = M(1,:,:,2,1);
    dd1 = M(1,:,:,2,2);

for n = 2:1:Nlayers % Get [aa1,bb1;cc1,dd1] = M(N) * M(N-1) * ... * M(2) * M(1)
    aa2 = M(n,:,:,1,1);
    bb2 = M(n,:,:,1,2);
    cc2 = M(n,:,:,2,1);
    dd2 = M(n,:,:,2,2);
    aa_new = aa1.*aa2 + bb2.*cc1;
    bb_new = aa2.*bb1 + bb2.*dd1;
    cc_new = aa1.*cc2 + cc1.*dd2;
    dd_new = bb1.*cc2 + dd1.*dd2;
    aa1 = aa_new;
    bb1 = bb_new;
    cc1 = cc_new;
    dd1 = dd_new;
end

%% Solve for G(k) with knowledge of paramaters of 2 channel substrate
% Doubly parallelized eigenvalue problem for substrate layer n = Nlayers+1
% For spinladders, first channel is isotropic (eta2 = 1), but 2nd channel 
% is 1D, so "infinitely anisotropic" (eta3 = 0).
TMatrix(1,1,:,:) = eta2*kvect2 + (ii*(ones(Nint,1) * omega)*C2 + g)/lambda2;
TMatrix(2,2,:,:) = eta3*kvect2 + (ii*(ones(Nint,1) * omega)*C3 + g)/lambda3;
TMatrix(1,2,:,:) = -g/lambda2;
TMatrix(2,1,:,:) = -g/lambda3;
%%

[LambdaEig,Lambda2] = AnalyticEig_2x2xkxf(TMatrix);
% Lambda2 is a (2,1,k,f) vector, LambdaEig is a (2,2,k,f) vector.

% Replaces the for loop
Kchannels = repmat([lambda2 0; 0 lambda3],[1,1,Nint,Nfreq]);

% Solve for Eigenfunctions/Eigenvectors of characteristic equation
% Take the square root
Lambda = sqrt(Lambda2);
Lambda = -sign(real(Lambda)).*Lambda;

% Need to relate Surface Flux to Surface Temperature
% [q1; q2] - [X2 Y2; X3 Y3]*[T2; T3] - need to find [X2 Y2; X3 Y3]
EigMat = [Lambda(1,1,:,:).*LambdaEig(1,1,:,:),LambdaEig(1,2,:,:).*Lambda(2,1,:,:);
          Lambda(1,1,:,:).*LambdaEig(2,1,:,:),LambdaEig(2,2,:,:).*Lambda(2,1,:,:)   ];

% n-dimensional matrix multiplication code packet, from Internets.
flux_BC = -mtimesx(Kchannels,EigMat);

% Analytic inversion of many 2x2 matrices
ConstMatrix = zeros(size(flux_BC));
detFluxBC = flux_BC(1,1,:,:).*flux_BC(2,2,:,:) - flux_BC(1,2,:,:).*flux_BC(2,1,:,:);
ConstMatrix(1,1,:,:) = flux_BC(2,2,:,:) ./ detFluxBC;
ConstMatrix(1,2,:,:) = -flux_BC(1,2,:,:) ./ detFluxBC;
ConstMatrix(2,1,:,:) = -flux_BC(2,1,:,:) ./ detFluxBC;
ConstMatrix(2,2,:,:) = flux_BC(1,1,:,:) ./ detFluxBC;

XY = mtimesx(LambdaEig,ConstMatrix);
%% T(1,2,3,S), Q(1,2,3,S) matrix:
% Seven equations, eight unknowns. Define new variables T(j)' = T(j) /
% q(S), so that we now have seven equations, seven unknowns.
% The [T1' T2' T3' TS' q1' q2' q3'] matrix is as follows:
% For [T2,T3] = X * [q2, q3], [T1, q1] = M * [TS, qS]
% and q2 = G2(T1 - T2), q3 = G3(T1 - T3), and q1 = q2 + q3, we have
% 1 0 0 -A 0 0 0          = B
% 0 1 0 0 0 -X11 -X12     = 0
% 0 0 1 0 0 -X21 -X22     = 0
% 0 0 0 -C 1 0 0          = D
% 0 0 0 0 1 -1 -1         = 0
% -G2 G2 0 0 0 1 0        = 0
% -G3 0 G3 0 0 0 1        = 0

% Hardcode the spinladder situation: g3 = 0. Rich's algebra.
X2 = ones(1,Nint,Nfreq);
X2(1,:,:) = XY(1,1,:,:);
GX = G2*ones(1,Nint,Nfreq)./(1+G2*X2);
G = squeeze((dd1-GX.*bb1)./(GX.*aa1-cc1));

%% Compute kernal, convolute with G, return integrand.
% Note: the form of the kernal does matter. Rich and Joe/David have
% different methods and different kernals.

Kernal=A_pump/(2*pi)*exp(     -(r_pump^2+r_probe^2)/8*kvect2).*kvect; %The rest of the integrand
if(length(G(:,1)) ~= length(Kernal(:,1)))
    G = G';
end
Integrand = G.*Kernal; % Integrand for Legendre-Gauss integration in REFL_V4B

%dT = G*A_pump;
%plot(kvect,G.*Kernal,'x')
end
