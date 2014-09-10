% Author: Gregory Hohensee, 4/7/2014
% Usage: called within errorbars_TTM_vH2 to initialize various (cell)
% matrices for error bar calculation. Specifies which parameters to
% consider in the uncertainty calculations.
%
% INPUTS: LCTE,LCTEG,Xijc
% OUTPUTS: Mcell_consider, r_probe_consider, r_pump_consider, phase_consider;
%          Mcell_err, Mcell_err_temp, r_err, degphase;
% INITIALIZES: Mcell_abs_err, r_probe_err, r_pump_err, phase_err

% Revision history: 14-July-2014: vH2. No change.
%%
% Booleans: TRUE if considering uncertainty from these parameters.
LCTE_consider = ones(4,nl);
LCTEG_consider = ones(5,nc);
gg_consider = ones(nc,nc);
r_probe_consider=1;
r_pump_consider=1;
phase_consider=1;

% initialize dimensions
nl = length(LCTE(1,:)); % number of layers
nf = length(Xijc(:,1)); % number of fit parameters
nc = length(LCTEG(1,:)); % number of channels
ng = 0; for i = 1:nc, ng = ng + (i-1); end % number of unique couplings (1 for TTM)

% Initialize fractional uncertainty in Xguess due to each parameter.
LCTE_err     = zeros(4,nl);
LCTE_err_temp = zeros(4,nl);
LCTE_abs_err = cell(4,nl);

LCTEG_err     = zeros(5,nc);
LCTEG_err_temp = zeros(5,nc);
LCTEG_abs_err = cell(5,nc);

gg_err     = zeros(nc,nc);
gg_err_temp = zeros(nc,nc);
gg_abs_err = cell(nc,nc);

r_probe_err = zeros(1,nf);
r_pump_err  = zeros(1,nf);
phase_err   = zeros(1,nf);


%% define parameters NOT to consider in error analysis (saves time)
% LCTE %
LCTE_consider(3,nl)=0; %last layer is semi-infinite
for i = 1:nl
    if LCTE(2,i) == 1e5 && LCTE(3,i) == 1e-9 % reliable signs of an interface layer
        LCTE_consider(2,i) = 0; %thermal interface layer has no capacitance
        LCTE_consider(3,i) = 0; %thermal interface conductance at fixed t/vary lambda
    end
    
    if ~aniso(i), LCTE_consider(4,i) = 0; end % skip eta for isotropic layers
end

for i = 1:nf
    switch Xijc(i,3) % c
        case 1, LCTE_consider(Xijc(i,1),Xijc(i,2)) = 0;  %solving for these!
        case 2, LCTEG_consider(Xijc(i,1),Xijc(i,2)) = 0; %solving for these!
        case 3, gg_consider(Xijc(i,1),Xijc(i,2)) = 0;    %solving for these!
    end
end

% LCTEG, gg %
LCTEG_consider(3,:) = 0; % skip semi-infinite substrate
for n = 1:nc, 
    % skip 1D and isotropic channels
    if LCTEG(4,n) == 0 || LCTEG(4,n) == 1, LCTEG_consider(4,n) = 0; end
    
    % skip adiabatic interfaces
    if LCTEG(5,n) == 0, LCTEG_consider(5,n) = 0; end 
    gg_consider(i,1:i) = 0;  % skip non-unique couplings
end

%% define percent uncertainty in each layer/parameter.
% LCTE %
% Example LCTE_err for a [*, Al, Al/substrate, substrate]
% model of Al/substrate sample, fitting for G and L(substrate),
% assuming an anisotropic substrate with in-plane conductivity
% uncertain to 5 percent.
% LCTE_err_example = [0.05 0.05 0   0;
%                     0.03 0.03 0   0.03;
%                     0.04 0.04 0   0;
%                     0    0    0   0.05];

% Implement default, generic uncertainties for the remaining LCTE...
for i = 1:nl
    if LCTE_consider(1,i) ~= 0, LCTE_err(1,i) = 0.05; end
    if LCTE_consider(2,i) ~= 0, LCTE_err(2,i) = 0.03; end
    if LCTE_consider(3,i) ~= 0, LCTE_err(3,i) = 0.04; end
    if LCTE_consider(4,i) ~= 0 && aniso(i), LCTE_err(4,i) = 0.05; end
end

% LCTEG %
for i = 1:nc
    if LCTEG_consider(1,i) ~= 0, LCTEG_err(1,i) = 0.05; end
    if LCTEG_consider(2,i) ~= 0, LCTEG_err(2,i) = 0.03; end
    if LCTEG_consider(3,i) ~= 0, LCTEG_err(3,i) = 0.04; end
    if LCTEG_consider(4,i) ~= 0, LCTEG_err(4,i) = 0.05; end
    if LCTEG_consider(5,i) ~= 0, LCTEG_err(5,i) = 0.05; end
end

% gg %
for i = 1:nc
    for j = i+1:nc
        if gg_consider(i,j) ~= 0, gg_err(i,j) = 0.10; end
    end
end

r_err=0.05; % 5% uncertainty in pump/probe beam spot sizes
degphase=0.15;  %phase uncertainty in degrees.
                %You can estimate this by starting with delphase
                %from SetPhase_vH.m, and considering the effect
                %of having many data points to average over.
                %(uncertainty is less than RMS noise)
