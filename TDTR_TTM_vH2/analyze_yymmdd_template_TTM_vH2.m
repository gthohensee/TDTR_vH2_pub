%analyze_yymmdd_templateTTM - Template for analyzing TDTR data of TTM
%                             substrate samples.
% After preparing the data, use this template to efficiently handle
% thermal model fits to a series of TDTR data files.
%
% Basic use of this template:
%    See analyze_yymmdd_template.m
%
% Other m-files required: the TDTR_TTM_vH package.
% Subfunctions: none
% MAT-files required: none
%
% See also: process_yymmdd_template.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 2-April-2014 - templateTTM, nonspecific
%                   7-April-2014 - template_TTM, based on analyze_110920,
%                                  tested prior to reorganizing senseplot.
%                   14-July-2014 - vH2, including twofit.
%------------- BEGIN CODE --------------
% Re-analysis of 110920 sp10rp low-T data to test TDTR_TTM_vH1 package. %
%% Define the directory and filenames
directory = '/Users/gregoryhohensee/Documents/MATLAB/projectbacon/exportbacon/TTM'; % pwd is the current directory
yymmdd = '110920'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('/',yymmdd,'_edit/'); % folder in directory where processed data is kept
savefolder = strcat('/',yymmdd,'_edit/'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);

%tagline = 's10rp_201'; % a shared prefix for your series of TDTR data files.
%% Repeated taglines because I need matched data at two different frequencies for rorfit
tagline = 's10rp_98'; 
% get filenames in alphanumerical order
filenames = dir(strcat(datadir,tagline,'*'));
filenames.name;

% separate low-T from high-T filenames (wouldn't need to do this if I had
% uniquely prefixed my TDTR data files)
qq = 1; jj = 1;
for i = 1:length(filenames)
    if strfind(filenames(i).name,'C') > 0 % filename refers to high temperature data
        filenamesC1{jj} = filenames(i).name;
        jj = jj+1;
    else % low-T data
        filenamesK1{qq} = filenames(i).name;
        qq = qq+1;
    end
end
clear qq; clear jj;
filenamesK1(:)

tagline = 's10rp_1608'; 
% get filenames in alphanumerical order
filenames = dir(strcat(datadir,tagline,'*'));
filenames.name;

% separate low-T from high-T filenames (wouldn't need to do this if I had
% uniquely prefixed my TDTR data files)
qq = 1; jj = 1;
for i = 1:length(filenames)
    if strfind(filenames(i).name,'C') > 0 % filename refers to high temperature data
        filenamesC2{jj} = filenames(i).name;
        jj = jj+1;
    else % low-T data
        filenamesK2{qq} = filenames(i).name;
        qq = qq+1;
    end
end
clear qq; clear jj;
filenamesK2(:)
        
nfiles = length(filenamesK1); % number of (paired) [low-T] data files
%% Type and conditions of modeling
P0 = -1;        % GPa initial. P = -1 assumes P = 0, ignores P-dependence.
BI = 0;         % boolean: TRUE if using bidirectional heat flow model
n_toplayer = 0; % # of model layers above the point where heat is deposited.

%% Common sysparams - TDTR system settings
f=9.8e6; % pump laser modulation frequency, Hz
tau_rep=1/80e6; % laser repetition period, s.
                % The repetition frequency is 80 MHz for TDTR-1, 
                % approximately 74.8 MHz for TDTR-2.
TCR=1e-4; % default coefficient of thermal reflectance
pm = 1.08; % TDTR-1: optical power is 1.08x the reading of the model 835 power meter.
cryofactor = 1 - 0.16; % factor of 16% to account for 4 cryostat window interfaces.

%% calparams - calculation parameters
Zdelay = 100; % ps starting time. In automatic fitting, goodness-of-fit
             % is calculated between Zdelay and tdelay_max. In V(in)
             % fitting, Zdelay indicates the time at which to
             % normalize the V(in) signal.

rorfit = 1; % TRUE if fitting the ratio1/ratio2 of two TDTR data sets
            % taken at two different frequencies.

twofit = 0; % TRUE if fitting ratio1, ratio2 simultaneously.
            % (only applies to two modulation frequencies; must edit
            % the rest of vH2 package to do this for two spot sizes,
            % for example.)

intscheme = 0; % integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
               % 2 = Simpson integration.
nnodes = 35;  % number of nodes for Legendre-Gauss or Simpson integration;
              % affects numerical accuracy. Don't go below 35 nodes
              % unless you know what you're doing.

%% Tparams - parameters for temperature-dependent samples
T0 = 111.3; % nominal (K). Set T0 = -1 to ignore laser heating and assume room temperature.

% optional parameters:
absC = [295 0.13]; % Transducer absorption coefficient (0.13 for room temperature Al)
perpulse = 0; % TRUE if accounting for per-pulse heating as well as steady-state

%% other operational parameters for MAIN
% Time delay boundaries for thermal modeling AND data presentation.
tdelay_min= 100e-12; % warning: the smaller this is, the longer the
                    % computation time of the thermal model.
tdelay_max= 4000e-9; % approx. max extent of our delay stage.

options = optimset('TolFun',1e-2,'TolX',1e-2); % tolerances for autofitting
                                               % and for error bars.
            
senseplot = 0; %Generate Sensitivity Plot? This option is available
               %dynamically in MANUALFIT.
              
ebar = 0; %Calculate Errorbars? 0 for no, 1 for yes (takes longer)

importdata = 1; % TRUE if fitting data. FALSE if just running sensitivity
                % plots or error bar calculations.
manualfit = 1;  % TRUE if fitting manually, FALSE if auto-fitting.

%% Thickness, resistivity, offset information
hAl = 96.1; % (nm)

rhoAl = 1.37e-8 + 2.678e-8; % Al electrical resistivity; unimportant if not using Al.
                            % 3.90e-8 ==> 186 W/m-K, typical of sputter 1.
 
%Vout_offset = 0.5; % uV, or uV/mW.

%% Arrays for dataset-dependent information

thickness = ones(nfiles,4);
thickness(:,1) = 1;    % absorption layer
thickness(:,2) = hAl; % transducer
thickness(:,3) = 1;
thickness(:,4) = 500e3; % spinladder

stack = {}; % initialize
stack(true(1,nfiles),1) = {'*'}; % first model layer
stack(true(1,nfiles),2) = {'Al'}; % second model layer
substack(true(1,nfiles),1) = {'sp95'}; % TTM substrate

aniso = [0 0]; % Al is isotropic. 
% TDTR_TEMP_TTM_vH.m is currently hardcoded for the spinladder dimensionality. %

% A 1D array of data sets requires 1D arrays of measurement parameters.
f_list = 1.608 * ones(1,nfiles); % pump modulation frequency
r_list = 10   * ones(1,nfiles); % objective lens ID: 5x, 10x, 20x, ...
jabs_list   = 1 * ones(1,nfiles); % column index of absorption layer
jtrans_list = 2 * ones(1,nfiles); % column index of transducer layer
measured_pump = 10e-3 * ones(1,nfiles);  % raw powermeter reading, W
measured_probe = 5e-3 * ones(1,nfiles); % raw powermeter reading, W

% measured temperatures (indices are matched for rorfit)
TT1608 = [111.3 111.3 131.3 131.3 154.9 154.9 182.8 182.8 215.6 215.6 254.3 254.3 300 300 80 94.4 94.4];
TT98   = [111.3 111.3 131.3 131.3 154.9 154.9 182.8 182.8 215.6 215.6 254.3 254.3 300 300 80 94.4 94.4];

% Xijc(1,1:3) = [i j c]. c = 1,2,3 indicates LCTE,LCTEG,gg.
Xijc(1,1:3) = [1,2,3]; % coupling g
Xijc(2,1:3) = [5,1,2]; % phonon conductance G
% Xijc(3,1:3) = [1,2,2]; % magnon K
% Xijc(4,1:3) = [1,1,2]; % phonon K
%% Initialize output and prev_output; careful not to erase your fit results!
output = zeros(nfiles,5);
prev_output = dlmread(strcat(datadir, yymmdd, '_solutions_',tagline,'_man.txt'));

%% Perform fitting.
nf = length(Xijc(:,1)); % number of fit parameters
Xguess = zeros(1,nf);
XguessIJC = zeros(nf,3);

for ii = 1:nfiles
   fname = filenamesK1{ii};
   fname2 = filenamesK2{ii};
   datain = strcat(datadir, fname);
   datain2 = strcat(datadir,fname2);
   sprintf('Datafile1: %s\n\n',fname)
   
   if exist('f_list','var') % overrides earlier "f" assignment
       f = f_list(ii) * 1e6;     %laser Modulation frequency, Hz
   end
   
   jabs = jabs_list(ii);     % index of absorption layer
   jtrans = jtrans_list(ii); % index of transducer layer
   
   f1 = 9.8e6; f2 = 1.608e6;
   
   % Spot sizes for TDTR-1:
   % According to Xiaojia, as of 1/23/2014:
   % 11.3 um for 5X, 6.1 um for 10X, Land 2.7 um for 20X
   % Spot sizes for TDTR-2:
   % As of March 24 2014, re: Rich Wilson:
   % 13.3 um, 6.7 um, 3.4 um, and 1.4 um  for 5x, 10x, 20x, and 50x.
   switch r_list(ii) % set pump 1/e^2 intensity focused radius, in meters.
        %case 2, r_pump = 2*11.3e-6; tc = 0.87; % "tc": transmission coefficient for this objective
        %case 5, r_pump = 11.3e-6; tc = 0.90;
        case 10, r_pump = 6.5e-6; tc = 0.80; % 6.5um was the 10x spot size for my 110920 data set.
        %case 20, r_pump = 2.7e-6; tc = 0.70;
        %case 50, r_pump = ?e-6; tc = 0.??;
        %case 100, r_pump = ?e-6; tc = 0.??;
   end 
   r_probe = r_pump;         % probe 1/e^2 radius, meters.
   
   % Laser powers that reach a typical in-air sample.
   A_pump = measured_pump(ii)*pm*tc*cryofactor;     
   A_probe = measured_probe(ii)*2*pm*tc*cryofactor; % 2x for optical chopper

   %Construct the thermal material properties matrix
   T0 = TT98(ii);
   refdir = '/Users/gregoryhohensee/Documents/MATLAB/projectbacon/exportbacon/TTM/refdir'; 
   [LCTEG,gg,T_LCTEG,T_gg] = writeLCTEGg_vH1(substack{ii,:},T0,Xijc,refdir);
   [LCTE,T_LCTE] = writeLCTE_vH1(stack(ii,:),thickness(ii,:),rhoAl,T0,Xijc,refdir);
   Mcell = {LCTE,LCTEG,gg};
   T_Mcell = {T_LCTE,T_LCTEG,T_gg};
   
   fprintf('Temperature is %0.1f K\n', T0);
   fprintf('transducer thickness is %f nm\n', thickness(ii,2));
   
   %Which variable(s) are you fitting for?
   for m = 1:nf
       switch Xijc(m,3) % c
           case 1, XguessIJC(m,1:4) = [LCTE(Xijc(m,1),Xijc(m,2)), Xijc(m,:)];
           case 2, XguessIJC(m,1:4) = [LCTEG(Xijc(m,1),Xijc(m,2)), Xijc(m,:)];
           case 3, XguessIJC(m,1:4) = [gg(Xijc(m,1),Xijc(m,2)), Xijc(m,:)];
       end
   end
   XguessIJC(1,1) = 0.54;
   XguessIJC(2,1) = .0654;
   Xguess = XguessIJC(:,1); % just the fit parameters
   
   % Set initial guess
   %XguessIJC(1:2,1) = prev_output(ii,5:6);
   %Xguess(1:2) = prev_output(ii,5:6);
   %gg(1,2) = Xguess(1,1);
   %LCTEG(5,1) = Xguess(2,1);
   manualfit = 0;
   TDTR_MAIN_TTM_vH2;
   
   % record sample and fit parameters as desired.
   output(ii,1:length(Xsol)+4) = [f(1)/1e6,T0,T_adj,Xsol',Z];%,kErr_abs];
end
%% Write fit results to a text file in the processed data folder.
caltag = '';
if rorfit, caltag = strcat(caltag,'_rorfit'); end
if manualfit, caltag = strcat(caltag,'_man');
else          caltag = strcat(caltag,'_auto'); end

dlmwrite(strcat(datadir, yymmdd,'_solutions_',tagline,caltag,'.txt'),output);

%------------- END CODE --------------