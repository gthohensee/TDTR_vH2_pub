function [partext] = genLCTEtext_TTM(Mcell,XsolIJC,stack,substack)
%genLCTEtext_TTM - Formats LCTE solution with (stack) layer IDs into a
%text block with a format similar to the standard .par file for the
%original FORTRAN TDTR thermal models. Include the two-temperature
%substrate model layer as well.
%
% Syntax:  [partext] = genLCTEtext_TTM(Mcell,XsolIJC,stack,substack)
%                                    
% Inputs:
%    Mcell    - {LCTE,LCTEG,gg}. Cell array of matrices for thermal parameters
%               in two-temperature model. LCTEG is the TTM substrate.
%               gg are the channel coupling parameters.
%    XsolIJC  - [V, i, j, c], where V = LCTE(i,j) for c = 1, LCTEG(i,j) 
%               for c = 2, and gg(i,j) for c = 3.
%    stack    - 1xN string array of model layer identifiers. See
%               analyze_yymmdd_template.m and writeLCT.m.
%    substack - 1xN cell array of strings: channel identifiers. Use labels 
%               like "phonon", "magnon", "electron", "momeraths".
%
% Outputs:
%    partext -  Array of text formatted for printing onto a MATLAB figure
%               text box. Intended for a fit results figure, showing the
%               data, the fitted curve, and the thermal parameters.
%
% Example: 
%    ...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_MAIN_vH.m, analyze_yymmdd_template.m, writeLCT.m

% Author: Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Apr 2014; Last revision: 3-Apr-2014
%------------- BEGIN CODE --------------
% Deconstruct input
LCTEsol = Mcell{1};
LCTEGsol = Mcell{2};
gg = Mcell{3}; % matrix of g(i-j) coupling. i < j. 2x2 in TTM.
Xijc = XsolIJC(:,2:4); % Xijc indexes X across (i,j) of the Mcell{c}.
nf = length(Xijc(:,1)); % total number of fit parameters
nc = length(gg(1,:));

XsolIJ = [];     % indices of XsolIJ will be [LCTE(i,j),i,j]
XsolIJ_ttm = []; % indices of XsolIJ_ttm will be[LCTEG(i,j),i,j]
gsolIJ = [];     % indices of gsolIJ will be [gg(ij), i, j].

for i = 1:nf
    switch Xijc(i,3); % c index
        case 1
            XsolIJ = vertcat(XsolIJ,XsolIJC(i,1:3));
        case 2
            XsolIJ_ttm = vertcat(XsolIJ,XsolIJC(i,1:3));
        case 3
            gsolIJ = vertcat(XsolIJ,XsolIJC(i,1:3));
    end
end

%% Iterate through LCTE 1-channel fit parameters
fitstr = [];
if isempty(XsolIJ), nrows = 0; else nrows = length(XsolIJ(:,1)); end
for j = 1:nrows % loop through fit parameters
    if nrows == 0, break; end
    if j > 1, fitstr = strcat(fitstr, ' and'); end % good grammar!

    LCTEsol(XsolIJ(j,2),XsolIJ(j,3)) = XsolIJ(j,1);

    switch XsolIJ(j,2) % switch over [L,C,T,E] rows
        case 1 % thermal conductivity row of LCT
            if ~isempty(strfind(stack{XsolIJ(jj,3)},'/')) % interface
                fitstr1 = ' G(';
                units = ' MW/m^2-K';
                unitfactor = 1e3;
            else
                fitstr1 = ' L(';
                units = ' W/cm-K';
                unitfactor = 1e-2;
            end
        case 2 % heat capacity row of LCT
            fitstr1 = ' C(';
            units = ' J/m^3K';
            unitfactor = 1e-6;
        case 3 % thickness row of LCT
            fitstr1 = ' t(';
            units = ' nm';
            unitfactor = 1e9;
        case 4
            fitstr1 = ' eta(';
            units = 'L_x/L_z';
            unitfactor = 1;
    end
    fitstr2 = stack{XsolIJ(j,3)}; % get layer ID
    
    fitstr = strcat(fitstr, fitstr1, fitstr2, ')=', num2str(XsolIJ(j,1)*unitfactor,3),units);
end

if isempty(fitstr)
    Xsol_text = [];
else
    Xsol_text = strcat('Fit to', fitstr);
end
    

%% Switch to standard .par format from Fortran TDTR. 
% For visibility, use SI-cm units for L and C, and nm units for T.
LCTEvis = LCTEsol; 
LCTEvis(1,:) = LCTEvis(1,:) * 1e-2; % to W/cm K
LCTEvis(2,:) = LCTEvis(2,:) * 1e-6; % from J/m^3K to J/cm^3K
LCTEvis(3,:) = LCTEvis(3,:) * 1e9;  % from m to nm
LCTEvis(4,:) = LCTEvis(4,:); % dimensionless

% Create a cell of strings from LCTvis
K = length(LCTEvis(1,:));
LCTEtext = cell(1+K+1,1);
% standard format puts the layers into rows, not columns.
LCTEtext{1} = '[Lx/Lz, W/cm-K, J/cm^3-K, nm]...layerID';
for k = 1:K % iterating through the 1-channel layers
    temp = strcat('[',mat2str(LCTEvis(4,k),3), ', ', ...
                      mat2str(LCTEvis(1,k),3), ', ',...
                      mat2str(LCTEvis(2,k),3), ', ',...
                      mat2str(LCTEvis(3,k),3), ']...', stack{k});
    LCTEtext{1+k} = temp;
end
LCTEtext{1+K+1} = Xsol_text;

%% Two-temperature substrate [G2, L2, C2, G3, L3, C3]
% How should I index XsolIJ_ttm? One row per parameter.
% Each row has [V,i,j], where LCTEG(i,j) = V. But what about "g"?
% Or the possibly (N-1)! g's for the N-temp extension? Well, "i"
% refers to [L,C,T,E,G], ok. "j" refers to jth channel. OK.
% I will have a separate gsol for the (N-1)! coupling parameters.
fitstr_ttm = [];
if isempty(XsolIJ_ttm), nrows = 0; else nrows = length(XsolIJ_ttm(:,1)); end
for j = 1:nrows % loop through fit parameters
    if nrows == 0, break; end
    if j > 1, fitstr_ttm = strcat(fitstr_ttm, ' and'); end % good grammar!

    LCTEGsol(XsolIJ_ttm(j,2),XsolIJ_ttm(j,3)) = XsolIJ_ttm(j,1);

    switch XsolIJ_ttm(j,2)
        case 1 % thermal conductivity row of LCTEG
            fitstr1 = ' L(';
            units = ' W/cm-K';
            unitfactor = 1e-2;
        case 2 % heat capacity row of LCTEG
            fitstr1 = ' C(';
            units = ' J/m^3K';
            unitfactor = 1e-6;
        case 3 % thickness row of LCTEG
            fitstr1 = ' t(';
            units = ' nm';
            unitfactor = 1e9;
        case 4
            fitstr1 = ' eta(';
            units = 'L_x/L_z';
            unitfactor = 1;
        case 5
            fitstr1 = ' G(';
            units = ' MW/m^2-K';
            unitfactor = 1e3; % G in LCTEG has same units as a G-layer in LCTE.
    end
    fitstr2 = substack{XsolIJ_ttm(j,3)}; % get channel ID
    
    fitstr_ttm = strcat(fitstr_ttm, fitstr1, fitstr2, ')=', ...
                        num2str(XsolIJ_ttm(j,1)*unitfactor,3),units);
end % if XsolIJ_ttm is empty, fitstr_ttm remains empty.

%% coupling g
if isempty(fitstr_ttm) && isempty(gsolIJ)
    Xsolttm_text = [];
else
    if isempty(gsolIJ)
        fitg = []; 
    else % the next line assumes two-temperature model
        unitfactor = 1; % gg is both stored and presented in pW/nm^3-K.
        fitg = strcat(' g=', mat2str(gsolIJ(1,1)*unitfactor,3),' pW/nm^3-K and');
    end
    
    if isempty(fitstr_ttm)
        Xsolttm_text = strcat('Fit to', fitg);
    else
        Xsolttm_text = strcat('Fit to', fitg, ' and', fitstr_ttm);
    end
end
       
%% Standard .par format text, except just for the N-channel substrate
LCTEGvis = LCTEGsol;
LCTEGvis(1,:) = LCTEGvis(1,:) * 1e-2; % to W/cm K
LCTEGvis(2,:) = LCTEGvis(2,:) * 1e-6; % from J/m^3K to J/cm^3K
LCTEGvis(3,:) = LCTEGvis(3,:) * 1e9;  % from m to nm
LCTEGvis(4,:) = LCTEGvis(4,:);        % dimensionless anisotropy
LCTEGvis(5,:) = LCTEGvis(5,:) * 1e6;  % from W/m^2-K to MW/m^2-K


LCTEGtext = cell(2+nc+1,1);
LCTEGtext{1} = 'Two-channel substrate parameters:';
LCTEGtext{2} = '[MW/m^2-K, Lx/Lz, W/cm-K, J/cm^3-K]...channelID';
for n = 1:nc
    temp = strcat('[',mat2str(LCTEGvis(5,n),3), ', ',...
                      mat2str(LCTEGvis(4,n),3), ', ', ...
                      mat2str(LCTEGvis(1,n),3), ', ', ...
                      mat2str(LCTEGvis(2,n),3), ']...', substack{1});
    LCTEGtext{2+n} = temp;
end
         
LCTEGtext{2+nc+1} = Xsolttm_text;

partext = vertcat(LCTEtext, LCTEGtext);
end

