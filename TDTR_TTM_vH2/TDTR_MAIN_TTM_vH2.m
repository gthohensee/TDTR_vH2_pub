%TDTR_MAIN_TTM_vH2 - Perform thermal modeling of TTM substrate sample.
% This script will do thermal modeling and print/save a figure representing
% the fit. It can handle manual and automatic fitting, to either the ratio
% or the "ratio of ratios" of two TDTR datasets.
% It will also make sensitivity plots and error bar calculations upon 
% request.
%
% Syntax:  run this script through an analyze_yymmdd_template_TTM script, 
% not independently, unless you define the expected variables first. I'd 
% make it a function, but it expects a lot of variables.
%
% List of variables expected to exist in the workspace:
% [Double-check against variables defined in analyze_yymmdd_template_TTM]
%    ii              - index from the analyze script's for loop.
%                      Used to pick out stack(ii,:) for individual samples.
%    r_pump, r_probe - pump and probe spot sizes, in meters
%    tau_rep         - Ti:sapphire pulse repetition period
%    f               - pump/EOM modulation frequency
%    TCR             - thermoreflectance coefficient
%    fname           - string indicating data file name.
%    datain          - string indicating path to data file.
%    datadir         - string indicating path to data folder.
%    stack           - cell string matrix from the analyze script.
%    LCTE            - 4xM thermal parameter matrix
%    LCTEG           - 5xN thermal parameters for N-channel substrate
%    gg              - NxN matrix containing coupling parameters
%    Xguess          - used in FIT and MANUALFIT functions
%    XguessIJC       - used in FIT and MANUALFIT functions
%    Xijc            - used in FIT and MANUALFIT functions
%    Zdelay          - starting time delay in ps from which to 
%                      plot data and evaluate goodness of fit.
%    tdelay_min      - minimum of delay time to model
%    tdelay_max      - maximum of delay time to model
%    P0,T0           - pressure and temperature. Set to -1 to assume ambient.
%    BI              - TRUE if bidirectional heat flow
%    n_toplayer      - thermal modeling parameter for bidirectional heat flow.
%                      indicates number of layers above the heat deposition
%                      plane in the sample.
%    intscheme       - integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
%                      2 = Simpson integration.
%    nnodes          - number of nodes for Legendre-Gauss integration;
%                      affects numerical accuracy. Don't go below 35 nodes
%                      unless you know what you're doing.
%    options         - tolerances for auto-fitting
%    senseplot       - Generate Sensitivity Plot? This option is available
%                      dynamically in MANUALFIT.
%    ebar            - TRUE if calculating Errorbars (takes longer).
%    importdata      - TRUE if fitting data. FALSE if just running 
%                      sensitivity plots or error bar calculations.
%    manualfit       - TRUE if fitting manually, FALSE if auto-fitting.
%    rorfit          - TRUE if fitting by the ratio-of-ratios signal.
%    twofit          - TRUE if fitting data at two modulation frequencies
%                      simultaneously.
%
%    g - two-channel linear coupling constant, W/m^3-K
%WARNING: g > 10e25 may result in ROUNDING ERRORS in eta-dependent terms,
%such that the sensitivity to eta does NOT track the "effective thermal model"
%prediction (fully coupled one-channel substrate).
%Make sure you fit for log(g) or divide by 1e14 or something.
%Otherwise fitting (esp. autofitting) takes ages.
%
% Important products of this script:
%    Xsol - thermal parameter fit
%    fit result figure in .fig and .eps formats, saved to datadir.
%
% Other m-files required: TDTR_REFL_TTM_vH2.m, TDTR_FIT_TTM_vH2.m,
%                         TDTR_MANUALFIT_TTM_vH2.m, TDTR_TDEP_TTM_vH2.m,
%                         SS_Heating_TTM_vH2.m, TDTR_REFL_TTM_vH2.m,
%                         TDTR_TEMP_TTM_vH2.m, INITIALIZE_CELLPARAMS_TTM_vH2.m,
%                         senseplot_TTM_vH2.m, errorbars_TTM_vH2.m,
%                         errorbar_conditions_TTM_vH2.m, genLCTEtext_TTM.m,
%                         AnalyticEig_2x2xkxf.m, extract_interior.m,
%                         rombint_VV3.m, SimpsonInt.m, lgwt_V4.m,
%                         mtimesx package.
% Subfunctions: none
% MAT-files required: none
%
% See also: The TDTR_vH2 package

% Author: Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 3/20/2013 - TDTR_MAIN_vTTM.m
%                   3/27/2014 - TDTR_MAIN_TTMsubstrate_vH.m
%                   4/7/2014  - TDTR_MAIN_TTM_vH1.m
%                   4/8/2014  - comments, harmonized with TDTR_MAIN_vH1.m
%                   7/14/2014 - vH2. Twofit functionality.
%------------- BEGIN CODE --------------
%%
% I removed keepdir functionality for this version of MAIN, because this
% MAIN assumes the data folder and fnames are all generated by the
% analyze script.

%______PROGRAM OPTIONS______________________
% These are specified in the analyze script.
%___________________________________________
%%
    
% If you want no laser heating in the model, but to still have T
% not at room temperature, give writeLCTE your specific T, and set
% T0 = -1.

Mcell = {LCTE, LCTEG, gg};
T_Mcell = {T_LCTE, T_LCTEG, T_gg};
matparams = {Mcell aniso BI n_toplayer TCR};

if rorfit || twofit, f = [f1 f2]; end
sysparams = {tau_rep f r_pump r_probe};

Tparams = {T0, T_Mcell, A_pump, A_probe, absC, perpulse, jabs, jtrans};
% calparams assignment needs to wait until Zind is defined.
% datparams assignment needs to wait until tdelay and data are defined.

% initialize variables relating to errorbar calculation
consider_error = 0;           % MAIN does not do error bars directly.
LCTE_err = zeros(size(LCTE)); % See errorbars_TTM_vH2.m for error bars.
LCTEG_err = zeros(size(LCTEG));
gg_err = zeros(size(gg));
Mcell_err = {LCTE_err, LCTEG_err, gg_err};

% define tolerances for automatic fitting in FIT / errorbar calculations
if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end

if senseplot || ebar
    % DESIGN CHOICE: I calculate sensitivities and error bars using the
    % thermal model generated best fit to the data, as specified by LCTE,
    % NOT with respect to the data itself. This ensures that these 
    % calculations do not second-guess my judgement of the fit, are not 
    % influenced by random noise, and represent sensitivities and error 
    % bars with respect ONLY to systematic or thermal parameters.
    
    %vector of time delays (used to generate sensitivity plots)
    tdelay=logspace(log10(tdelay_min),log10(tdelay_max),20)';
    [~,Zind] = min(abs(tdelay_data - Zdelay*1e-12));
    calparams = {Zind rorfit intscheme nnodes consider_error Mcell_err};
    
    % calculate data array.
    sysparams{2} = f(1); % first (or only) frequency
    [deltaR_data,ratio_data]=...
            TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);                 
    if rorfit || twofit
        sysparams{2} = f(2); % second frequency
        [deltaR_data2,ratio_data2]=...
            TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,A_pump,intscheme,nnodes);
        if rorfit
            dr_data = ratio_data ./ ratio_data2;
            datparams = {tdelay_data dr_data datadir};
        else % twofit must be true...
            two_data = {ratio_data,ratio_data2};
            deltaR_twodata = {deltaR_data,deltaR_data2};
            datparams = {tdelay_data two_data datadir};
        end
        sysparams{2} = f{1}; % return to default.                
    else % neither rorfit nor twofit
        dr_data = 0;
        datparams = {tdelay_data ratio_data datadir};
    end
end
%% -------------Generate Sensitivity Plot--------------
if senseplot
    tic
    fprintf('Calculating Sensitivities...Please Wait.\n')  
    
    % calls REFL, requres cellparams
    senseplot_TTM_vH2(datparams,sysparams, calparams, matparams, Tparams); 
    toc
    fprintf('Sensitivities calculated.\n')
end
%% --------------Compute Errorbars---------------------
if ebar
    tic
    %options = optimset('TolFun',1e-1,'TolX',1e-1); % default tolerances for errorbars
    fprintf('Calculating Errorbar...Please Wait.\n')
    
    % errorbars_TTM_vH2 calls FIT, requires all cellparams
    [kErr_perc_TTM, kErr_abs_TTM, ...
        ErrSummary_perc_TTM, ErrSummary_abs_TTM] = ...
        errorbars_TTM_vH2(Xijc,datparams,sysparams,...
                          calparams,matparams,Tparams,options); 
    
    % save Xsol and errorbar result to file
    dlmwrite(strcat(datadir,'Xsol_ErrSummary_perc_',filename,'.txt'),...
                    vertcat(Xsol,ErrSummary_perc));
    toc
    fprintf('Errorbar calculated.\n')
end

%% --------------Import Data---------------------------
% reads and extracts data matrix, generates Zind, checks that ratio is 
% positive, linearizes V(out) if desired.
if importdata
    DM1 = dlmread(datain);
    tdelay_raw  =DM1(:,2)*1e-12; %imported in picoseconds.  Change to seconds.
    Vin_raw     =DM1(:,3); 
    Vout_raw    =DM1(:,4);
    ratio_raw   =DM1(:,5);
    
    [~,ratio_data]          =extract_interior(tdelay_raw,ratio_raw, tdelay_min,tdelay_max);
    [~,Vin_data]            =extract_interior(tdelay_raw,Vin_raw,   tdelay_min,tdelay_max);
    [tdelay_data,Vout_data] =extract_interior(tdelay_raw,Vout_raw,  tdelay_min,tdelay_max);
    
    % Define Zind: Zdelay = tdelay(Zind), approximately.
    [~,Zind] = min(abs(tdelay_data - Zdelay*1e-12));
    
    if rorfit || twofit
        DM2 = dlmread(datain2); % second frequency data set
        tdelay_raw2 =DM2(:,2)*1e-12; %imported in picoseconds.  Change to seconds.
        Vin_raw2    =DM2(:,3); 
        Vout_raw2   =DM2(:,4);
        ratio_raw2  =DM2(:,5);
        [~,ratio_data2]          =extract_interior(tdelay_raw2,ratio_raw2, tdelay_min,tdelay_max);
        [~,Vin_data2]            =extract_interior(tdelay_raw2,Vin_raw2,   tdelay_min,tdelay_max);
        [tdelay_data2,Vout_data2] =extract_interior(tdelay_raw2,Vout_raw2,  tdelay_min,tdelay_max);
    end
    
    if rorfit
        % If the data points don't match up for a given pair of datasets,
        % we need to interpolate. If they do line up, ROR_interp just
        % returns dr_data, with tdelay_data unchanged.
        [tdelay_data,dr_data,tdelay_deviation] = ...
            ROR_interp(tdelay_data,ratio_data,...
                       tdelay_data2,ratio_data2);
    end
    
    % If the analysis comes with a reasonable linear fit to V_out, 
    % use that to smooth ratio_data. [NOT TESTED]
    if exist('Voutlinfit','var') && exist('fitOK','var')
        if Voutlinfit && fitOK && ~rorfit
            Vout_lin = m1 * tdelay_data + m2;
            r_lin = -Vin_data ./ Vout_lin;
            ratio_data = r_lin;
        end
    end
            
    %% Compose remaining parameter cells, 
    % now that Zind, tdelay, and data are defined.
    if rorfit
        datparams = {tdelay_data dr_data datadir}; 
        f = [f1 f2];
        sysparams{2} = f;
    elseif twofit 
        two_data = {ratio_data,ratio_data2};
        datparams = {tdelay_data two_data datadir}; 
        f = [f1 f2];
        sysparams{2} = f;
    else
        datparams = {tdelay_data ratio_data datadir}; 
    end
    
    calparams = {Zind rorfit twofit intscheme nnodes consider_error Mcell_err};

%% --------------Perform Fit--------------------------
   % Fitting assigns values to Xsol, Z, deltaR_model, ratio_model, LCTE.
   % LCTE changes if eta changes by changing L, or if the entire LCTE
   % changes from self-consistent temperature changes.
    if manualfit
        if rorfit, fprintf('Manual fitting to R(%0.2f MHz)/R(%0.2f MHz)...\n',f(1)/1e6,f(2)/1e6);
        elseif twofit, fprintf('Manual fitting to R(%0.2f MHz) and R(%0.2f MHz)...\n',f(1)/1e6,f(2)/1e6);
        else fprintf('Manual fitting to r(t)...\n'); 
        end
        
        [Xsol,Z,deltaR_model,ratio_model,Mcell,T_adj,dr_model] = ...
                TDTR_MANUALFIT_TTM_vH2(XguessIJC, datparams,...
                sysparams, calparams, matparams, Tparams);
        
        XguessIJC(:,1) = Xsol;
    else
        Xguess = XguessIJC(:,1);
        
        fprintf('Please wait for automatic fitting...\n');
        tic
        Xsol = fminsearch(@(X) TDTR_FIT_TTM_vH2(X, Xijc, datparams,...
                               sysparams, calparams, matparams, Tparams),...
                               Xguess,options);
                           
        % fminsearch doesn't output the fitted model, so get it here.                
        [Z,deltaR_model,ratio_model,Mcell,T_adj,dr_model]=...
            TDTR_FIT_TTM_vH2(Xsol,Xijc,datparams,sysparams,calparams,matparams,Tparams);
        toc
    end
    fprintf('Data fit completed\n');
else % not importing data, just doing errorbars, sensitivities, or making figures.
    Xsol = Xguess;
    Z = 0;
    
    % Assign tdelay, deltaR, ratio models to the model reference "data"
    % from the sensitivity and error bar calculations.
    tdelay_data = tdelay;
    
    if twofit
        ratio_model = two_data;
        deltaR_model = deltaR_twodata;
    else
        ratio_model = ratio_data;
        deltaR_model = deltaR_data;
    end
    T_adj = 0;
    dr_model = dr_data;
    
end
% From here on, ratio_model and deltaR_model are 2-element cells
% if twofit was true. "datparams{2}" contains the 2-element cell of ratio data.
% deltaR_data1 and deltaR_data2 remain separate and not used in twofit.

%% --------------Generate save data --------------------
% uses saveres
% uses fname for results
% uses ratio_data
% saves a sensitivity plot, if you asked for it in the analyze script.

saveres=input('Want to save results? (0=no, 1=yes): ');
%saveres = 1;  % TRUE to save fit result figure and/or sensitivity figure.

if saveres
    figsens = 202;
    figfit = 203;
    
    % Save the workspace
    save(strcat(datadir,'Results_', fname(1:end-4),'.mat'))
    
    % define the data and model arrays for the results figure
    if rorfit
        plot_data = dr_data;
        plot_model = dr_model;
        ytext = sprintf('R(%0.1f MHz)/R(%0.1f MHz)\n',f(1)/1e6,f(2)/1e6);
        fittext = 'FIT_TTMs_ROR_';
    elseif twofit
        plot_data = ratio_data{1};
        plot_model = ratio_model{1};
        plot_data2 = ratio_data{2};
        plot_model2 = ratio_model{2};
        ytext = sprintf('R(%0.1f MHz)and R(%0.1f MHz)\n',f(1)/1e6,f(2)/1e6);
        fittext = 'FIT_TTMs_2Freq_';    
    else
        plot_data = ratio_data;
        plot_model = ratio_model;
        ytext = '-V(in)/V(out)';
        fittext = 'FIT_TTMs_';
    end
    
    % Create a figure for the last fit
    figure(figfit)
    clf;
    
    if importdata
        hd = loglog(tdelay_data,plot_data,'ko');
        hold on;
        hm = loglog(tdelay_data,plot_model,'k-');
        
        if twofit
            hd2 = loglog(tdelay_data,plot_data2,'bs');
            hm2 = loglog(tdelay_data,plot_model2,'b-');
            set(hd2,'LineWidth',2)
            set(hm2,'LineWidth',2)
        end
        hold off;
    else % if importdata is false, then we didn't fit anything,
         % so there's no data points to plot -- just the model.
        hm = loglog(tdelay_data,plot_model,'k-');
        if twofit
            hm2 = loglog(tdelay_data,plot_model2,'b-');
            set(hm2,'LineWidth',2)
        end
    end
    set(hm,'LineWidth',2)
    fontsize = 16;
    condtext = sprintf('T = %0.1f K, Z = %0.2d, t(fit) = %i ps',T_adj,Z,Zdelay);
    xlabel('time delay (ps)','FontSize',fontsize)
    ylabel(ytext, 'FontSize',fontsize)
    title(condtext,'FontSize',fontsize)
    axis([tdelay_min 5e-9 min(0.1,min([plot_data; plot_data2])) ...
                          2*max([plot_data;plot_data2])])
    set(gca, 'XTick', [1e-11 2e-11 5e-11 1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
    set(gca, 'XTickLabel', [10 20 50 100, 200, 500, 1000, 2000, 5000, 1e4]);
    set(gca, 'YTick', [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 30 50]);
    set(gca, 'XMinorTick', 'off');
    set(gca, 'YMinorTick', 'off');
    set(gca, 'TickLength' , [.02 .02]);
    set(gca,'FontSize',fontsize);
    %% Provide a summary of final fit parameters in the figure.
   
    % Warning: uses "ii" index from the analyze script's for loop!
    XsolIJC = XguessIJC;
    Mcell = {LCTE, LCTEG, gg};
    partext = genLCTEtext_TTM(Mcell,XsolIJC,stack(ii,:),substack(ii,:));
    
    % write LCTEstr contents to a text box in the figure
    pBox = annotation('textbox',[0.15,0.2,0.8,0.3]);
    set(pBox,'Units','characters')
    set(pBox,'HorizontalAlignment','left')
    set(pBox,'FontSize',12)
    set(pBox,'String',partext);
    set(pBox,'LineStyle','none');
    
    % also record the data file used
    dBox = annotation('textbox',[0.15,0.85,0.8,0.07]);
    set(dBox,'Units','characters')
    set(dBox,'Interpreter','none')
    set(dBox,'HorizontalAlignment','left')
    set(dBox,'FontSize',12)
    set(dBox,'String',{'Data file:'; fname(1:end-4)});
    set(dBox,'LineStyle','none');
    
    %% Save the figure to .fig and .eps files
    saveas(figfit, strcat(datadir,fittext,fname(1:end-4),'.fig'))
    print(figfit,'-depsc',strcat(datadir,fittext,fname(1:end-4),'.eps'))
    
    if senseplot
        tag = input('Label the sensitivity plot: ','s');
        save(strcat(datadir,'/SENS_', tag, '.mat'))
        figure(figsens)
        saveas(figsens, strcat(datadir,'/SENS_',tag,'.fig'))
        print(figsens,'-depsc',strcat(datadir,'/SENS_',tag,'.eps'))
    end
end
%----------------------------------------------------
fprintf('Program Completed\n')
beep
pause(0.1);
beep
%---------------- END OF CODE -----------------------