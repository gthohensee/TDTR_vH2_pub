function [Xsol,Z,deltaR_model,ratio_model,LCTE,T_adj,N]=...
    TDTR_MANUALFIT_vH2(XguessIJ,datparams,...
                       sysparams, calparams, matparams, Tparams)
%TDTR_MANUALFIT_vH2 - Manually fit thermal model to TDTR ratio data.
%This program lets you iteratively fit the thermal model to the ratio data
%by varying the LCT thermal parameters specified in XguessIJ. It will tell
%you the residual deviation (crude goodness-of-fit) Z, allow you to
%generate a sensitivity plot at any time, and do thermal modeling with
%anisotropic unidirectional or bidirectional heat flow. It also handles
%temperature dependence in a self-consistent manner with steady-state and
%per-pulse heating through calls to TDTR_TDEP_vH2.m.
%
%It will NOT allow you to:
% ** reduce the transducer thickness below the absorption length if you
%    model an absorption layer, 
% ** compute the goodness-of-fit weighted by the sensitivities. 
%    It'd take too long to generate sensitivities every time you change 
%    the fit.
%
% Syntax:  [Xsol,Z,deltaR_model,ratio_model,LCTE,T_adj]=...
%    TDTR_MANUALFIT_vH2(XguessIJ,datparams,...
%                       sysparams, calparams, matparams, Tparams)
%
% Inputs:
%    XguessIJ  - Mx3 matrix: each row represents a fit parameter Xguess,
%                so [Xguess i j]. (i,j) designates Xguess => LCTE(i,j).
%    datparams - {tdelay ratio_data datadir}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind sigfit intscheme nnodes consider_error LCTE_err T0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_vH2.m for details on the params inputs.]
%
% Outputs:
%    Xsol        - Final values for the fitted parameters.
%    Z           - Goodness-of-fit. In Joe's words: "Typically, this is 
%                  the sum of the squares of the residuals, but you might 
%                  want to weight these by the sensitivity, particularly 
%                  if you don't intend to calculate the errorbars!"
%    deltaR      - Complex number array. Real part is the model V(in),
%                  imaginary part is the model V(out). Represents change
%                  in reflectance from pump heating.
%    ratio_model - Ratio signal -V(in)/V(out) from the thermal model.
%    LCTE        - The values in Xsol may be tied to the anisotropy eta,
%                  or they could affect T_adj, which affects values in
%                  LCTE generally. So, output a new LCTE.
%    T_adj       - T_adj = T0 + dTss + dTpp, the "actual" temperature
%                  in Kelvin adjusted for steady-state (and per-pulse
%                  if perpulse is TRUE) heating.
%
% Example:
%    --
%
% Other m-files required: TDTR_REFL_vH2.m, TDTR_TEMP_vH2.m, 
%                         TDTR_TDEP_vH2.m, SS_Heating_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_FIT_TTM_vH2.m, TDTR_MANUALFIT_TTM_vH2.m, TDTR_FIT_vH2.m

% Author: Gregory Hohensee
% Acknowledgement: built from TDTR_FIT_V4B, my bi-directional tweak to the
% original TDTR_FIT_V4.m, written by the great and powerful Joseph P.
% Feser.
% University of Illinois at Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history:  25-Mar-2014 - written
%                    8-Apr-2014 - more comments, harmonized with TTM
%                    14-July-2014 - vH2. No changes.
%------------- BEGIN CODE --------------
%% Check input parameters, assign defaults, errors, warnings as necessary
INITIALIZE_CELLPARAMS_vH2;

%% Assign fit variables X according to their XguessIJ index.
% XguessIJ = [X,i,j] indexes X across LCTE(i,j)
Xij = XguessIJ(:,2:3); % for easier comparison to TDTR_FIT_T_vH.m
nrows = length(Xij(:,1));

% iterate over fit parameters
for x = 1:nrows, LCTE(Xij(x,1),Xij(x,2)) = XguessIJ(x,1); end 
%% User input loop
done = 0; N = 1; % N is for adjusting the normalization factor when fitting by V(in) or V(out)
while done ~= 1
    %% Self-consistent steady-state (and optionally per-pulse) heating 
    if T0 ~= -1
        [dTss, dTpp, LCTE] = TDTR_TDEP_vH2(matparams,sysparams,...
                                               Tparams,intscheme,nnodes);
        matparams{1} = LCTE; % update                                   
        fprintf('T0 = %0.2f K, dTss = %0.2f, dTpp = %0.2f\n',T0,dTss,dTpp)
    end
    
    %% Compute model
    [deltaR_model,ratio_model]=TDTR_REFL_vH2(tdelay,matparams,sysparams,...
                                            A_pump,intscheme,nnodes);
    switch sigfit
        case 1 % V(in) fit
            % Construct normalized V(in) model and data,
            % relative to its value at Zdelay picoseconds.
            Vin_model = real(deltaR_model);
            Vin_model_Zdelay = Vin_model(Zind) / N;
            nVin_model = Vin_model / Vin_model_Zdelay;
    
            Vin_data_Zdelay = Vin_data(Zind);
            nVin_data = Vin_data / Vin_data_Zdelay;
        case 2 % V(out) fit
            % Construct normalized V(out) model and data,
            % relative to its mean value near Zdelay picoseconds.
            Vout_model = imag(deltaR_model);
            Vout_model_Zdelay = Vout_model(Zind) / N;
            nVout_model = Vout_model / Vout_model_Zdelay;

            Vout_data_Zdelay = Vout_data(Zind);
            nVout_data = Vout_data / Vout_data_Zdelay;
        otherwise % ratio fit
            % do nothing here
    end
    
    %% Update the data and fit comparison figure
    figure(10)
     
    switch sigfit
        case 1
            loglog(tdelay,nVin_data,'ob',tdelay,nVin_model,'r');
            
            ylabel('normalized V(in)','FontSize',16);
            axis([1e-10 10e-9 min([nVin_data;nVin_model])/2 2])
        case 2
            semilogx(tdelay,nVout_data,'ob',tdelay,nVout_model,'r');
        
            ylabel('normalized V(out)','FontSize',16);
            axis([1e-10 10e-9 min([0;nVout_data;nVout_model]) 2])
            set(gca,'YMinorTick','on')
            set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]);
            
            
        otherwise
            loglog(tdelay,ratio_data,'ob',tdelay,ratio_model,'r');
            ylabel('Ratio','FontSize',16); 
            axis([100e-12 10e-9 min([ratio_data;ratio_model])/2 max([ratio_data;ratio_model])*2])
            
    end
    if sigfit ~= 2
        set(gca, 'YTick', [0.5, 1, 2, 5, 10, 20, 50, 100]);
        set(gca, 'YTickLabel', [0.5, 1, 2, 5, 10, 20, 50, 100]);
        set(gca,'YMinorTick','off')
    end
    xlabel('Time delay (ps)','FontSize',16);
    set(gca, 'XTick', [1e-11, 2e-11, 5e-11, 1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
    set(gca, 'XTickLabel', [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4]);
    set(gca, 'TickLength' , [.02 .02]);
    set(gca,'XMinorTick','off')
    set(gca,'FontSize',16);
    
    % goodness-of-fit Z from fractional residuals
    switch sigfit
        % For V(in) and V(out) fits, it may be wiser to weight the 
        % residuals by the amount of time delay away from Zdelay, 
        % where V(in) is pinned to unity.
        case 1
            res=(1-(nVin_model(Zind:length(nVin_model)) ...
                    ./nVin_data(Zind:length(nVin_data))) ).^2;
        case 2
            res=(1-(nVout_model(Zind:length(nVout_model)) ...
                    ./nVout_data(Zind:length(nVout_data))) ).^2;
        otherwise % ratiofit
        res=(1-(ratio_model(Zind:length(ratio_model)) ...
                ./ratio_data(Zind:length(ratio_data))) ).^2;
    end 
    Z = sum(res)
    
    %% Inform user of current fit parameters
    
    % Define symbols and units for all possible fit parameters
    tag = {'L','C','t','eta'};
    units = {'W/m-K', 'J/cm^3-K', 'nm', '(Lx/Lz)'};
    scale = [1 1e-6 1e9 1];
    
    % Report current fit parameters
    Mat = sprintf('Current material fit parameters LCTE(i,j):');
    for x = 1:nrows % iterate through fit parameters
        itag = Xij(x,1);
        Mat = char(Mat,sprintf('%s(%i) = %0.4f %s',...
                   tag{itag},Xij(x,2),scale(itag)*LCTE(Xij(x,1),Xij(x,2)),units{itag}));
    end
    Mat
    
    %% get and clean input
    done = input('Enter 1 if done, 2 for a sensitivity plot, 3 to rescale; else hit "Enter": ');
    if isempty(done) 
        done = 0;
    else
        if sum(done == [0 1 2 3]) == 0
            done = 0; fprintf('Hey! Invalid input. Go home, you are drunk.\n')
        end
    end
    %% execute user choice
    switch done
        case 3 % tweak the normalization factor
            if sigfit ~= 1 && sigfit ~= 2
                fprintf('You cannot normalize the ratio!\n')
            else
                fprintf('Pick the normalization point on the plot...\n')
                [~,N] = ginput(1);
                fprintf('You picked %0.3f. Calculating new fit...\n',N);
            end
        case 2 % sense plot
            figure(202)
            clf
            
            senseplot_vH2(datparams,sysparams, calparams, matparams, Tparams);
            
            savesens = input('Enter 1 to save sensitivity plot to datadir: ');
            if savesens == 1
                tag = input('Name the sensitivity plot: ','s');
                save(strcat(datadir,'/SENS_', tag, '.mat'))
                fignum=202;
                figure(fignum)
                saveas(fignum, strcat(datadir,'/SENS_',tag,'.fig'))
                print(fignum,'-depsc',strcat(datadir,'/SENS_',tag,'.eps'))
            end
        case 0 % manual fitting
            for x = 1:nrows % iterate through fit parameters
                itag = Xij(x,1);
                temp = input(sprintf('Adjust parameter %s(%i): ',...
                                        tag{itag},Xij(x,2)));
                if ~isempty(temp)
                    XguessIJ(x,1) = temp / scale(itag);
                    
                    % if adjusting conductivity, then eta may
                    % change, depending on aniso.
                    if Xij(x,1) == 1 && aniso(Xij(x,2))
                        LCTE(4,Xij(x,2)) = LCTE(4,Xij(x,2)) * LCTE(1,Xij(x,2)) / temp;
                    end

                    % having used the old LCTE(1,j) to update eta if necessary,
                    % can now update LCTE(i,j).
                    LCTE(Xij(x,1),Xij(x,2)) = temp / scale(itag);
                end
                
            end % end iteration through fit parameters
    end % end switch for execution of user's decision
    
    matparams{1} = LCTE; % update matparams.
end % end user input loop

% Assign final temperature adjusted for SS and PP heating.
if T0 ~= -1, T_adj = T0 + dTss + dTpp; else T_adj = T0; end

% Assign Xsol values
Xsol = XguessIJ(:,1);
end % end program
