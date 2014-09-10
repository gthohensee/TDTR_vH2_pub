function [Xsol,Z,deltaR_model,ratio_model,Mcell,T_adj,dr_model]=...
         TDTR_MANUALFIT_TTM_vH2(XguessIJC,datparams,...
                               sysparams, calparams, matparams, Tparams)
%TDTR_MANUALFIT_TTM_vH2 - Manually fit thermal model to TDTR ratio data.
%This program lets you iteratively fit the thermal model to the ratio data
%by varying the LCT thermal parameters specified in XguessIJC. It will tell
%you the residual deviation (goodness-of-fit) Z, allow you to generate
%a sensitivity plot at any time, and do thermal modeling with anisotropic
%unidirectional or bidirectional heat flow. It also handles temperature
%dependence in a self-consistent manner with steady-state and per-pulse
%heating through calls to TDTR_TDEP_TTM_vH2.m.
%
%It will NOT allow you to:
% ** reduce the transducer thickness below the absorption length if you
%    model an absorption layer, 
% ** compute the goodness-of-fit weighted by the sensitivities. 
%    It'd take too long to generate sensitivities every time you change 
%    the fit.
%
% Syntax:  [Xsol,Z,deltaR_model,ratio_model,Mcell,T_adj,dr_model]=...
%         TDTR_MANUALFIT_TTM_vH2(XguessIJC,datparams,...
%                               sysparams, calparams, matparams, Tparams)
%
% Inputs:
%    XguessIJC - Mx4 matrix: each row represents a fit parameter Xguess,
%                so [Xguess i j c]. "c" designates either LCTE, LCTEG,
%                or gg, and (i,j) the index. So for c = 2, Xguess
%                represents LCTEG(i,j).
%    datparams - {tdelay ratio_data datadir}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind rorfit intscheme nnodes consider_error Mcell_err T0_err}
%    matparams - {Mcell aniso BI n_toplayer TCR}
%      Mcell     - {LCTE LCTEG gg}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_TTM_vH2.m for details on the params inputs.]
%
% Outputs:
%    Xsol        - Final values for the fitted parameters.
%    Z           - Goodness-of-fit. In Joe's words: "Typically, this is 
%                  the sum of the squares of the residuals, but you might 
%                  want to weight these by the sensitivity, particularly 
%                  if you don't intend to calculate the errorbars!"
%    deltaR_model - Complex number array. Real part is the model V(in),
%                   imaginary part is the model V(out). Represents change
%                   in reflectance from pump heating.
%                  If twofit is true, this is a 2-element cell array.
%    ratio_model - Ratio signal -V(in)/V(out) from the thermal model.
%                  If twofit is true, this is a 2-element cell array.
%    Mcell       - Self-consistent laser heating can change the thermal
%                  parameters depending on the fit parameters, so the
%                  final Mcell is also an output of this function.
%    T_adj       - T_adj = T0 + dTss + dTpp, the "actual" temperature
%                  in Kelvin adjusted for steady-state (and per-pulse
%                  if perpulse is TRUE) heating.
%    dr_model    - model data for ratio-over-ratio fitting (rorfit).
%
% Example:
%    --
%
% Other m-files required: TDTR_REFL_TTM_vH2.m, TDTR_TEMP_TTM_vH2.m, 
%                         TDTR_TDEP_TTM_vH2.m, SS_Heating_TTM_vH2.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_FIT_TTM_vH2.m, TDTR_MANUALFIT_vH2.m, TDTR_FIT_vH2.m

% Author: Gregory Hohensee
% Acknowledgement: built from TDTR_FIT_V4B, my bi-directional tweak to the
% original TDTR_FIT_V4.m, courtesy of Joseph P. Feser.
% University of Illinois at Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 7-Apr-2014 - version vH1
%                   8-Apr-2014 - more comments, harmonized with TTM
%                   14-July-2014 - vH2, twofit functionality
%------------- BEGIN CODE --------------
%% Check input parameters, assign defaults, errors, warnings as necessary
INITIALIZE_CELLPARAMS_TTM_vH2;

%% Assign fit variables X according to their XguessIJC index.
% XguessIJC = [X,i,j,c] indexes X across (i,j) of the Mcell{c}.
Xijc = XguessIJC(:,2:4);
nrows = length(Xijc(:,1));

for x = 1:nrows % iterate over fit parameters
    c = Xijc(x,3);
    Xij = Xijc(:,1:2);
    switch c % switch over Mcell structures
        case 1, LCTE(Xij(x,1),Xij(x,2)) = XguessIJC(x,1); 
        case 2, LCTEG(Xij(x,1),Xij(x,2)) = XguessIJC(x,1); 
        case 3, gg(Xij(x,1),Xij(x,2)) = XguessIJC(x,1);
    end
end
Mcell = {LCTE LCTEG gg}; % update
matparams{1} = Mcell; % update
%% User input loop
done = 0;
while done ~= 1    
    %% Self-consistent steady-state (and optionally per-pulse) heating 
    if T0 ~= -1
        [dTss, dTpp, Mcell] = TDTR_TDEP_TTM_vH2(matparams,sysparams,...
                                               Tparams,intscheme,nnodes);
        fprintf('T0 = %0.2f K, dTss = %0.2f, dTpp = %0.2f\n',T0,dTss,dTpp)
    end
    matparams{1} = Mcell; % update matparams after TDEP
    %% Compute model
    sysparams{2} = f(1); % first frequency
    [deltaR_model,ratio_model]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,...
                                            A_pump,intscheme,nnodes);
    dr_model = 0; % default
    if rorfit
        sysparams{2} = f(2); % second frequency
        [~,ratio_model2]=TDTR_REFL_TTM_vH2(tdelay,matparams,sysparams,...
                                            A_pump,intscheme,nnodes);
        sysparams{2} = f; % return to reference
        dr_model = ratio_model ./ ratio_model2;
    end
    
    %% Update the data and fit comparison figure
    figure(10)
    if rorfit
        semilogx(tdelay,dr_data,'ok',tdelay,dr_model,'r'); 
        xlabel('Time delay (ps)','FontSize',16);
        ylabel('R_1 / R_2','FontSize',16);
        axis([1e-10 10e-9 0 2]);
    elseif twofit
        semilogx(tdelay,ratio_data, 'ok', tdelay,ratio_model,'r',...
                 tdelay,ratio_data2,'or',tdelay,ratio_model2,'b');  

        ylabel('R_1 and R_2','FontSize',16);
        axis([1e-10 10e-9 min([ratio_data;ratio_model;ratio_data2;ratio_model2])*0.8 ...
                          max([ratio_data;ratio_model;ratio_data2;ratio_model2])*1.2]);
        legend('D1','M1','D2','M2');    
    else % ratio fit
        loglog(tdelay,ratio_data,'ob',tdelay,ratio_model,'r');
        xlabel('Time delay (ps)','FontSize',16);
        ylabel('Ratio','FontSize',16);
        axis([1e-10 10e-9 min([ratio_data;ratio_model])/5 max([ratio_data;ratio_model])*5])
        set(gca, 'YTick', [0.5, 1, 2, 5, 10, 20, 50, 100]);
        set(gca, 'YTickLabel', [0.5, 1, 2, 5, 10, 20, 50, 100]);
    end
    set(gca, 'XTick', [1e-11, 2e-11, 5e-11, 1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
    set(gca, 'XTickLabel', [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4]);
    set(gca, 'TickLength' , [.02 .02]);
    set(gca,'XMinorTick','off')
    set(gca,'YMinorTick','off')
    set(gca,'FontSize',16);
    
    %% Goodness-of-fit by fractional residuals, starting from Zdelay ps
    if rorfit
        Lr = length(dr_model);
        res=(1-(dr_model(Zind:Lr) ...
                ./dr_data(Zind:Lr)) ).^2;
    elseif twofit % res = (1 - (R1_model / R1_data))^2 + (1 - (R2_model / R2_data))^2
        L1 = length(ratio_model);
        L2 = length(ratio_model2);
        res=(1-(ratio_model(Zind:L1)./ratio_data(Zind:L1))).^2 ...
                + (1-(ratio_model2(Zind:L2) ./ratio_data2(Zind:L2))).^2;
            
    else % ordinary ratio fitting
        L1 = length(ratio_model);
        res=(1-(ratio_model(Zind:L1) ...
                ./ratio_data(Zind:L1) )).^2;
    end
    Z=sum(res)
    
    %% Inform user of current fit parameters
    
    % Define symbols and units for all possible fit parameters
    tag = {'L','C','t','eta','L*','C*','t*','eta*','G*','g'};
    units = {'W/m-K', 'J/cm^3-K', 'nm', '(Lx/Lz)',...
            'W/m-K', 'J/cm^3-K', 'mm', '(Lx/Lz)','MW/m^3-K',...
            'pW/nm^3-K'};
    scale = [1 1e-6 1e9 1 1 1e-6 1e3 1 1e-3 1];
    
    % Report current fit parameters
    Mat = sprintf('Current material fit parameters (*: %i-channel substrate)',length(LCTEG(1,:)));
    for x = 1:nrows % iterate through fit parameters
        c = Xijc(x,3);
        switch c
            case 3 % gg
                itag = 10;
                Mat = char(Mat,sprintf('%s(%i,%i) = %0.3f %s',...
                      tag{itag},Xijc(x,1),Xijc(x,2),...
                      gg(Xijc(x,1),Xijc(x,2))/scale(itag),units{itag}));
            otherwise % LCTE or LCTEG
                itag = Xijc(x,1) + 4*(c-1);
                if itag <= 4 % LCTE parameters
                    Mat = char(Mat,sprintf('%s(%i) = %0.4f %s',...
                      tag{itag},Xijc(x,1),...
                      LCTE(Xijc(x,1),Xijc(x,2))/scale(itag),units{itag}));
                else % itag between 5 and 9, LCTEG
                    Mat = char(Mat,sprintf('%s(%i) = %0.4f %s',...
                      tag{itag},Xijc(x,1),...
                      LCTEG(Xijc(x,1),Xijc(x,2))/scale(itag),units{itag}));
                end
        end
    end
    Mat
    
    %% get and clean input
    done = input('Enter 1 if done, 2 for a sensitivity plot; else hit "Enter": ');
    if isempty(done) 
        done = 0;
    else
        if done ~= 1 && done ~= 2 && done ~= 0
            done = 0; fprintf('Hey! Invalid input. Go home, you are drunk.\n')
        end
    end
    %% execute user choice
    switch done
        case 2 % sense plot
            if twofit
                % For twofit scenario, just produce two sensitivity
                % plots, one for each frequency.
                sysparams{2} = f(1);
                fignum = 301;
                fprintf('Calculating sensitivities at first frequency...\n')
                senseplot_TTM_vH2(datparams,sysparams, calparams, matparams, Tparams, fignum);
                
                sysparams{2} = f(2);
                fignum = 302;
                fprintf('Calculating sensitivities at second frequency...\n')
                senseplot_TTM_vH2(datparams,sysparams, calparams, matparams, Tparams, fignum);
            else
                fignum = 202;
                senseplot_TTM_vH2(datparams,sysparams, calparams, matparams, Tparams, fignum);
            end
            
            savesens = input('Enter 1 to save sensitivity plot(s) to datadir: ');
            if savesens == 1
                tag = input('Name the sensitivity plot: ','s');
                if twofit
                    % For twofit scenario, just produce two sensitivity
                    % plots, one for each frequency.
                    save(strcat(datadir,'/SENS_', tag, '_F1.mat'))
                    fignum=301;
                    figure(fignum)
                    saveas(fignum, strcat(datadir,'/SENS_',tag,'_F1.fig'))
                    print(fignum,'-depsc',strcat(datadir,'/SENS_',tag,'_F1.eps'))
                    
                    save(strcat(datadir,'/SENS_', tag, '_F2.mat'))
                    fignum=302;
                    figure(fignum)
                    saveas(fignum, strcat(datadir,'/SENS_',tag,'_F2.fig'))
                    print(fignum,'-depsc',strcat(datadir,'/SENS_',tag,'_F2.eps'))
                else
                    save(strcat(datadir,'/SENS_', tag, '.mat'))
                    fignum=202;
                    figure(fignum)
                    saveas(fignum, strcat(datadir,'/SENS_',tag,'.fig'))
                    print(fignum,'-depsc',strcat(datadir,'/SENS_',tag,'.eps'))
                end
            end
        case 0 % manual fitting
            for x = 1:nrows % iterate through fit parameters
                c = Xijc(x,3);
                
                switch c
                    case 3 % coupling parameter g = gg(i,j).
                        itag = 10;
                        temp = input(sprintf('Adjust coupling %s(%i,%i): ',...
                                       tag{itag},Xijc(x,1),Xijc(x,2)));
                        if ~isempty(temp)
                            XguessIJC(x,1) = temp * scale(itag);
                            gg(Xijc(x,1),Xijc(x,2)) = temp * scale(itag);
                        end
                    otherwise % LCTE or LCTEG
                        itag = Xijc(x,1) + 4*(c-1);
                        temp = input(sprintf('Adjust parameter %s(%i): ',...
                                       tag{itag},Xijc(x,3)));
                        if ~isempty(temp)
                            XguessIJC(x,1) = temp * scale(itag);
                            
                            % if adjusting conductivity, then eta may
                            % change, depending on aniso.
                            if Xijc(x,1) == 1 && aniso(Xijc(x,2))
                                if itag <= 4
                                    LCTE(4,Xijc(x,2)) = LCTE(4,Xijc(x,2)) * LCTE(1,Xijc(x,2)) / temp;
                                else
                                    LCTEG(4,Xijc(x,2)) = LCTEG(4,Xijc(x,2)) * LCTEG(1,Xijc(x,2)) / temp;
                                end
                            end
                            
                            % having used LCTE(1,j) to update eta, can now
                            % update LCTE(1,j).
                            if itag <= 4 % LCTE parameters
                                LCTE(Xijc(x,1),Xijc(x,2)) = temp * scale(itag);
                            else % itag between 5 and 9, LCTEG
                                LCTEG(Xijc(x,1),Xijc(x,2)) = temp * scale(itag);
                            end
                        end
                end % end switch over Mcell
            end % end iteration through fit parameters
    end % end switch for execution of user's decision
    
    Mcell = {LCTE LCTEG gg}; % update Mcell now.
    matparams{1} = Mcell; % update matparams.
end % end user input loop

if twofit
    % combine models into cell arrays for convenient output.
    ratio_model = {ratio_model,ratio_model2};
    deltaR_model = {deltaR_model,deltaR_model2};
end

% Assign final temperature adjusted for SS and PP heating.
if T0 ~= -1, T_adj = T0 + dTss + dTpp; else T_adj = T0; end

% Assign Xsol values
Xsol = XguessIJC(:,1);
end % end program
