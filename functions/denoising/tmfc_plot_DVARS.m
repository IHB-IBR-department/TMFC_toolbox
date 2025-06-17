function tmfc_plot_DVARS(preDVARS,postDVARS,FD,options,SPM_paths,subject_paths,struct_paths,funct_paths,masks)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Opens a GUI with FD/DVARS plots before and after noise regression.
%
% -------------------------------------------------------------------------
% FORMAT: tmfc_plot_DVARS(preDVARS,postDVARS,FD)
% Allows to save group FD/DVARS statistics only.
%
% FORMAT: tmfc_plot_DVARS(preDVARS,postDVARS,FD,options,SPM_paths,subject_paths,struct_paths,funct_paths,masks)
% Allows to save group FD statistics anf TMFC denoise settings.
% -------------------------------------------------------------------------
%
% =========================================================================
%
% Copyright (C) 2025 Ruslan Masharipov
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% Contact email: masharipov@ihb.spb.ru


% FD/DVARS input only
if nargin == 3
    options = []; SPM_paths = []; subject_paths = []; struct_paths = []; funct_paths = []; masks = [];
end

% GUI elements
DVARS_MW = figure('Name','Framewise displacement','NumberTitle','off','Units','normalized','Position',[0.25 0.09 0.50 0.80],'MenuBar','none','ToolBar','none','color','w','CloseRequestFcn',@DVARS_MW_exit);
DVARS_MW_txt = uicontrol(DVARS_MW,'Style','text','String','Select subject:','Units','normalized','Position',[0.075 0.94 0.85 0.038],'fontunits','normalized','FontSize',0.55,'HorizontalAlignment','Left','backgroundcolor','w');
DVARS_MW_LB1 = uicontrol(DVARS_MW,'Style','listbox','String',[],'Max',1,'Value',1,'Units','normalized','Position',[0.075 0.76 0.85 0.180],'FontUnits','points','FontSize',12,'callback',@update_plot);
movegui(DVARS_MW,'center');

% Create axes for FD plot
ax_frame_1 = axes(DVARS_MW,'Units','normalized','Position',[0.075 0.53 .85 .2],'color',[0.75 0.75 0.75]);
box on; xlabel(ax_frame_1,'Scans','FontSize',9); ylabel(ax_frame_1,'FD, [mm]','FontSize',9); S1 = plot([0]);

% Create axes for DVARS plot
ax_frame_2 = axes(DVARS_MW,'Units','normalized','Position',[0.075 0.23 .85 .25],'color',[0.75 0.75 0.75]);
box on; xlabel(ax_frame_2,'Scans','FontSize',9); ylabel(ax_frame_2,'FD, [mm]','FontSize',9); S2 = plot([0]);

% Statistics, Save & Ok buttons
DVARS_MW_LB2 = uicontrol(DVARS_MW,'Style','listbox','String',[],'min', 1, 'max', 3,'Value',[],'enable', 'inactive','Units','normalized','Position',[0.075 0.090 0.85 0.075],'FontUnits','points','FontSize',12,'callback',@update_plot);
DVARS_MW_OK = uicontrol(DVARS_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.30 0.03 0.15 0.045],'FontUnits','normalized','FontSize',0.34,'callback',@DVARS_MW_exit);
DVARS_MW_SAVE = uicontrol(DVARS_MW,'Style','pushbutton','String','Save','Units','normalized','Position',[0.55 0.03 0.15 0.045],'FontUnits','normalized','FontSize',0.34,'callback',@save_data);
group_mean_pre_FD_DVARS_corr = 0; group_mean_post_FD_DVARS_corr = 0;
group_SD_pre_FD_DVARS_corr = 0; group_SD_post_FD_DVARS_corr = 0;

% Intial plotting for subject 1
plot_data(1);

% Close GUI
%--------------------------------------------------------------------------
function DVARS_MW_exit(~,~)
    uiresume(DVARS_MW);
end

% Update DVARS plots
%--------------------------------------------------------------------------
function update_plot(~,~)
    % Function to update plot w.r.t selected subject
    selected_subject = get(DVARS_MW_LB1,'Value');
    plot_data(selected_subject); 
end

% Update plot w.r.t Subject
%--------------------------------------------------------------------------
function plot_data(iSub)
    
    delete([ax_frame_1,S1]);
    delete([ax_frame_2,S2]);
    
    % Prepare FD and DVARS time series -------------------------------------
    FD_ts = []; preDVARS_ts = []; postDVARS_ts = [];
    sess_sum = 0; sess = 0;
    for jSess = 1:length(FD(iSub).Sess)
        FD_ts = [FD_ts; FD(iSub).Sess(jSess).FD_ts];
        if options.rWLS == 1
            preDVARS_ts = [preDVARS_ts; NaN(3,1); zscore(preDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
        else
            preDVARS_ts = [preDVARS_ts; NaN(3,1); spm_detrend(preDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
        end
        if ~isempty(postDVARS)
            if options.rWLS == 1
                postDVARS_ts = [postDVARS_ts; NaN(3,1); zscore(postDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
            else
                postDVARS_ts = [postDVARS_ts; NaN(3,1); spm_detrend(postDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
            end
        end
        sess_sum = sess_sum + length(FD(iSub).Sess(jSess).FD_ts) + 1;
        sess = [sess; sess_sum];
    end
    % Plotting DVARS
    ax_frame_1 = axes(DVARS_MW,'Units','normalized','Position',[0.075 0.53 .85 .2],'color',[0.75 0.75 0.75]);
    box on; xlabel(ax_frame_1,'Scans','FontSize',9); ylabel(ax_frame_1,'FD, [mm]','FontSize',9); S1 = plot([0]);
    S1 = plot(FD_ts,'color',[0 0.447 0.7410 0.8]); xlabel('Scans'); ylabel('FD, [mm]'); xlim tight; x = xlim; y = ylim; hold on;
    
    % Plot sessions
    for jSess = 1:length(FD(iSub).Sess)
        if jSess>1
            S1 = plot([sess(jSess) sess(jSess)],[y(1) y(2)],'-k');
        end
        text(sess(jSess)+10,y(2),{['Sess ' num2str(jSess)]},'VerticalAlignment','top'); hold on;
    end
    
    
    % DVARS Plot ----------------------------------------------------------
    if isempty(postDVARS)
        tmp = preDVARS_ts;
    else
        tmp = [preDVARS_ts; postDVARS_ts];
    end
    tmp(isnan(tmp)) = [];
    y1 = max(tmp); y2 = min(tmp);
    
    if options.rWLS == 1
        DVARS_label = 'DVARS_{GM} (normalized)';        
    else
        DVARS_label = 'DVARS_{GM} (de-meaned)';
    end

    ax_frame_2 = axes(DVARS_MW,'Units','normalized','Position',[0.075 0.23 .85 .25],'color',[0.75 0.75 0.75]);
    box on; xlabel(ax_frame_2,'Scans','FontSize',9); xlabel(ax_frame_2,DVARS_label,'FontSize',9); S2 = plot([0]);
    plot(preDVARS_ts,'color',[0 0.447 0.7410 0.8]); hold on; xlabel('Scans'); ylabel(DVARS_label); xlim tight;  ylim([y2*1.3 y1*1.1]); x = xlim;
    if ~isempty(postDVARS)
        plot(postDVARS_ts,'color',[0.8500 0.3250 0.0980]); hold on;
    end
    
    for jSess = 1:length(FD(iSub).Sess)
        if jSess>1
            plot([sess(jSess) sess(jSess)],[y2*1.3 y1*1.1],'-k');
        end
        text(sess(jSess)+10,y1*1.1,{['Sess ' num2str(jSess)]},'VerticalAlignment','top'); hold on;
    end
    
    text(x(2)+5,y1*0.5,{'Before'},'color',[0 0.447 0.7410 0.8]); hold on;
    if ~isempty(postDVARS)
        text(x(2)+5,0,{'After'},'color',[0.8500 0.3250 0.0980]);
    end
    
    update_txt();
end

% Update text info
%--------------------------------------------------------------------------
function update_txt(~,~)
    % Function to update subjects list and mean corr text
    LB1_str = {}; temp_str = {};
    for iSub = 1:length(FD)
        if ~isempty(postDVARS)
            temp_str = [FD(iSub).Subject ' :: Mean FD/DVARS correlation across sessions: [before/after denoising] = [' ...
                              num2str(round(preDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') '/' ...
                              num2str(round(postDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') ']']; 
        else
            temp_str = [FD(iSub).Subject ' :: Mean FD/DVARS correlation across sessions: [before denoising] = [' ...
                              num2str(round(preDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') ']']; 
        end
        LB1_str =  vertcat(LB1_str,temp_str);
    end
    set(DVARS_MW_LB1,'String',LB1_str);
    
    % Mean correlation calculation and text generation
    pre_FD_DVARS_corr = []; post_FD_DVARS_corr = [];
    for iSub = 1:length(FD)
        pre_FD_DVARS_corr = [pre_FD_DVARS_corr preDVARS(iSub).DVARS.Mean_FD_DVARS_corr];
        if ~isempty(postDVARS)
            post_FD_DVARS_corr = [post_FD_DVARS_corr postDVARS(iSub).DVARS.Mean_FD_DVARS_corr];
        end
    end
    group_mean_pre_FD_DVARS_corr = mean(pre_FD_DVARS_corr);
    group_SD_pre_FD_DVARS_corr = std(pre_FD_DVARS_corr);
    if ~isempty(postDVARS)
        group_mean_post_FD_DVARS_corr = mean(post_FD_DVARS_corr);
        group_SD_post_FD_DVARS_corr = std(post_FD_DVARS_corr);
    else
        group_mean_post_FD_DVARS_corr = [];
        group_SD_post_FD_DVARS_corr = [];
    end

    text_info{1,1} = ['Mean (SD) FD/DVARS correlation across subjects before denoising: ' num2str(round(group_mean_pre_FD_DVARS_corr,2),'%.2f') ' (' num2str(round(group_SD_pre_FD_DVARS_corr,2),'%.2f') ')'];
    if ~isempty(postDVARS)
        text_info{1,2} = ['Mean (SD) FD/DVARS correlation across subjects after denoising: ' num2str(round(group_mean_post_FD_DVARS_corr,2),'%.2f') ' (' num2str(round(group_SD_post_FD_DVARS_corr,2),'%.2f') ')'];
    end

    set(DVARS_MW_LB2,'String',text_info);    
end

% Save group statistics & user-specified tmfc_denoise settings
%--------------------------------------------------------------------------
function save_data(~,~)
    if isempty(SPM_paths)
        [filename, pathname] = uiputfile('*.mat', 'Save FD/DVARS group statistics');
        if isequal(filename,0) || isequal(pathname,0)
            fprintf(2,'FD/DVARS group statistics not saved. File name or path not selected.\n'); 
        else
            fullpath = fullfile(pathname, filename);
            save(fullpath,'FD','preDVARS','postDVARS', ...
                'group_mean_pre_FD_DVARS_corr','group_mean_post_FD_DVARS_corr', ...
                'group_SD_pre_FD_DVARS_corr','group_SD_post_FD_DVARS_corr');
            fprintf('FD/DVARS group statistics saved: %s\n', fullpath);
        end
    else   
        [filename, pathname] = uiputfile('*.mat', 'Save FD/DVARS group statistics and TMFC denoise settings');
        if isequal(filename,0) || isequal(pathname,0)
            fprintf(2,'FD/DVARS statistics not saved. File name or path not selected.\n'); 
        else
            fullpath = fullfile(pathname, filename);
            denoising_settings.SPM_paths = SPM_paths;
            denoising_settings.subject_paths = subject_paths;
            denoising_settings.options = options;
            denoising_settings.struct_paths = struct_paths;
            denoising_settings.funct_paths = funct_paths; 
            denoising_settings.masks = masks; 
            save(fullpath,'denoising_settings','FD','preDVARS','postDVARS', ...
                'group_mean_pre_FD_DVARS_corr','group_mean_post_FD_DVARS_corr', ...
                'group_SD_pre_FD_DVARS_corr','group_SD_post_FD_DVARS_corr');
            fprintf('FD/DVARS group statistics and TMFC denoise settings saved: %s\n', fullpath);
        end
    end
end

uiwait(DVARS_MW);
delete(DVARS_MW);
end