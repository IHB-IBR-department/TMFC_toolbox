function tmfc_plot_DVARS(preDVARS,postDVARS,FD,options,SPM_paths,subject_paths,anat_paths,func_paths,masks)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Opens a GUI with FD-DVARS plots before and after noise regression.
%
% FORMAT: tmfc_plot_DVARS(preDVARS,postDVARS,FD)
% Allows saving group FD-DVARS statistics only.
%
% FORMAT: tmfc_plot_DVARS(preDVARS,postDVARS,FD,options,SPM_paths,subject_paths,anat_paths,func_paths,masks)
% Allows saving group FD-DVARS statistics and TMFC denoise settings.
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru


% FD-DVARS input only
if nargin == 3
    options = []; SPM_paths = []; subject_paths = []; anat_paths = []; func_paths = []; masks = [];
end

hasRWLS = isstruct(options) && isfield(options,'rWLS') && options.rWLS==1;

if hasRWLS
    DVARS_label = 'DVARS_{GM} (normalized)';        
else
    DVARS_label = 'DVARS_{GM} (de-meaned)';
end

% GUI elements
DVARS_MW = figure('Name','Framewise displacement and DVARS','NumberTitle','off','Units','normalized','Position',[0.25 0.09 0.50 0.80],'MenuBar','none','ToolBar','none','Color','w','CloseRequestFcn',@DVARS_MW_exit);
DVARS_MW_txt = uicontrol(DVARS_MW,'Style','text','String','Select subject:','Units','normalized','Position',[0.075 0.94 0.85 0.038],'fontunits','normalized','FontSize',0.55,'HorizontalAlignment','Left','backgroundcolor','w');
DVARS_MW_LB1 = uicontrol(DVARS_MW,'Style','listbox','String',[],'Max',1,'Value',1,'Units','normalized','Position',[0.075 0.76 0.85 0.180],'FontUnits','points','FontSize',12,'callback',@update_plot);
movegui(DVARS_MW,'center');

% Create axes for FD plot
ax_frame_1 = axes('Parent',DVARS_MW,'Units','normalized','Position',[0.075 0.53 .85 .2],'Color',[0.75 0.75 0.75]);
box(ax_frame_1,'on'); xlabel(ax_frame_1,'Scans','FontSize',9); ylabel(ax_frame_1,'FD, [mm]','FontSize',9); 

% Create axes for DVARS plot
ax_frame_2 = axes('Parent',DVARS_MW,'Units','normalized','Position',[0.075 0.23 .85 .25],'Color',[0.75 0.75 0.75]);
box(ax_frame_2,'on');  xlabel(ax_frame_2,'Scans','FontSize',9); ylabel(ax_frame_2,DVARS_label,'FontSize',9);

% Statistics, Save & Ok buttons
DVARS_MW_LB2 = uicontrol(DVARS_MW,'Style','listbox','Enable','inactive','String',[],'min', 1, 'max', 3,'Value',[],'Units','normalized','Position',[0.075 0.090 0.85 0.075],'FontUnits','points','FontSize',12,'callback',@update_plot);
DVARS_MW_OK = uicontrol(DVARS_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.30 0.03 0.15 0.045],'FontUnits','normalized','FontSize',0.34,'callback',@DVARS_MW_exit);
DVARS_MW_SAVE = uicontrol(DVARS_MW,'Style','pushbutton','String','Save','Units','normalized','Position',[0.55 0.03 0.15 0.045],'FontUnits','normalized','FontSize',0.34,'callback',@save_data);
group_mean_pre_FD_DVARS_corr = 0; group_mean_post_FD_DVARS_corr = 0;
group_SD_pre_FD_DVARS_corr = 0; group_SD_post_FD_DVARS_corr = 0;

% Initial plotting for subject 1
plot_data(1);

% Close GUI
%--------------------------------------------------------------------------
function DVARS_MW_exit(~,~)
    uiresume(DVARS_MW);
end

% Select subject from the list
%--------------------------------------------------------------------------
function update_plot(~,~)
    selected_subject = get(DVARS_MW_LB1,'Value');
    plot_data(selected_subject); 
end

% Update plot for selected subject
%--------------------------------------------------------------------------
function plot_data(iSub)
    
    if exist('ax_frame_1','var') && ishandle(ax_frame_1), delete(ax_frame_1); end
    if exist('ax_frame_2','var') && ishandle(ax_frame_2), delete(ax_frame_2); end
    
    % Prepare FD and DVARS time series -------------------------------------
    FD_ts = []; preDVARS_ts = []; postDVARS_ts = [];
    sess_sum = 0; sess = 0;
    for jSess = 1:length(FD(iSub).Sess)
        FD_ts = [FD_ts; FD(iSub).Sess(jSess).FD_ts];
        if hasRWLS
            preDVARS_ts = [preDVARS_ts; NaN(3,1); tmfc_zscore(preDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
        else
            preDVARS_ts = [preDVARS_ts; NaN(3,1); spm_detrend(preDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
        end
        if ~isempty(postDVARS)
            if hasRWLS
                postDVARS_ts = [postDVARS_ts; NaN(3,1); tmfc_zscore(postDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
            else
                postDVARS_ts = [postDVARS_ts; NaN(3,1); spm_detrend(postDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(4:end-1)); NaN];
            end
        end
        sess_sum = sess_sum + length(FD(iSub).Sess(jSess).FD_ts) + 1;
        sess = [sess; sess_sum];
    end

    % Plotting FD
    ax_frame_1 = axes('Parent',DVARS_MW,'Units','normalized','Position',[0.075 0.53 .85 .2],'Color',[0.75 0.75 0.75]);
    box(ax_frame_1,'on'); xlabel(ax_frame_1,'Scans','FontSize',9); ylabel(ax_frame_1,'FD, [mm]','FontSize',9);
    S1 = plot(ax_frame_1, FD_ts,'Color',[0 0.447 0.7410]); xlabel(ax_frame_1,'Scans'); ylabel(ax_frame_1,'FD, [mm]'); xlim(ax_frame_1,'tight'); x = xlim(ax_frame_1); y = ylim(ax_frame_1); hold(ax_frame_1,'on');
    
    % Plot sessions
    for jSess = 1:length(FD(iSub).Sess)
        if jSess>1
            S1 = plot(ax_frame_1, [sess(jSess) sess(jSess)],[y(1) y(2)],'-k');
        end
        text(ax_frame_1,sess(jSess)+10,y(2),{['Sess ' num2str(jSess)]},'VerticalAlignment','top'); hold(ax_frame_1,'on');
    end
    
    
    % DVARS Plot ----------------------------------------------------------
    if isempty(postDVARS)
        tmp = preDVARS_ts;
    else
        tmp = [preDVARS_ts; postDVARS_ts];
    end
    tmp(isnan(tmp)) = [];
    y1 = max(tmp); y2 = min(tmp);

    ax_frame_2 = axes('Parent',DVARS_MW,'Units','normalized','Position',[0.075 0.23 .85 .25],'Color',[0.75 0.75 0.75]);
    box(ax_frame_2,'on'); xlabel(ax_frame_2,'Scans','FontSize',9); ylabel(ax_frame_2,DVARS_label,'FontSize',9); 
    plot(ax_frame_2,preDVARS_ts,'Color',[0 0.447 0.7410]); hold(ax_frame_2,'on'); xlabel(ax_frame_2,'Scans'); ylabel(ax_frame_2,DVARS_label); xlim(ax_frame_2,'tight');  ylim(ax_frame_2,[y2*1.3 y1*1.1]); x = xlim(ax_frame_2);
    if ~isempty(postDVARS)
        plot(ax_frame_2,postDVARS_ts,'Color',[0.8500 0.3250 0.0980]); hold(ax_frame_2,'on');
    end
    
    for jSess = 1:length(FD(iSub).Sess)
        if jSess>1
            plot(ax_frame_2, [sess(jSess) sess(jSess)],[y2*1.3 y1*1.1],'-k');
        end
        text(ax_frame_2,sess(jSess)+10,y1*1.1,{['Sess ' num2str(jSess)]},'VerticalAlignment','top'); hold(ax_frame_2,'on');
    end
    
    text(ax_frame_2,x(2)+2,y1*0.5,{'Before'},'Color',[0 0.447 0.7410]); hold(ax_frame_2,'on');
    if ~isempty(postDVARS)
        text(ax_frame_2,x(2)+2,0,{'After'},'Color',[0.8500 0.3250 0.0980]);
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
            temp_str = [FD(iSub).Subject ' :: Mean FD-DVARS correlation across sessions: [before/after denoising] = [' ...
                              num2str(round(preDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') '/' ...
                              num2str(round(postDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') ']']; 
        else
            temp_str = [FD(iSub).Subject ' :: Mean FD-DVARS correlation across sessions: [before denoising] = [' ...
                              num2str(round(preDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') ']']; 
        end
        LB1_str = [LB1_str; {temp_str}];
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

    text_info{1,1} = ['Mean (SD) FD-DVARS correlation across subjects before denoising: ' num2str(round(group_mean_pre_FD_DVARS_corr,2),'%.2f') ' (' num2str(round(group_SD_pre_FD_DVARS_corr,2),'%.2f') ')'];
    if ~isempty(postDVARS)
        text_info{1,2} = ['Mean (SD) FD-DVARS correlation across subjects after denoising: ' num2str(round(group_mean_post_FD_DVARS_corr,2),'%.2f') ' (' num2str(round(group_SD_post_FD_DVARS_corr,2),'%.2f') ')'];
    end

    set(DVARS_MW_LB2,'String',text_info);    
end

% Save group statistics & user-specified TMFC denoise settings
%--------------------------------------------------------------------------
function save_data(~,~)
    if isempty(SPM_paths)
        [filename, pathname] = uiputfile('*.mat', 'Save FD-DVARS group statistics');
        if isequal(filename,0) || isequal(pathname,0)
            fprintf(2,'FD-DVARS group statistics not saved: file name or path not selected.\n'); 
        else
            fullpath = fullfile(pathname, filename);
            save(fullpath,'FD','preDVARS','postDVARS', ...
                'group_mean_pre_FD_DVARS_corr','group_mean_post_FD_DVARS_corr', ...
                'group_SD_pre_FD_DVARS_corr','group_SD_post_FD_DVARS_corr');
            fprintf('FD-DVARS group statistics saved: %s\n', fullpath);
        end
    else   
        [filename, pathname] = uiputfile('*.mat', 'Save FD-DVARS group statistics and TMFC denoise settings');
        if isequal(filename,0) || isequal(pathname,0)
            fprintf(2,'FD-DVARS group statistics not saved: file name or path not selected.\n'); 
        else
            fullpath = fullfile(pathname, filename);
            denoising_settings.SPM_paths = SPM_paths;
            denoising_settings.subject_paths = subject_paths;
            denoising_settings.options = options;
            denoising_settings.anat_paths = anat_paths;
            denoising_settings.func_paths = func_paths; 
            denoising_settings.masks = masks; 
            save(fullpath,'denoising_settings','FD','preDVARS','postDVARS', ...
                'group_mean_pre_FD_DVARS_corr','group_mean_post_FD_DVARS_corr', ...
                'group_SD_pre_FD_DVARS_corr','group_SD_post_FD_DVARS_corr');
            fprintf('FD-DVARS group statistics and TMFC denoise settings saved: %s\n', fullpath);
        end
    end
end

uiwait(DVARS_MW);
delete(DVARS_MW);
end