function tmfc_plot_DVARS(preDVARS,postDVARS,FD)

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

% Create panel for strings
DVARS_MW_panel = uipanel(DVARS_MW,'Units','normalized','Position',[0.075 0.10 0.85 0.08],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
DVARS_MW_stat_str = uicontrol(DVARS_MW,'Style','text','String',[],'Units','normalized','Position',[0.08 0.101 0.84 0.075],'fontunits','normalized','FontSize',0.28,'HorizontalAlignment','left','backgroundcolor','w');
DVARS_MW_OK = uicontrol(DVARS_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.42 0.03 0.15 0.045],'FontUnits','normalized','FontSize',0.34,'callback',@DVARS_MW_exit);

% Intial plotting for subject 1
plot_data(1);

function DVARS_MW_exit(~,~)
    delete(DVARS_MW);
end

function update_plot(~,~)
    % Function to update plot w.r.t selected subject
   selected_subject = get(DVARS_MW_LB1,'Value');
   plot_data(selected_subject); 
end

function plot_data(iSub)
   
   delete([ax_frame_1,S1]);
   delete([ax_frame_2,S2]);
          
   % Prepare FD and DVARS time series ---------------------------------
   FD_ts = []; preDVARS_ts = []; postDVARS_ts = [];
   sess_sum = 0; sess = 0;
   for jSess = 1:length(FD(iSub).Sess)
       FD_ts = [FD_ts; FD(iSub).Sess(jSess).FD_ts];
       preDVARS_ts = [preDVARS_ts; NaN; spm_detrend(preDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(2:end))];
       postDVARS_ts = [postDVARS_ts; NaN; spm_detrend(postDVARS(iSub).DVARS.Sess(jSess).DVARS_ts(2:end))];
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
   
   
    % DVARS Plot ------------------------------------------------------
    tmp = [preDVARS_ts; postDVARS_ts]; tmp(isnan(tmp)) = [];
    y1 = max(tmp); y2 = min(tmp);
   
    ax_frame_2 = axes(DVARS_MW,'Units','normalized','Position',[0.075 0.23 .85 .25],'color',[0.75 0.75 0.75]);
    box on; xlabel(ax_frame_2,'Scans','FontSize',9); xlabel(ax_frame_2,'DVARS_{GM} (de-meaned)','FontSize',9); S2 = plot([0]);
    plot(preDVARS_ts,'color',[0 0.447 0.7410 0.8]); hold on; xlabel('Scans'); ylabel('DVARS_{GM} (de-meaned)'); xlim tight;  ylim([y2*1.3 y1*1.1]); x = xlim;
    plot(postDVARS_ts,'color',[0.8500 0.3250 0.0980]); hold on; 
    
    for jSess = 1:length(FD(iSub).Sess)
        if jSess>1
            plot([sess(jSess) sess(jSess)],[y2*1.3 y1*1.1],'-k'); 
        end
        text(sess(jSess)+10,y1*1.1,{['Sess ' num2str(jSess)]},'VerticalAlignment','top'); hold on;
    end 

    text(x(2)+5,y1*0.5,{'Before'},'color',[0 0.447 0.7410 0.8]); hold on;
    text(x(2)+5,0,{'After'},'color',[0.8500 0.3250 0.0980]);

    update_txt();        
end

function update_txt(~,~)
    % Function to update subjects list and mean corr text
    LB1_str = {}; temp_str = {};
    for iSub = 1:length(FD)
        temp_str = [FD(iSub).Subject ' :: Mean FD/DVARS correlation across sessions: [before/after denoising] = [' ...
                          num2str(round(preDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') '/' ...
                          num2str(round(postDVARS(iSub).DVARS.Mean_FD_DVARS_corr,2),'%.2f') ']']; 
        LB1_str =  vertcat(LB1_str,temp_str);
    end
    set(DVARS_MW_LB1,'String',LB1_str);
    
    % Mean correlation calculation and text generation
    pre_FD_DVARS_corr = []; post_FD_DVARS_corr = [];
    for iSub = 1:length(FD)
        pre_FD_DVARS_corr = [pre_FD_DVARS_corr preDVARS(iSub).DVARS.Mean_FD_DVARS_corr];
        post_FD_DVARS_corr = [post_FD_DVARS_corr postDVARS(iSub).DVARS.Mean_FD_DVARS_corr];
    end
    mean_pre_FD_DVARS_corr = mean(pre_FD_DVARS_corr);
    SD_pre_FD_DVARS_corr = std(pre_FD_DVARS_corr);
    mean_post_FD_DVARS_corr = mean(post_FD_DVARS_corr);
    SD_post_FD_DVARS_corr = std(post_FD_DVARS_corr);

    text_info{1,1} = ['Mean (SD) FD/DVARS correlation across subjects before denoising: ' num2str(round(mean_pre_FD_DVARS_corr,2),'%.2f') '(' num2str(round(SD_pre_FD_DVARS_corr,2),'%.2f') ')'];
    text_info{1,2} = ['Mean (SD) FD/DVARS correlation across subjects after denoising: ' num2str(round(mean_post_FD_DVARS_corr,2),'%.2f') '(' num2str(round(SD_post_FD_DVARS_corr,2),'%.2f') ')'];
    
    set(DVARS_MW_stat_str,'String',text_info);
end
end