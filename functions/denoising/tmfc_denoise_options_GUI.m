function options = tmfc_denoise_options_GUI()

%-Default options
%--------------------------------------------------------------------------
options = struct; 
options.motion = '24HMP';
options.rotation_indx = [4 5 6];
options.rotation_unit = 'rad';
options.head_radius = 50;
options.DVARS = 1;
options.aCompCor = [5 5];
options.aCompCor_ort = 1;
options.spikereg = 0;
options.spikeregFDthr = 0.5;
options.WM_CSF = 'none';
options.GSR = 'none';
% Can be selected via tmfc_denoise_masks_GUI.m:
% options.GMmask.prob = 0.95;
% options.WMmask.prob = 0.99;
% options.CSFmask.prob = 0.99;
% options.GMmask.dilate = 2;
% options.WMmask.erode = 2;
% options.CSFmask.erode = 2;

%-Options GUI
%--------------------------------------------------------------------------
set_HMP = {'Add 6 temporal derivatives and 12 quadratic terms (24HMP)','Add 6 temporal derivatives (12HMP)','Use standard 6 head motion parameters (6HMP)'};
set_FD_reg = {'First three regressors – translation (e.g., SPM12, HCP, fMRIPrep)','First three regressors – rotation (e.g., FSL, AFNI)'};
set_FD_rot = {'Radians (e.g., SPM12, FSL, fMRIPrep)','Degrees (e.g., HCP, AFNI)'};
set_DVARS = {'Calculate DVARS and FD/DVARS correlations','None'};
set_ACC = {'Add fixed number of aCompCor regressors','Add aCompCor regressors explaining 50% of the signal variability in WM and CSF masks (aCompCor50)','None'};
set_SR = {'None','Add spike regressors'};
set_WM_CSM = {'None','Add WM and CSF signals (2Phys)','Add WM and CSF signals along with their temporal derivatives (4Phys)','Add WM and CSF signals, 2 derivatives and 4 quadratic terms (8Phys)'};
set_GSR = {'None','Add whole-brain signal (GSR)','Add whole-brain signal and its temporal derivative (2GSR)','Add whole-brain signal, its temporal derivative and 2 quadratic terms (4GSR)'};
set_ACC_PO = {'w.r.t. HMP and HPF', 'None'};

HMP_str = {'Motion parameters are taken from SPM.mat file (user-specified regressors of no interest, see SPM.Sess.C.C). Temporal derivatives are calculated as backwards differences (Van Dijk et al., 2012). Quadratic terms represent 6 squared motion parameters and 6 squared temporal derivatives (Satterthwaite et al., 2013).'};
FD_str = {'FD is calculated at each time point as the sum of the absolute values of the derivatives of translational and rotational motion parameters (Power et al., 2012).'};
DVARS_str = {'DVARS is computed as the root mean square of the differentiated BOLD time series within the GM mask before and after denoising (Muschelli et al., 2014).'};
AA_str = {'Extract non-neuronal noise-related principal components (PCs) from WM and CSF signals (Behzadi et al., 2007; Muschelli et al., 2014). Perform well in relatively low-motion samples (Parkes et al., 2017). WM and CSF signals can be pre-orthogonalized w.r.t. high-pass filter (HPF) and HMP (Mascali et al., 2020).'};
SR_str = {'Censor high-motion volumes. For each flagged time point, a unit impulse function that had a value of 1 at that time point and 0 elsewhere included as a spike regressor (Satterthwaite et al., 2013). Spike regression with WM/CSF regression perform well in high-motion samples (Parkes et al., 2017).'};
WM_CSM_str = {'Extract averaged WM and CSF signals to account for physiological fluctuations of non-neuronal origin (Fox and Raichle, 2007). Optionally calculate derivatives, squares, and squares of derivatives (Parkes et al., 2017).'};
GSR_str = {'Extract averaged whole-brain signal to account for head motion and physiological fluctuations of non-neuronal origin (Fox et al, 2009). Optionally calculate derivatives, squares, and squares of derivatives (Parkes et al., 2017). GSR may also remove BOLD signal fluctuations of neuronal origin (Chen et al., 2012) and cause the emergence of negative correlations which may be artefactual (Murphy et al., 2009).'};

tmfc_DN_GUI = figure('Name','TMFC denoise: Options','MenuBar','none','ToolBar','none','NumberTitle','off','Units','norm','Position',[0.224 0.07 0.59 0.850],'color','w','Tag','TMFC_DN_GUI','resize','on','CloseRequestFcn',@close_options_GUI);
movegui(tmfc_DN_GUI,'center');

DN_MP_1 = uipanel(tmfc_DN_GUI,'Units','normalized','Position',[0.03 0.88 0.94 0.115],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
DN_MP_2 = uipanel(tmfc_DN_GUI,'Units','normalized','Position',[0.03 0.743 0.94 0.13],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
DN_MP_3 = uipanel(tmfc_DN_GUI,'Units','normalized','Position',[0.03 0.628 0.94 0.11],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
DN_MP_4 = uipanel(tmfc_DN_GUI,'Units','normalized','Position',[0.03 0.458 0.94 0.162],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
DN_MP_5 = uipanel(tmfc_DN_GUI,'Units','normalized','Position',[0.03 0.323 0.94 0.128],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
DN_MP_6 = uipanel(tmfc_DN_GUI,'Units','normalized','Position',[0.03 0.192 0.94 0.124],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');
DN_MP_7 = uipanel(tmfc_DN_GUI,'Units','normalized','Position',[0.03 0.054 0.94 0.132],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');

% Head motion parameters (HMP) GUI elements
DN_HMP_txt_1 = uicontrol(tmfc_DN_GUI,'Style','text','String','Head motion parameters (HMP)','Units','normalized','Position',[0.048 0.964 0.4 0.021],'fontunits','normalized','fontsize',0.80,'HorizontalAlignment','left','fontweight','bold','backgroundcolor','w');
DN_HMP_pop = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_HMP,'Units','normalized','Position',[0.048 0.887 0.90 0.075],'fontunits','normalized','fontsize',0.207);
DN_HMP_txt_2 = uicontrol(tmfc_DN_GUI,'Style','text','String',HMP_str,'Units','normalized','Position',[0.048 0.885 0.90 0.043],'fontunits','normalized','fontsize',0.35,'HorizontalAlignment','left','backgroundcolor','w');

% Framewise Displacement (FD) GUI elements
DN_FD_txt_1 = uicontrol(tmfc_DN_GUI,'Style','text','String','Framewise displacement (FD)','Units','normalized','Position',[0.048 0.842 0.4 0.021],'fontunits','normalized','fontsize',0.80,'HorizontalAlignment','left','fontweight','bold','backgroundcolor','w');
DN_FD_txt_2 = uicontrol(tmfc_DN_GUI,'Style','text','String',FD_str,'Units','normalized','Position',[0.048 0.815 0.90 0.022],'fontunits','normalized','fontsize',0.7,'HorizontalAlignment','left','backgroundcolor','w');
DN_FD_txt_3 = uicontrol(tmfc_DN_GUI,'Style','text','String','The order of motion regressors:','Units','normalized','Position',[0.048 0.788 0.40 0.024],'fontunits','normalized','fontsize',0.625,'HorizontalAlignment','left','backgroundcolor','w');
DN_FD_txt_4 = uicontrol(tmfc_DN_GUI,'Style','text','String','Rotation units:','Units','normalized','Position',[0.525 0.788 0.40 0.024],'fontunits','normalized','fontsize',0.625,'HorizontalAlignment','left','backgroundcolor','w');
DN_FD_pop_1 = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_FD_reg,'Units','normalized','Position',[0.048 0.714 0.42 0.075],'fontunits','normalized','fontsize',0.208);
DN_FD_pop_2 = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_FD_rot,'Units','normalized','Position',[0.525 0.714 0.42 0.075],'fontunits','normalized','fontsize',0.208);

% DVARS GUI elements
DN_DVARS_txt_1 = uicontrol(tmfc_DN_GUI,'Style','text','String','Derivative of root mean square variance over voxels (DVARS)','Units','normalized','Position',[0.048 0.705 0.45 0.021],'fontunits','normalized','fontsize',0.80,'HorizontalAlignment','left','fontweight','bold','backgroundcolor','w');
DN_DVARS_txt_2 = uicontrol(tmfc_DN_GUI,'Style','text','String',DVARS_str,'Units','normalized','Position',[0.048 0.675 0.90 0.022],'fontunits','normalized','fontsize',0.7,'HorizontalAlignment','left','backgroundcolor','w');
DN_DVARS_pop = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_DVARS,'Units','normalized','Position',[0.048 0.595 0.90 0.075],'fontunits','normalized','fontsize',0.208);

% Anatomical CompCor (ACC) GUI elements
DN_ACC_txt_1 = uicontrol(tmfc_DN_GUI,'Style','text','String','Anatomical component correction (aCompCor)','Units','normalized','Position',[0.048 0.588 0.45 0.021],'fontunits','normalized','fontsize',0.80,'HorizontalAlignment','left','fontweight','bold','backgroundcolor','w');
DN_ACC_txt_2 = uicontrol(tmfc_DN_GUI,'Style','text','String',AA_str,'Units','normalized','Position',[0.048 0.54 0.90 0.043],'fontunits','normalized','fontsize',0.35,'HorizontalAlignment','left','backgroundcolor','w');
DN_ACC_pop = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_ACC,'Units','normalized','Position',[0.048 0.464 0.90 0.075],'fontunits','normalized','fontsize',0.208,'callback',@ACC_value);
DN_ACC_txt_3 = uicontrol(tmfc_DN_GUI,'Style','text','String','Number of PCs for WM:','Units','normalized','Position',[0.048 0.47 0.15 0.024],'fontunits','normalized','fontsize',0.625,'HorizontalAlignment','left','backgroundcolor','w');
DN_ACC_txt_4 = uicontrol(tmfc_DN_GUI,'Style','text','String','Number of PCs for CSF:','Units','normalized','Position',[0.3 0.47 0.15 0.024],'fontunits','normalized','fontsize',0.625,'HorizontalAlignment','left','backgroundcolor','w');
DN_ACC_txt_5 = uicontrol(tmfc_DN_GUI,'Style','text','String','Pre-orthogonalize:','Units','normalized','Position',[0.55 0.47 0.13 0.024],'fontunits','normalized','fontsize',0.625,'HorizontalAlignment','left','backgroundcolor','w');
DN_ACC_PO_pop = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_ACC_PO,'Units','normalized','Position',[0.67 0.422 0.278 0.075],'fontunits','normalized','fontsize',0.208,'callback',@ACC_PO_value);
DN_ACC_E1 = uicontrol(tmfc_DN_GUI,'Style','edit','String',options.aCompCor(1),'Units','normalized','HorizontalAlignment','center','Position',[0.200 0.468 0.08 0.03],'fontunits','normalized','fontsize',0.55);
DN_ACC_E2 = uicontrol(tmfc_DN_GUI,'Style','edit','String',options.aCompCor(2),'Units','normalized','HorizontalAlignment','center','Position',[0.45 0.468 0.08 0.03],'fontunits','normalized','fontsize',0.55);

% Spike regression (SR) GUI elements
DN_SR_txt_1 = uicontrol(tmfc_DN_GUI,'Style','text','String','Spike regression (SpikeReg)','Units','normalized','Position',[0.048 0.42 0.45 0.021],'fontunits','normalized','fontsize',0.80,'HorizontalAlignment','left','fontweight','bold','backgroundcolor','w');
DN_SR_txt_2 = uicontrol(tmfc_DN_GUI,'Style','text','String',SR_str,'Units','normalized','Position',[0.048 0.370 0.90 0.043],'fontunits','normalized','fontsize',0.35,'HorizontalAlignment','left','backgroundcolor','w');
DN_SR_pop = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_SR,'Units','normalized','Position',[0.048 0.29 0.90 0.075],'fontunits','normalized','fontsize',0.208);

% White matter (WM) and  Cerebral spinal fuild (CSF) signal regression
DN_WM_CSM_txt_1 = uicontrol(tmfc_DN_GUI,'Style','text','String','WM and CSF signal regression (Phys)','Units','normalized','Position',[0.048 0.284 0.55 0.021],'fontunits','normalized','fontsize',0.80,'HorizontalAlignment','left','fontweight','bold','backgroundcolor','w');
DN_WM_CSM_txt_2 = uicontrol(tmfc_DN_GUI,'Style','text','String',WM_CSM_str,'Units','normalized','Position',[0.048 0.236 0.90 0.043],'fontunits','normalized','fontsize',0.35,'HorizontalAlignment','left','backgroundcolor','w');
DN_WM_CSM_pop = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_WM_CSM,'Units','normalized','Position',[0.048 0.16 0.90 0.075],'fontunits','normalized','fontsize',0.208);
%  
% Global signal regressions (GSR)
DN_GSR_txt_1 = uicontrol(tmfc_DN_GUI,'Style','text','String','Global signal regression (GSR)','Units','normalized','Position',[0.048 0.155 0.55 0.021],'fontunits','normalized','fontsize',0.80,'HorizontalAlignment','left','fontweight','bold','backgroundcolor','w');
DN_GSR_txt_2 = uicontrol(tmfc_DN_GUI,'Style','text','String',GSR_str,'Units','normalized','Position',[0.048 0.096 0.90 0.058],'fontunits','normalized','fontsize',0.25,'HorizontalAlignment','left','backgroundcolor','w');
DN_GSR_pop = uicontrol(tmfc_DN_GUI,'Style','popupmenu','String',set_GSR,'Units','normalized','Position',[0.048 0.020 0.90 0.075],'fontunits','normalized','fontsize',0.208);

% OK Button 
DN_OK = uicontrol(tmfc_DN_GUI,'Style','pushbutton','String','OK','Units','normalized','Position',[0.4 0.0080 0.180 0.04],'FontUnits','normalized','fontsize',0.38,'callback',@export_options);

% Enter the number of PCs for aCompCor
function ACC_value(~,~)
    approach = (DN_ACC_pop.String{DN_ACC_pop.Value});
    if strcmp(approach,'Add fixed number of aCompCor regressors')
        set(DN_ACC_E1,'enable','on');
        set(DN_ACC_E2,'enable','on');
        set(DN_ACC_txt_3,'enable','on');
        set(DN_ACC_txt_4,'enable','on');
        set(DN_ACC_E1,'String','5');
        set(DN_ACC_E2,'String','5');
        options.aCompCor = [5 5];
    elseif strcmp(approach, 'Add aCompCor regressors explaining 50% of the signal variability in WM and CSF masks (aCompCor50)')
        set(DN_ACC_E1,'enable','off');
        set(DN_ACC_E2,'enable','off');
        set(DN_ACC_txt_3,'enable','off');
        set(DN_ACC_txt_4,'enable','off');
        set(DN_ACC_E1,'String','0.5');
        set(DN_ACC_E2,'String','0.5');
        options.aCompCor = [0.5 0.5];
    elseif strcmp(approach,'None')
        set(DN_ACC_E1,'enable','off');
        set(DN_ACC_E2,'enable','off');
        set(DN_ACC_txt_3,'enable','off');
        set(DN_ACC_txt_4,'enable','off');
        set(DN_ACC_E1,'String','0');
        set(DN_ACC_E2,'String','0');
        options.aCompCor = [0 0];
    end
end


function ACC_PO_value(~,~)
    approach = (DN_ACC_PO_pop.String{DN_ACC_PO_pop.Value});
    if strcmp(approach, 'w.r.t. HMP and HPF')
        options.aCompCor_ort = 1;
    else
        options.aCompCor_ort = 0;
    end    
end

% Close GUI
function close_options_GUI(~,~)
    options = [];
    delete(tmfc_DN_GUI);
end

% Export options
function export_options(~,~)
    
    % Head motion parameters 
    HMP_select{1} = get(DN_HMP_pop,'String');
    HMP_select{2} = get(DN_HMP_pop,'Value');    
    if strcmp(HMP_select{1}(HMP_select{2}),'Add 6 temporal derivatives and 12 quadratic terms (24HMP)')
    	options.motion = '24HMP';
    elseif strcmp(HMP_select{1}(HMP_select{2}),'Add 6 temporal derivatives (12HMP)')
    	options.motion = '12HMP';       
    elseif strcmp(HMP_select{1}(HMP_select{2}),'Use standard 6 head motion parameters (6HMP)')
    	options.motion = '6HMP';       
    end
    clear HMP_select
    
    % Framewise displacement
    FD_select_1{1} = get(DN_FD_pop_1,'String');
    FD_select_1{2} = get(DN_FD_pop_1,'Value');            
    if strcmp(FD_select_1{1}(FD_select_1{2}),'First three regressors – translation (e.g., SPM12, HCP, fMRIPrep)')
    	options.rotation_indx = [4 5 6];
    elseif strcmp(FD_select_1{1}(FD_select_1{2}),'First three regressors – rotation (e.g., FSL, AFNI)')
    	options.rotation_indx = [1 2 3];                   
    end
    clear FD_select_1

    FD_select_2{1} = get(DN_FD_pop_2,'String');
    FD_select_2{2} = get(DN_FD_pop_2,'Value');            
    if strcmp(FD_select_2{1}(FD_select_2{2}),'Radians (e.g., SPM12, FSL, fMRIPrep)')
    	options.rotation_unit = 'rad';
    elseif strcmp(FD_select_2{1}(FD_select_2{2}),'Degrees (e.g., HCP, AFNI)')
    	options.rotation_unit = 'deg';                   
    end
    clear FD_select_2    
    
    % DVARS
    DVARS_select{1} = get(DN_DVARS_pop, 'String');
    DVARS_select{2} = get(DN_DVARS_pop, 'Value');
    if strcmp(DVARS_select{1}(DVARS_select{2}), 'Calculate DVARS and FD/DVARS correlations')
        options.DVARS = 1;
    else
        options.DVARS = 0;
    end
    
    
    % aCompCor
    ACC_select{1} = get(DN_ACC_pop,'String');
    ACC_select{2} = get(DN_ACC_pop,'Value');            
    if strcmp(ACC_select{1}(ACC_select{2}),'Add fixed number of aCompCor regressors')
        nPC_WM = str2double(get(DN_ACC_E1,'String'));
        nPC_CSF = str2double(get(DN_ACC_E2,'String'));                  
        if isnan(nPC_WM) || isnan(nPC_CSF)
            error('Please enter natural numbers for number of principle components for WM/CSF.');
        elseif ~(nPC_WM > 0 && floor(nPC_WM) == nPC_WM) || ~(nPC_CSF > 0 && floor(nPC_CSF) == nPC_CSF)
            error('Please enter natural numbers for number of principle components for WM/CSF.');
        elseif (nPC_WM > 100) || (nPC_CSF > 100)
            error('The number of principle components must be between 0 to 100, please re-enter.'); 
        else
            options.aCompCor = [nPC_WM nPC_CSF];
        end
    elseif strcmp(ACC_select{1}(ACC_select{2}),'None')
    	options.aCompCor = [0 0];                   
    end
    clear ACC_select    

    % Spike Regression 
    SR_select{1} = get(DN_SR_pop,'String');
    SR_select{2} = get(DN_SR_pop,'Value');            
    if strcmp(SR_select{1}(SR_select{2}),'None')
    	options.spikereg = 0;
    elseif strcmp(SR_select{1}(SR_select{2}),'Add spike regressors')
    	options.spikereg = 1;                   
    end
    clear SR_select    
    
    % WM/CSF regression
    WM_CSM_select{1} = get(DN_WM_CSM_pop,'String');
    WM_CSM_select{2} = get(DN_WM_CSM_pop,'Value');            
    if strcmp(WM_CSM_select{1}(WM_CSM_select{2}),'None')
    	options.WM_CSF = 'none';
    elseif strcmp(WM_CSM_select{1}(WM_CSM_select{2}),'Add WM and CSF signals (2Phys)')
    	options.WM_CSF = '2Phys';        
    elseif strcmp(WM_CSM_select{1}(WM_CSM_select{2}),'Add WM and CSF signals along with their temporal derivatives (4Phys)')
    	options.WM_CSF = '4Phys';        
    elseif strcmp(WM_CSM_select{1}(WM_CSM_select{2}),'Add WM and CSF signals, 2 derivatives and 4 quadratic terms (8Phys)')
    	options.WM_CSF = '8Phys';                   
    end
    clear WM_CSM_select    
    
    % GSR
    GSR_select{1} = get(DN_GSR_pop,'String');
    GSR_select{2} = get(DN_GSR_pop,'Value');
    if strcmp(GSR_select{1}(GSR_select{2}),'None')
    	options.GSR = 'none';
    elseif strcmp(GSR_select{1}(GSR_select{2}),'Add whole-brain signal (GSR)')
    	options.GSR = 'GSR';        
    elseif strcmp(GSR_select{1}(GSR_select{2}),'Add whole-brain signal and its temporal derivative (2GSR)')
    	options.GSR = '2GSR';        
    elseif strcmp(GSR_select{1}(GSR_select{2}),'Add whole-brain signal, its temporal derivative and 2 quadratic terms (4GSR)')
    	options.GSR = '4GSR';                   
    end
    clear GSR_select
    
    disp('Denoising options selected.');
    delete(tmfc_DN_GUI);
end

uiwait();
end