function [GMprob,WMprob,CSFprob,GMdilate,WMerode,CSFerode] = tmfc_masks_GUI()

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Opens a GUI for defining binary mask parameters: 
% (1) Probability thresholds for GM, WM, and CSF maps.
% (2) Number of dilation/erosion cycles for GM, WM, and CSF masks.
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

GMprob = 0.95;
WMprob = 0.99;
CSFprob = 0.99;
GMdilate = 2;
WMerode = 3;
CSFerode = 2;

CSF_MW = figure('Name','Create whole-brain, GM, WM, and CSF masks','NumberTitle','off','Units','normalized','Position',[0.33 0.08 0.40 0.85],'MenuBar','none','ToolBar','none','color','w','CloseRequestFcn',@CSF_MW_exit,'WindowStyle','modal');
movegui(CSF_MW,'center'); 

CSF_str = {'Whole-brain, GM, WM, and CSF masks are calculated based on tissue probability maps generated after segmentation. All masks are normalized to the MNI coordinate space.','',...
    '(1) The whole-brain mask contains voxels that have a non-zero GM probability or WM/CSF probabilities greater than user-specified probability thresholds (99% by default).','',...
    '(2) The GM is obtained by applying a user-specified threshold (95% by default). This mask is used to calculate DVARS (Mascali et al., 2020).','',...
    '(3) The WM mask is obtained by applying a user-specified threshold (99% by default). Next, the WM mask is eroded several times (3 cycles by default). One erosion cycle removes one voxel from the edge of a binary mask. Finally, the WM mask is normalized to the MNI space and deprived of the brainstem (Mascali et al., 2020).','',...
    '(4) The CSF mask is obtained by applying a user-specified threshold (99% by default). To exclude voxels in close proximity to gray matter, the liberal GM mask is subtracted from the CSF mask. The liberal GM mask is obtained by applying a user-specified threshold (95% by default) and dilation (2 cycles by default). Next, the CSF mask is eroded (2 cycles by default). Finally, the CSF mask is normalized to the MNI space and deprived of any non-ventricle structure (Mascali et al., 2020).','',...
    'NOTE: BOLD signals are extracted from the intersection of the whole-brain or GM/WM/CSF masks with the implicit SPM mask created during model estimation. The implicit SPM mask includes voxels with sufficient BOLD signal intensity. By default, SPM includes voxels in the analysis if their signal is above 0.8 of the global mean signal. To create a more liberal SPM mask, the default "Masking threshold" should be reduced (e.g., to 0.4) and the models should be re-estimated prior to TMFC denoise application.'};

if isunix; fontscale = 0.85; else; fontscale = 1; end

CSF_MW_txt = uicontrol(CSF_MW,'Style','text','String',CSF_str,'Units','normalized','Position',[0.04 0.51 0.92 0.46],'FontUnits','normalized','FontSize',0.0332*fontscale,'HorizontalAlignment','Left','backgroundcolor','w');
CSF_MW_P = uipanel(CSF_MW,'Units','normalized','Position',[0.21 0.095 0.55 0.40],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');

GM_prob_TXT = uicontrol(CSF_MW,'Style','text','String','GM probability threshold:','Units','normalized','Position',[0.26 0.437 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
GM_prob_ED = uicontrol(CSF_MW ,'Style','Edit','String',num2str(GMprob),'Units','normalized','Position',[0.55 0.433 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');

WM_prob_TXT = uicontrol(CSF_MW,'Style','text','String','WM probability threshold:','Units','normalized','Position',[0.26 0.376 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
WM_prob_ED = uicontrol(CSF_MW ,'Style','Edit','String',num2str(WMprob),'Units','normalized','Position',[0.55 0.371 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
CSF_prob_TXT = uicontrol(CSF_MW,'Style','text','String','CSF probability threshold:','Units','normalized','Position',[0.26 0.312 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
CSF_prob_ED = uicontrol(CSF_MW ,'Style','Edit','String',num2str(CSFprob),'Units','normalized','Position',[0.55 0.307 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');

GMdilate_cyc_TXT = uicontrol(CSF_MW,'Style','text','String','GM dilation cycles:','Units','normalized','Position',[0.26 0.248 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
GMdilate_cyc_ED = uicontrol(CSF_MW ,'Style','Edit','String',num2str(GMdilate),'Units','normalized','Position',[0.55 0.243 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
WM_erode_TXT = uicontrol(CSF_MW,'Style','text','String','WM erosion cycles:','Units','normalized','Position',[0.26 0.183 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
WM_erode_ED = uicontrol(CSF_MW ,'Style','Edit','String',num2str(WMerode),'Units','normalized','Position',[0.55 0.179 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
CSF_erode_TXT = uicontrol(CSF_MW,'Style','text','String','CSF erosion cycles:','Units','normalized','Position',[0.26 0.120 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
CSF_erode_ED = uicontrol(CSF_MW ,'Style','Edit','String',num2str(CSFerode),'Units','normalized','Position',[0.55 0.115 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
CSF_OK = uicontrol(CSF_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.40 0.02 0.2 0.05],'FontUnits','normalized','FontSize',0.34,'callback',@read_values);

% --- Close (select defaults) ---
function CSF_MW_exit(~,~)
    GMprob = 0.95;
    WMprob = 0.99;
    CSFprob = 0.99;
    GMdilate = 2;
    WMerode = 3;
    CSFerode = 2;
    disp('Default settings have been selected.');
    uiresume(CSF_MW);
end

% --- OK: read & validate ---
function read_values(~,~)
    % helper validators
    toScalar = @(x) isfinite(x) && isscalar(x);
    isProb   = @(x) toScalar(x) && x >= 0 && x <= 1;
    isNNInt  = @(x) toScalar(x) && x >= 0 && x == floor(x);
    
    % GM probability
    temp = str2double(get(GM_prob_ED,'String'));
    if isProb(temp), GMprob = temp; else
        error('The GM probability threshold must be between 0 and 1.');
    end

    % WM probability
    temp = str2double(get(WM_prob_ED,'String'));
    if isProb(temp), WMprob = temp; else
        error('The WM probability threshold must be between 0 and 1.');
    end

    % CSF probability
    temp = str2double(get(CSF_prob_ED,'String'));
    if isProb(temp), CSFprob = temp; else
        error('The CSF probability threshold must be between 0 and 1.');
    end

    % GM dilation cycles 
    temp = str2double(get(GMdilate_cyc_ED,'String'));
    if isNNInt(temp), GMdilate = temp; else
        error('The number of dilation cycles for GM must be a non-negative integer.');
    end

    % WM erosion cycles
    temp = str2double(get(WM_erode_ED,'String'));
    if isNNInt(temp), WMerode = temp; else
        error('The number of erosion cycles for WM must be a non-negative integer.');
    end

    % CSF erosion cycles
    temp = str2double(get(CSF_erode_ED,'String'));
    if isNNInt(temp), CSFerode = temp; else
        error('The number of erosion cycles for CSF must be a non-negative integer.');
    end
    
    uiresume(CSF_MW);
end
uiwait(CSF_MW);
delete(CSF_MW);
end