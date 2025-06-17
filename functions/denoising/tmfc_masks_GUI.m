function [GMprob,WMprob,CSFprob,GMdilate,WMerode,CSFerode] = tmfc_masks_GUI()

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Opens a GUI for defininf binary mask parameters: 
% (1) Probability thresholds for GM, WM and CSF maps.
% (2) Number of dilation/erosion cycles for GM, WM and CSF masks.
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

GMprob = 0.95;
WMprob = 0.99;
CSFprob = 0.99;
GMdilate = 2;
WMerode = 3;
CSFerode = 2;

CSF_MW = figure('Name','Create whole-brain, GM, WM, and CSF masks','NumberTitle','off','Units','normalized','Position',[0.33 0.08 0.40 0.85],'MenuBar','none','ToolBar','none','color','w','CloseRequestFcn',@CSF_MW_exit);
movegui(CSF_MW,'center'); 

CSF_str = {'Whole brain, GM, WM and CSF masks are calculated based on tissue probability maps generated after segmentation. All masks are normalized to the MNI coordinate space.','',...
    '(1) The whole-brain mask contains voxels that have a non-zero GM probability or WM/CSF probabilities greater than user-specified probability thresholds (99% by default).','',...
    '(2) The GM is obtained by applying a user-specified threshold (95% by default). This mask is used to calculate DVARS (Mascali et al., 2020).','',...
    '(3) The WM mask is obtained by applying a user-specified threshold (99% by default). Next, the WM mask is eroded several times (3 cycles by default). One erosion cycle removes one voxel from the edge of a binary mask. Finally, the WM mask is normalized to the MNI space and deprived of the brainstem (Mascali et al., 2020).','',...
    '(4) The CSF mask is obtained by applying a user-specified threshold (99% by default). To exclude voxels in close proximity to gray matter, liberal GM mask is subtracted from the CSF mask. The liberal GM mask is obtained by applying a user-specified threshold (95% by default) and dilation (2 cycles by default). Next, the CSF mask is eroded (2 cycles by default). Finally, the CSF mask is normalized to the MNI space and deprived of any nonventricle structure (Mascali et al., 2020).','',...
    'NOTE: BOLD signals are extracted from the intersection of the whole-brain or GM/WM/CSF masks and the implicit SPM mask created during model estimation. The implicit SPM mask contains voxels with sufficient BOLD signal intensity. By default, SPM includes voxels in the analysis if the signal is above 0.8 of the global mean signal. To create more liberal SPM mask, the default "Masking threshold" should be reduced (e.g., to 0.4) and the models should be re-estimated prior to TMFC denoise application.'};

if isunix; fontscale = 0.85; else; fontscale = 1; end

CSF_MW_txt = uicontrol(CSF_MW,'Style','text','String',CSF_str,'Units','normalized','Position',[0.04 0.51 0.92 0.46],'FontUnits','normalized','FontSize',0.0332*fontscale,'HorizontalAlignment','Left','backgroundcolor','w');
CSF_MW_P = uipanel(CSF_MW,'Units','normalized','Position',[0.21 0.095 0.55 0.40],'HighLightColor',[0.78 0.78 0.78],'BackgroundColor','w','BorderType','line');

GM_prob_TXT = uicontrol(CSF_MW,'Style','text','String','GM probability threshold:','Units','normalized','Position',[0.26 0.437 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
GM_prob_ED = uicontrol(CSF_MW ,'Style','Edit','String',GMprob,'Units','normalized','Position',[0.55 0.433 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');

WM_prob_TXT = uicontrol(CSF_MW,'Style','text','String','WM probability threshold:','Units','normalized','Position',[0.26 0.376 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
WM_prob_ED = uicontrol(CSF_MW ,'Style','Edit','String',WMprob,'Units','normalized','Position',[0.55 0.371 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
CSF_prob_TXT = uicontrol(CSF_MW,'Style','text','String','CSF probability threshold:','Units','normalized','Position',[0.26 0.312 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
CSF_prob_ED = uicontrol(CSF_MW ,'Style','Edit','String',CSFprob,'Units','normalized','Position',[0.55 0.307 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');

GMdilate_cyc_TXT = uicontrol(CSF_MW,'Style','text','String','GM dilation cycles:','Units','normalized','Position',[0.26 0.248 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
GMdilate_cyc_ED = uicontrol(CSF_MW ,'Style','Edit','String',GMdilate,'Units','normalized','Position',[0.55 0.243 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
WM_erode_TXT = uicontrol(CSF_MW,'Style','text','String','WM erosion cycles:','Units','normalized','Position',[0.26 0.183 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
WM_erode_ED = uicontrol(CSF_MW ,'Style','Edit','String',WMerode,'Units','normalized','Position',[0.55 0.179 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
CSF_erode_TXT = uicontrol(CSF_MW,'Style','text','String','CSF erosion cycles:','Units','normalized','Position',[0.26 0.120 0.25 0.03],'FontUnits','normalized','FontSize',0.56,'HorizontalAlignment','right','backgroundcolor','w');
CSF_erode_ED = uicontrol(CSF_MW ,'Style','Edit','String',CSFerode,'Units','normalized','Position',[0.55 0.115 0.12 0.044],'FontUnits','normalized','FontSize',0.38,'backgroundcolor','w');
 
CSF_OK = uicontrol(CSF_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.40 0.02 0.2 0.05],'FontUnits','normalized','FontSize',0.34,'callback',@read_values);


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

function read_values(~,~)
    
    % GM probablity
    temp_GMprob = str2double(get(GM_prob_ED,'String'));
    if (temp_GMprob >= 0 && temp_GMprob <= 1)
       GMprob = temp_GMprob;
    else
        error('The GM probablity threshold value must be between 0.0 and 1.0. Please try again.');
    end
    
    % WM probablity
    temp_WMprob = str2double(get(WM_prob_ED,'String'));
    if (temp_WMprob >= 0 && temp_WMprob <= 1)
       WMprob = temp_WMprob;
    else
        error('The WM probablity threshold value must be between 0.0 and 1.0. Please try again.');
    end
    
    % CSF probablity
    temp_CSFprob = str2double(get(CSF_prob_ED,'String'));
    if (temp_CSFprob >= 0 && temp_CSFprob <= 1)
        CSFprob = temp_CSFprob;
    else
        error('The CSF probablity threshold value must be between 0.0 and 1.0. Please try again.');
    end
    
    % GM dialtion cycles
    temp_GMdilate = str2double(get(GMdilate_cyc_ED,'String'));
    if (temp_GMdilate >= 0 && floor(temp_GMdilate) == temp_GMdilate)
        GMdilate = temp_GMdilate;
    else
        error('The number of dilation cycles for GM must be a natural number. Please try again.');
    end
    
    % WM erosion cycles
    temp_WMerode = str2double(get(WM_erode_ED,'String'));
    if (temp_WMerode >= 0 && floor(temp_WMerode) == temp_WMerode)
        WMerode = temp_WMerode;
    else
        error('The number of erosion cycles for WM must be a natural number. Please try again.');
    end
    
    % CSF erosion cycles
    temp_CSFerode = str2double(get(CSF_erode_ED,'String'));
    if (temp_CSFerode >= 0 && floor(temp_CSFerode) == temp_CSFerode)
        CSFerode = temp_CSFerode;
    else
        error('The number of erosion cycles for CSF must be a natural number. Please try again.');
    end
    
    uiresume(CSF_MW);
end
uiwait(CSF_MW);
delete(CSF_MW);
end