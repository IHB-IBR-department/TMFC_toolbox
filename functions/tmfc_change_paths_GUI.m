function tmfc_change_paths_GUI(paths)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Opens a GUI window for changing paths. Calls the spm_changepath.m function 
% to change all paths specified in SPM.mat files. Overwrites SPM.mat files
% and saves the backup of the original SPM.mat files as SPM.mat.old.
% 
% FORMAT tmfc_change_paths_GUI(paths)
%
% paths - cell array containing paths to SPM.mat files to update 
%  
% If a tmfc structure containing subject paths is already defined:
% tmfc_change_paths_GUI({tmfc.subjects.path})
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

tmfc_CP_MW = figure('Name', 'Change paths', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.30 0.40 0.35 0.26],'Resize','off','color','w','MenuBar', 'none','ToolBar', 'none');

if isunix; fontscale1 = 0.8; xloc = 1.2; else; fontscale1 = 1; xloc = 1; end

tmfc_CP_MW_S1 = uicontrol(tmfc_CP_MW,'Style','text','String', 'Change paths in SPM.mat files','Units', 'normalized','fontunits','normalized', 'fontSize', 0.20, 'Position', [0.18 0.66 0.660 0.300],'backgroundcolor',get(tmfc_CP_MW,'color'));

% Old pattern 
tmfc_CP_MW_S2 = uicontrol(tmfc_CP_MW,'Style','text','String', 'Old pattern (e.g., C:\Project_folder\Subjects):','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.20*fontscale1, ...
    'Position', [0.063*xloc 0.55 0.450 0.260],'backgroundcolor',get(tmfc_CP_MW,'color'));

tmfc_CP_MW_E1 = uicontrol(tmfc_CP_MW,'Style','edit','String', '','Units', 'normalized','fontunits','normalized', 'fontSize', 0.50,'HorizontalAlignment', 'left', 'Position', [0.05 0.60 0.880 0.110]);

% New pattern
tmfc_CP_MW_S3 = uicontrol(tmfc_CP_MW,'Style','text','String', 'New pattern (e.g., E:\All_Projects\Project_folder\Subjects):','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.20*fontscale1, ...
    'Position', [0.058 0.29 0.590 0.260],'backgroundcolor',get(tmfc_CP_MW,'color'));

tmfc_CP_MW_E2 = uicontrol(tmfc_CP_MW,'Style','edit','String', '','Units', 'normalized','fontunits','normalized', 'fontSize', 0.50,'HorizontalAlignment', 'left', 'Position', [0.05 0.35 0.880 0.110]);

% Backups info
tmfc_CP_MW_S4 = uicontrol(tmfc_CP_MW,'Style','text','String', 'Backups of original SPM.mat files will be saved with the ''.old'' suffix.','Units', 'normalized','fontunits', ...
    'normalized', 'fontSize', 0.20*fontscale1, ...
    'Position', [0.15 0.06 0.70 0.260],'backgroundcolor',get(tmfc_CP_MW,'color'));

tmfc_CP_MW_OK = uicontrol(tmfc_CP_MW,'Style','pushbutton', 'String', 'OK','Units', 'normalized','fontunits','normalized', 'fontSize', 0.35, 'Position', [0.26 0.05 0.180 0.130],'callback', @execute_change);
tmfc_CP_MW_Help = uicontrol(tmfc_CP_MW,'Style','pushbutton', 'String', 'Help','Units', 'normalized','fontunits','normalized', 'fontSize', 0.35, 'Position', [0.54 0.05 0.180 0.130],'callback', @help_window);
movegui(tmfc_CP_MW,'center');

function execute_change(~,~)
    old_path = get(tmfc_CP_MW_E1, 'String');
    new_path = get(tmfc_CP_MW_E2, 'String');
    try   
        spm_changepath(char(paths),char(old_path),char(new_path));        
        disp('Paths have been changed.');
    catch 
        disp('Paths have not been changed.');
    end
    uiresume(tmfc_CP_MW);
end


function help_window(~,~)
    tmfc_CP_HW = figure('Name', 'Change paths: Help', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.60 0.40 0.35 0.45],'Resize','off','MenuBar', 'none','ToolBar', 'none','color','w');

    help_S1 = {'Suppose you moved the project folder after first-level model specification and/or estimation.','',...
        'Original SPM.mat file contains old paths to the model directory, functional images, etc.',...
        '',...
        'The TMFC toolbox uses paths recorded in SPM.mat files. To access functional files, you need to update these paths in SPM.mat files.'};
    help_S3 = {'Old pattern: C:\Project_folder\Subjects','New pattern: E:\All_Projects\Project_folder\Subjects','',...
        'Old path for the first functional image: ','C:\Project_folder\Subjects\Sub_01\func\swar_001.nii','',...
        'New path for the first functional image: ','E:\All_Projects\Project_folder\Subjects\Sub_01\func\swar_001.nii'};

    if isunix; fontscale = 0.8; else; fontscale = 1; end

    tmfc_CP_HW_S1 = uicontrol(tmfc_CP_HW,'Style','text','String', help_S1,'Units', 'normalized', 'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.09*fontscale,'backgroundcolor',get(tmfc_CP_HW,'color'), 'Position', [0.04 0.56 0.900 0.400]);
    tmfc_CP_HW_S2 = uicontrol(tmfc_CP_HW,'Style','text','String', 'Example:','Units', 'normalized', 'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.35,'fontweight', 'bold','backgroundcolor',get(tmfc_CP_HW,'color'), 'Position', [0.04 0.48 0.200 0.100]);
    tmfc_CP_HW_S3 = uicontrol(tmfc_CP_HW,'Style','text','String', help_S3,'Units', 'normalized', 'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.09*fontscale,'backgroundcolor',get(tmfc_CP_HW,'color'), 'Position', [0.04 0.08 0.900 0.400]);
    tmfc_CP_HW_OK = uicontrol(tmfc_CP_HW,'Style','pushbutton', 'String', 'OK','Units', 'normalized','fontunits','normalized', 'fontSize', 0.35, 'Position', [0.45 0.04 0.140 0.090],'callback', @close_tmfc_CP_HW);
    movegui(tmfc_CP_HW,'center');
    
    function close_tmfc_CP_HW(~,~)
        close(tmfc_CP_HW);
    end
end

uiwait(tmfc_CP_MW);
delete(tmfc_CP_MW);
end
