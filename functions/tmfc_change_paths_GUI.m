function tmfc_change_paths_GUI(SPM_paths,old_path,new_path)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Change paths inside SPM.mat files using SPM's spm_changepath.
% 
% Open a GUI window for changing paths. 
% FORMAT: tmfc_change_paths_GUI()
% FORMAT: tmfc_change_paths_GUI(SPM_paths)
%
% Command-line (no GUI) mode.
% FORMAT: tmfc_change_paths_GUI(paths,old_path,new_path)
%
% INPUTS:
%   paths     - cell array containing paths to SPM.mat files to update
%   old_path  - string; old root/pattern to replace
%   new_path  - string; new root/pattern to insert
%
% EXAMPLES:
%   Non-GUI:
%   tmfc_change_paths_GUI(SPM_paths, 'C:\Project\Subjects', 'E:\Archive\Project\Subjects');
%
%   GUI, with pre-supplied list:
%   tmfc_change_paths_GUI(SPM_paths);
%
%   GUI, derive from tmfc struct:
%   tmfc_change_paths_GUI({tmfc.subjects.path});
%
%   GUI, ask user to pick SPM.mat files:
%   tmfc_change_paths_GUI();
%
% NOTE: SPM will backup original files as SPM.mat.old
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

% Select SPM.mat files
if nargin<1 || isempty(SPM_paths)
    try
        [SPM_paths, ~] = tmfc_select_subjects_GUI(0);
    catch
        error('Could not open subject selection GUI. Provide ''paths'' explicitly.');
    end
end

% Accept single string
if ischar(SPM_paths) || (isstring(SPM_paths) && isscalar(SPM_paths))
    SPM_paths = cellstr(SPM_paths);
end

% Ensure cellstr
if isstring(SPM_paths)
    SPM_paths = cellstr(SPM_paths);
end
if ~iscellstr(SPM_paths)
    error('SPM.mat paths must be a cell array of character vectors.');
end

% Check SPM paths
SPM_paths = strtrim(SPM_paths(:));
for i = 1:numel(SPM_paths)
    if exist(SPM_paths{i},'dir')
        SPM_paths{i} = fullfile(SPM_paths{i},'SPM.mat');
    end
end
exists_mask = cellfun(@(p) exist(p,'file')==2, SPM_paths);
if ~all(exists_mask)
    warning('Missing SPM.mat files:\n%s', strjoin(SPM_paths(~exists_mask), sprintf('\n')));
    SPM_paths = SPM_paths(exists_mask);
end
if isempty(SPM_paths)
    error('No valid SPM.mat files found.');
end
SPM_paths_char = char(SPM_paths);

%-Command-line mode--------------------------------------------------------
if nargin >= 3 && ~isempty(old_path) && ~isempty(new_path)
    try
        spm_changepath(SPM_paths_char,char(old_path),char(new_path));   
        disp('Paths have been changed.');
    catch 
        warning('Paths have not been changed.');
    end
    return;
end

%-GUI mode-----------------------------------------------------------------
tmfc_CP_MW = figure('Name', 'Change paths', 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [0.30 0.40 0.35 0.26], ...
    'Resize','off','color','w','MenuBar', 'none','ToolBar', 'none');

if isunix; fontscale1 = 0.8; xloc = 1.2; else; fontscale1 = 1; xloc = 1; end

tmfc_CP_MW_S1 = uicontrol(tmfc_CP_MW,'Style','text','String', ...
    'Change paths in SPM.mat files','Units', 'normalized','fontunits','normalized', 'fontSize', 0.20, ...
    'Position', [0.18 0.66 0.660 0.300],'backgroundcolor',get(tmfc_CP_MW,'color'));

% Old pattern 
tmfc_CP_MW_S2 = uicontrol(tmfc_CP_MW,'Style','text','String', ...
    'Old pattern (e.g., C:\Project_folder\Subjects):','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.20*fontscale1, ...
    'Position', [0.063*xloc 0.55 0.450 0.260],'backgroundcolor',get(tmfc_CP_MW,'color'));

tmfc_CP_MW_E1 = uicontrol(tmfc_CP_MW,'Style','edit','String', '','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.50,'HorizontalAlignment', 'left', ...
    'Position', [0.05 0.60 0.880 0.110]);

% New pattern
tmfc_CP_MW_S3 = uicontrol(tmfc_CP_MW,'Style','text','String', ...
    'New pattern (e.g., E:\All_Projects\Project_folder\Subjects):','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.20*fontscale1, ...
    'Position', [0.058 0.29 0.590 0.260],'backgroundcolor',get(tmfc_CP_MW,'color'));

tmfc_CP_MW_E2 = uicontrol(tmfc_CP_MW,'Style','edit','String', '','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.50,'HorizontalAlignment', 'left', ...
    'Position', [0.05 0.35 0.880 0.110]);

% Backups info
tmfc_CP_MW_S4 = uicontrol(tmfc_CP_MW,'Style','text','String', ...
    'Backups of original SPM.mat files will be saved with the ''.old'' suffix.', ...
    'Units', 'normalized','fontunits', ...
    'normalized', 'fontSize', 0.20*fontscale1, ...
    'Position', [0.15 0.06 0.70 0.260],'backgroundcolor',get(tmfc_CP_MW,'color'));

tmfc_CP_MW_OK = uicontrol(tmfc_CP_MW,'Style','pushbutton', 'String', 'OK','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.35, ...
    'Position', [0.26 0.05 0.180 0.130],'callback', @execute_change);
tmfc_CP_MW_Help = uicontrol(tmfc_CP_MW,'Style','pushbutton', 'String', 'Help','Units', 'normalized', ...
    'fontunits','normalized', 'fontSize', 0.35, ...
    'Position', [0.54 0.05 0.180 0.130],'callback', @help_window);

movegui(tmfc_CP_MW,'center');
uiwait(tmfc_CP_MW);
if isvalid(tmfc_CP_MW), delete(tmfc_CP_MW); end
return;

%-Nested functions---------------------------------------------------------
function execute_change(~,~)
    old_path_gui = strtrim(get(tmfc_CP_MW_E1, 'String'));
    new_path_gui = strtrim(get(tmfc_CP_MW_E2, 'String'));
    if isempty(old_path_gui) || isempty(new_path_gui)
        warndlg('Please provide both the OLD and NEW patterns.','TMFC: Change paths');
        return;
    end
    try   
        spm_changepath(SPM_paths_char,char(old_path_gui),char(new_path_gui));        
        disp('Paths have been changed.');
    catch 
        warning('Paths have not been changed.');
    end
    uiresume(tmfc_CP_MW);
end

function help_window(~,~)
    tmfc_CP_HW = figure('Name', 'Change paths: Help', 'NumberTitle', 'off', 'Units', 'normalized', ...
        'Position', [0.60 0.40 0.35 0.45],'Resize','off','MenuBar', 'none','ToolBar', 'none','color','w');

    help_S1 = {'Suppose you moved the project folder after first-level model specification and/or estimation.','',...
        'Original SPM.mat file contains old paths to the model directory, functional images, etc.',...
        '',...
        'The TMFC toolbox uses paths recorded in SPM.mat files. To access functional files, you need to update these paths in SPM.mat files.'};
    
    help_S3 = {'Old pattern: C:\Project_folder\Subjects','New pattern: E:\All_Projects\Project_folder\Subjects','',...
        'Old path for the first functional image: ','C:\Project_folder\Subjects\Sub_01\func\swar_001.nii','',...
        'New path for the first functional image: ','E:\All_Projects\Project_folder\Subjects\Sub_01\func\swar_001.nii'};

    if isunix; fontscale = 0.8; else; fontscale = 1; end

    tmfc_CP_HW_S1 = uicontrol(tmfc_CP_HW,'Style','text','String', help_S1,'Units', 'normalized', ...
        'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.09*fontscale, ...
        'backgroundcolor',get(tmfc_CP_HW,'color'), 'Position', [0.04 0.56 0.900 0.400]);

    tmfc_CP_HW_S2 = uicontrol(tmfc_CP_HW,'Style','text','String', 'Example:','Units', 'normalized', ...
        'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.35,'fontweight', 'bold', ...
        'backgroundcolor',get(tmfc_CP_HW,'color'), 'Position', [0.04 0.48 0.200 0.100]);

    tmfc_CP_HW_S3 = uicontrol(tmfc_CP_HW,'Style','text','String', help_S3,'Units', 'normalized', ...
        'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.09*fontscale, ...
        'backgroundcolor',get(tmfc_CP_HW,'color'), 'Position', [0.04 0.08 0.900 0.400]);
    tmfc_CP_HW_OK = uicontrol(tmfc_CP_HW,'Style','pushbutton', 'String', 'OK','Units', 'normalized', ...
        'fontunits','normalized', 'fontSize', 0.35, 'Position', [0.45 0.04 0.140 0.090], ...
        'callback', @close_tmfc_CP_HW);

    movegui(tmfc_CP_HW,'center');
    
    function close_tmfc_CP_HW(~,~)
        if isvalid(tmfc_CP_HW), close(tmfc_CP_HW); end
    end
end
end
