function [func_paths] = tmfc_select_func_GUI(SPM_paths,subject_paths)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
%
% Opens a GUI to select unsmoothed, normalized, realigned functional images.
% 
% First, select the parent folder that contains the subject folders
% with FUNC subfolders. This may differ from the parent folder
% that contains the subject folders with STAT subfolders, which is 
% selected by default.
%
% Then select the FUNC subfolder for the first subject and apply text
% filter (e.g., *war*.nii, *wr*.nii, or *preproc*.nii.gz) to match all
% unsmoothed, normalized, realigned functional images.
%
% Alternatively, you can preserve functional image paths from the SPM.mat 
% files (if the GLMs were specified with unsmoothed functional images).
%
% ----------------EXAMPLE #1 (SPM-like folder structure)-------------------
% There is no need to change the parent folder to select functional files.
%
% project/
% ├─ rawdata/   # DICOM
% └─ derivatives/ <---------------------------- [Parent folder with FUNC subfolders (BY DEFAULT)]
%    ├─ sub-01/   <---------------------------------------------------- [Selected subject folder]                                          
%    │  ├─ anat/              
%    │  │  ├─ *T1*.nii  
%    │  │  └─ *T1*.nii derivatives
%    │  ├─ func/  <---------------------------- [Select the FUNC subfolder for the first subject] (1)
%    │  │  ├─ sess-01/
%    │  │  │  ├─ Unprocessed functional files (*.nii)
%    │  │  │  └─ Preprocessed functional files:
%    │  │  │     ● smoothed + normalized + realigned (e.g., swar*.nii)
%    │  │  │     ● unsmoothed + normalized + realigned (e.g., war*.nii) <---- [Apply text filter] (2)
%    │  │  └─ sess-02/ ... 
%    │  └─ stat/                # first-level models (one folder per GLM)
%    │     ├─ GLM-01/
%    │     │  ├─ SPM.mat  <---------------------------------------------- [Selected SPM.mat file]
%    │     │  └─ TMFC_denoise/ <------------------------------------------------- [Output folder]
%    │     └─ GLM-02/ ...
%    └─ sub-02/ ...
%
%
%
% --------------EXAMPLE #2 (BIDS-like folder structure)--------------------
% "project/derivatives/firstlevel-spm" parent folder with STAT subfolders
% needs to be changed to "project/derivatives/fmriprep" — parent folder
% with FUNC subfolders
%
% project/  
% ├── sub-01/
% │   ├── ses-01/                
% │   │   ├── anat/  
% │   │   │   └── Structural file: *T1*.nii  
% │   │   └── func/
% │   │       └── Functional files (unprocessed)
% │   └── ses-02/ ...
% ├── sub-02/ ... 
% └── derivatives
%     ├── fmriprep/ <------------------ [Select parent folder (contains sub-*/ses-*/func)] (1)
%     │   ├── sub-01/
%     │   │   ├── ses-01/
%     │   │   │   └── func/   <--------- [Select the FUNC subfolder for the first subject] (2)
%     │   │   │       └── Preprocessed functional files:
%     │   │   │           ● smoothed + normalized + realigned
%     │   │   │           ● unsmoothed + normalized + realigned <----- [Apply text filter] (3)
%     │   │   └── ses-02/ ...
%     │   └── sub-02/ ...
%     └── firstlevel-spm/ <------------- [Parent folder with FUNC subfolders (BY DEFAULT)] (Needs to be changed!)
%         ├── sub-01/     <------------------------------------- [Selected subject folder]                       
%         │   ├── GLM-01/
%         │   │   ├── SPM.mat        <---------------------------- [Selected SPM.mat file]
%         │   │   └── TMFC_denoise/  <------------------------------------ [Output folder]
%         │   └── GLM-02/ ...
%         └── sub-02/ ...
%
%
%
% --------------EXAMPLE #3 (Other non-BIDS folder structure)---------------
% "project/firstlevel-spm" parent folder with FUNC subfolders needs to be
% changed to "project/nifti" — parent folder with FUNC subfolders
%
% project/
% ├─ rawdata/   # DICOM
% ├─ nifti/ <---------------- [Select parent folder (contains sub-*/ses-*/func)] (1)
% │  ├─ sub-01/                                         
% │  │  ├─ anat/             
% │  │  │  ├─ *T1*.nii 
% │  │  │  └─ *T1*.nii derivatives 
% │  │  └─ func/  <----------- [Select the FUNC subfolder for the first subject] (2) 
% │  │     ├─ sess-01/
% │  │     │  ├─ Unprocessed functional files (*.nii)
% │  │     │  └─ Preprocessed functional files (*.nii):
% │  │     │     ● smoothed + normalized + realigned 
% │  │     │     ● unsmoothed + normalized + realigned <---- [Apply text filter] (3)
% │  │     └─ sess-02/ ... 
% │  └─ sub-02/ ...   
% └─ firstlevel-spm/ <-------- [Parent folder with FUNC subfolders (BY DEFAULT)] (Needs to be changed!)
%    ├─ sub-01/   <----------------------------------- [Selected subject folder]  
%    │  ├─ GLM-01/
%    │  │  ├─ SPM.mat  <-------------------------------- [Selected SPM.mat file]
%    │  │  └─ TMFC_denoise/ <----------------------------------- [Output folder]
%    │  └─ GLM-02/ ...
%    └─ sub-02/ ...
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin < 2, error('Check inputs.'); end
if ischar(subject_paths), subject_paths = cellstr(subject_paths); end
subject_paths = subject_paths(:);  % column cell

func_paths = [];
func_root = '';                 
func_root_default = '';         
func_subfolder_full = '';       
func_subfolder_rel  = '';       
no_files = {};

preserve_spm_paths = false;   
spm_first_images   = {};     

% Mapping container built by Apply Filter:
% unsm_index(iSub).groups(g) with fields:
%   .orig_dir        (directory from SPM.xY.VY)
%   .target_dir      (mapped under FUNC root)
%   .unsm_files      (cellstr of matched files; 1 = 4D, >1 = 3D series)
%   .scan_indices    (indices of scans in SPM.xY.VY belonging to this dir)
unsm_index = struct([]);

% Track whether user changed root/subfolder (to disable fallback to orig_dir)
user_changed_root = false;
user_changed_subf = false;

% GUI elements
% -------------------------------------------------------------------------
SF_MW = figure('Name','Select unsmoothed functional images','NumberTitle','off', ...
    'Units','normalized','Position',[0.30 0.22 0.45 0.60], ...
    'MenuBar','none','ToolBar','none','Color','w','CloseRequestFcn',@SF_MW_exit);

SF_S1_str = {'Functional images must be realigned, normalized and unsmoothed.'};
uicontrol(SF_MW,'Style','text','String',SF_S1_str,'Units','normalized', ...
    'Position',[0.03 0.92 0.95 0.05],'FontUnits','normalized','FontSize',0.48, ...
    'HorizontalAlignment','left','BackgroundColor','w');

% Parent folder with FUNC subfolders (by default)
func_root_default = fileparts(subject_paths{1});
if isempty(func_root_default), func_root_default = filesep; end
func_root = func_root_default;

% Button: Select parent folder (contains sub-*/ses-*/func)
SF_MW_S0 = uicontrol(SF_MW,'Style','pushbutton','String','Select parent folder (contains sub-*/ses-*/func)', ...
    'TooltipString','Folder that directly contains subject folders (e.g., sub-01, sub-02) with FUNC subfolders (which include unsmoothed functional images).', ...
    'Units','normalized','Position',[0.025 0.82 0.41 0.080],'FontUnits','normalized','FontSize',0.295, ...
    'Callback',@select_func_root);

% Box panel: Select parent folder (contains sub-*/ses-*/func)
SF_MW_S0_panel = uipanel(SF_MW,'Units','normalized','Position',[0.45 0.823 0.525 0.076], ...
    'HighlightColor',[0.78 0.78 0.78],'BorderType','line','BackgroundColor','w');

% Box text: Select parent folder (contains sub-*/ses-*/func)
SF_MW_S0_txt = uicontrol('Parent',SF_MW_S0_panel,'Style','text','String',func_root, ...
    'Units','normalized','Position',[0.03 0.05 0.94 0.64],'FontUnits','normalized','FontSize',0.48, ...
    'HorizontalAlignment','center','BackgroundColor','w','ForegroundColor',[0 0 0]);

% Button: Select the FUNC subfolder for the first subject
SF_MW_S1 = uicontrol(SF_MW,'Style','pushbutton','String','Select the FUNC subfolder for the first subject', ...
    'Units','normalized','Position',[0.025 0.73 0.41 0.080],'FontUnits','normalized','FontSize',0.295, ...
    'Callback',@select_subfolder);

% Box panel: Select the FUNC subfolder for the first subject
SF_MW_S1_panel = uipanel(SF_MW,'Units','normalized','Position',[0.45 0.733 0.525 0.076], ...
    'HighlightColor',[0.78 0.78 0.78],'BorderType','line','BackgroundColor','w');

% Box text: Select the FUNC subfolder for the first subject
SF_MW_S1_txt = uicontrol('Parent',SF_MW_S1_panel,'Style','text','String','Not selected', ...
    'ForegroundColor','red','Units','normalized','Position',[0.03 0.14 0.94 0.64], ...
    'FontUnits','normalized','FontSize',0.55,'HorizontalAlignment','center','BackgroundColor','w');

% Button: Text filter
SF_MW_B1 = uicontrol(SF_MW,'Style','pushbutton','String','Apply text filter unique for functional images:', ...
    'Units','normalized','Position',[0.025 0.64 0.41 0.080],'FontUnits','normalized','FontSize',0.295, ...
    'Callback',@apply_filter);

% Edit: Text filter
SF_MW_B1_E = uicontrol(SF_MW,'Style','edit','String','^war*.nii', ...
    'Units','normalized','Position',[0.45 0.64 0.525 0.080], ...
    'FontUnits','normalized','FontSize',0.32,'BackgroundColor','w');

% Preview list (first matched file per subject)
SF_MW_LB1 = uicontrol(SF_MW,'Style','listbox','String','','Max',100000, ...
    'Units','normalized','Position',[0.025 0.25 0.95 0.37],'FontUnits','points', ...
    'FontSize',10,'Value',[]);

% Button: Preserve functional images' paths from SPM.mat files
SF_MW_KEEP = uicontrol(SF_MW,'Style','pushbutton',...
    'String','Preserve functional image paths from the SPM.mat files',...
    'TooltipString','(i.e., if the GLMs were specified with unsmoothed functional images)',...
    'Units','normalized','Position',[0.025 0.15 0.95 0.078],...
    'FontUnits','normalized','FontSize',0.3,...
    'Callback',@preserve_from_spm_cb);

% Button: OK
SF_MW_OK   = uicontrol(SF_MW,'Style','pushbutton','String','OK', ...
    'Units','normalized','Position',[0.300 0.045 0.199 0.078], ...
    'FontUnits','normalized','FontSize',0.30,'Callback',@export_paths);

% Button: Help
SF_MW_HELP = uicontrol(SF_MW,'Style','pushbutton','String','Help', ...
    'Units','normalized','Position',[0.510 0.045 0.200 0.078], ...
    'FontUnits','normalized','FontSize',0.30,'Callback',@open_help);

movegui(SF_MW,'center');

% Close GUI
function SF_MW_exit(~,~)
    func_paths = [];
    fprintf(2,'Functional images are not selected.\n');
    uiresume(SF_MW);
end

% Select parent folder (contains sub-*/ses-*/func)
function select_func_root(~,~)
    tmp = deblank(spm_select(1,'dir','Select parent folder (contains sub-*/ses-*/func)',{},func_root,'..'));
    if ~isempty(tmp)
        user_changed_root = ~strcmp(tmp, func_root);
        func_root = tmp;
        set(SF_MW_S0_txt,'String',func_root,'ForegroundColor',[0 0 0],'HorizontalAlignment','center');
        % Reset downstream choices & cached mappings
        func_subfolder_full = '';
        func_subfolder_rel  = '';
        user_changed_subf    = false;
        set(SF_MW_S1_txt,'String','Not selected','ForegroundColor','red','HorizontalAlignment','center');
        set(SF_MW_LB1,'String','');  % clear preview list
        unsm_index = struct([]);     % clear mappings
    end
end

% Select the FUNC subfolder for the first subject
function select_subfolder(~,~)
    set(SF_MW_S1_txt,'String','Not selected','ForegroundColor','red','HorizontalAlignment','center');

    [~, subj_id1] = fileparts(subject_paths{1});
    start_dir = fullfile(func_root, subj_id1);
    if ~exist(start_dir,'dir'), start_dir = func_root; end
    if ~exist(start_dir,'dir'), start_dir = subject_paths{1}; end

    sel_dir = deblank(spm_select(1,'dir','Select ''functional'' subfolder',{},start_dir,'..'));
    if ~isempty(sel_dir)
        user_changed_subf  = true;
        func_subfolder_full = sel_dir;
        update_rel_display();
        set(SF_MW_LB1,'String','');  % clear preview (must re-apply filter)
        unsm_index = struct([]);     % clear mappings
    else
        fprintf(2,'Please select ''functional'' subfolder for the first subject.\n');
    end
end

% Compute relative path & show it (relative to <func_root>/<sub-01>)
function update_rel_display()
    [~, subj_id1] = fileparts(subject_paths{1});
    base_first = [fullfile(func_root, subj_id1) filesep];
    if strncmpi(func_subfolder_full, base_first, numel(base_first))
        func_subfolder_rel = func_subfolder_full(numel(base_first)+1:end);
    else
        base_first2 = [subject_paths{1} filesep];
        func_subfolder_rel = strrep(func_subfolder_full, base_first2, '');
    end
    if isempty(func_subfolder_rel)
        set(SF_MW_S1_txt,'String','Not selected','ForegroundColor','red','HorizontalAlignment','center');
    else
        set(SF_MW_S1_txt,'String',func_subfolder_rel,'ForegroundColor',[0 0 0], ...
            'HorizontalAlignment','center');
    end
end

% Map original SPM.mat directory -> target directory under FUNC root.
% If a FUNC subfolder was selected, we FORCE searches under that subfolder,
% keeping only the tail (e.g., "sess-01") from the original path.
function target_dir = map_dir_to_target(orig_dir, subj_id)
    % remainder of path after ".../sub-XX/"
    low_dir = lower(orig_dir);
    low_key = [filesep lower(subj_id) filesep];
    remainder = '';
    idx = strfind(low_dir, low_key);
    if ~isempty(idx)
        cut = idx(1) + length(low_key) - 1;
        remainder = orig_dir(cut+1:end); % e.g. "func/sess-01" or "func"
    end

    if ~isempty(func_subfolder_rel)
        % drop the first component of the remainder (e.g., "func") and keep the tail (e.g., "sess-01")
        tail = '';
        if ~isempty(remainder)
            parts = strsplit(remainder, filesep);
            if numel(parts) >= 2
                tail = fullfile(parts{2:end});
            end
        end
        % Force under the subfolder the user selected
        target_dir = fullfile(func_root, subj_id, func_subfolder_rel, tail);
    else
        % No explicit subfolder selected — preserve the whole remainder
        if isempty(remainder)
            target_dir = fullfile(func_root, subj_id);
        else
            target_dir = fullfile(func_root, subj_id, remainder);
        end
    end
end

% Build scan groups by original directory in SPM.xY.VY
function groups = build_scan_groups(SPM)
    n = numel(SPM.xY.VY);
    groups = struct('orig_dir',{},'scan_indices',{});
    seen = {};
    for k=1:n
        [d,~,~] = fileparts(SPM.xY.VY(k).fname);
        gi = find(strcmp(d, seen), 1, 'first');
        if isempty(gi)
            seen{end+1} = d;
            gidx = numel(groups)+1;
            groups(gidx).orig_dir    = d;
            groups(gidx).scan_indices = k;
        else
            gidx = gi;
            groups(gidx).scan_indices(end+1) = k;
        end
    end
end

% Apply text filter: find actual unsmoothed files per session/run dir
function apply_filter(~,~)
    f = msgbox('Selecting functional images. Please wait . . .');
    set(SF_MW_LB1,'String','');    % clear preview
    unsm_index = struct([]);       % reset mappings
    no_files = {};

    txt_filter = get(SF_MW_B1_E,'String');
    txt_filter = strrep(txt_filter,' ','');
    if isempty(txt_filter)
        fprintf(2,'Filter is empty or invalid, please re-enter.\n');
        try; close(f); end
        return;
    end
    % 'dir' uses wildcards, not regex; allow '^war*.nii' in UI by stripping '^'
    dir_filter = regexprep(txt_filter,'^\^','');

    preview_first = {};

    for iSub = 1:numel(subject_paths)
        [~, subj_id] = fileparts(subject_paths{iSub});
        SPM = load(SPM_paths{iSub}).SPM;

        % Build groups of scans in order, per original session directory
        groups = build_scan_groups(SPM);

        out_groups = struct('orig_dir',{},'target_dir',{},'unsm_files',{},'scan_indices',{});
        first_added = '';

        for g = 1:numel(groups)
            orig_dir   = groups(g).orig_dir;
            scan_idx   = groups(g).scan_indices;

            target_dir = map_dir_to_target(orig_dir, subj_id);
            cand = dir(fullfile(target_dir, dir_filter));
            
            % Only allow fallback to the ORIGINAL SPM folder if the user did NOT change root/subfolder.
            if isempty(cand) && (~user_changed_root && ~user_changed_subf)
                cand = dir(fullfile(orig_dir, dir_filter));
                if ~isempty(cand)
                    target_dir = orig_dir;
                end
            end

            % Only allow fallback to original dir when user didn't change root/subfolder
            if isempty(cand) && (~user_changed_root && ~user_changed_subf)
                cand = dir(fullfile(orig_dir, dir_filter));
                if ~isempty(cand)
                    target_dir = orig_dir;
                end
            end

            if ~isempty(cand)
                % sort by name (assumes zero-padded numbering; common for NIfTI series)
                names = {cand.name}';
                names = sort(names);
                unsm_files = cellfun(@(z) fullfile(target_dir,z), names, 'UniformOutput', false);

                out_groups(end+1).orig_dir    = orig_dir;        
                out_groups(end).target_dir    = target_dir;
                out_groups(end).unsm_files    = unsm_files;      % could be 1 (4D) or many (3D series)
                out_groups(end).scan_indices  = scan_idx;

                if isempty(first_added)
                    first_added = unsm_files{1};
                end
            end
        end

        if isempty(out_groups)
            no_files = vertcat(no_files, subject_paths(iSub));
        else
            unsm_index(iSub).groups = out_groups; 
            preview_first{end+1,1} = first_added; 
        end
    end

    if isempty(preview_first)
        fprintf(2,'No images matched the filter "%s". Please check the FUNC root/subfolder and filter.\n', txt_filter);
    end

    set(SF_MW_LB1,'String',preview_first);
    try; close(f); end
end

% Preserve FUNC paths form SPM.mat files
function preserve_from_spm_cb(~,~)
    % Set mode and preview: first image per subject from SPM.xY.VY(1).fname
    fpr1 = msgbox('Selecting functional images. Please wait . . .');
    preserve_spm_paths = true;
    spm_first_images = cell(numel(SPM_paths),1);
    for iSub = 1:numel(SPM_paths)
        SPM = load(SPM_paths{iSub}).SPM;
        spm_first_images{iSub,1} = SPM.xY.VY(1).fname;
    end
    try; close(fpr1); end
    set(SF_MW_LB1,'String',spm_first_images,'Value',[]);
end

% Check and export paths (per scan; 4D or 3D series handled)
function export_paths(~,~)
    if preserve_spm_paths
        fpr2 = msgbox('Selecting functional images. Please wait . . .');
        func_paths = struct([]);
        for jSub = 1:numel(SPM_paths)
            SPM = load(SPM_paths{jSub}).SPM;
            nScan = numel(SPM.xY.VY);
            func_paths(jSub).fname = cell(nScan,1);
            for kScan = 1:nScan
                v = SPM.xY.VY(kScan);
                func_paths(jSub).fname{kScan,1} = sprintf('%s,%d', v.fname, v.n(1));
            end
        end
        disp('Functional images selected (preserved from SPM.mat).');
        try; close(fpr2); end
        uiresume(SF_MW);
        return;
    end

    f3 = msgbox('Checking functional images. Please wait . . .');

    if isempty(unsm_index) || ~isfield(unsm_index, 'groups') || isempty([unsm_index.groups])
        fprintf(2,'Functional images are not selected, please apply the filter first.\n');
        try; close(f3); end
        return;
    end

    func_paths = struct([]);
    miss_sub = false(1,numel(SPM_paths));

    for jSub = 1:numel(SPM_paths)
        SPM = load(SPM_paths{jSub}).SPM;

        if jSub>numel(unsm_index) || ~isfield(unsm_index(jSub),'groups') || isempty(unsm_index(jSub).groups)
            miss_sub(jSub) = true;
            continue;
        end

        G = unsm_index(jSub).groups;

        % Initialize cell big enough for all scans
        nScan = numel(SPM.xY.VY);
        func_paths(jSub).fname = cell(nScan,1);

        ok_this_sub = true;

        for g = 1:numel(G)
            scan_idx = G(g).scan_indices(:)';         % indices of scans for this session
            unsmf    = G(g).unsm_files(:)';           % matched unsmoothed files for this session

            if isempty(scan_idx)
                continue;
            end

            if numel(unsmf) == 1
                % Assume 4D file; use volume index from original SPM
                file4D = unsmf{1};
                for kk = scan_idx
                    vol_idx = SPM.xY.VY(kk).n(1);
                    func_paths(jSub).fname{kk,1} = [file4D ',' num2str(vol_idx)];
                end
            else
                % Assume series of 3D files
                if numel(unsmf) < numel(scan_idx)
                    fprintf(2,['Not enough files matched in "%s" for subject %d: ' ...
                               'found %d, need %d. Check filter or folder mapping.\n'], ...
                               G(g).target_dir, jSub, numel(unsmf), numel(scan_idx));
                    ok_this_sub = false;
                    break;
                end
                % Map scan #i to file #i (1-based), append ",1" for SPM compatibility
                for i = 1:numel(scan_idx)
                    kk = scan_idx(i);
                    func_paths(jSub).fname{kk,1} = [unsmf{i} ',1'];
                end
            end
        end

        if ~ok_this_sub
            miss_sub(jSub) = true;
        end
    end

    if any(miss_sub)
        missing_idx = find(miss_sub);
        miss = subject_paths(missing_idx);
        missing_images_GUI(miss);
    else
        disp('Functional images selected.');
        uiresume(SF_MW);
    end

    try; close(f3); end
end

% Warning window: missing images
function missing_images_GUI(no_files_loc)
    SF_WW = figure('Name','Select subjects','NumberTitle','off','Units','normalized', ...
        'Position',[0.32 0.30 0.35 0.28],'Color','w','MenuBar','none','ToolBar','none','WindowStyle','Modal');
    uicontrol(SF_WW,'Style','listbox','String',no_files_loc,'Max',1000,'Units','normalized', ...
        'Position',[0.032 0.250 0.940 0.520],'FontUnits','points','FontSize',10,'Value',[]);
    uicontrol(SF_WW,'Style','text','String', ...
        'Warning, functional images are missing for the following subjects:', ...
        'Units','normalized','Position',[0.15 0.820 0.720 0.095], ...
        'FontUnits','points','FontSize',11,'HorizontalAlignment','center','BackgroundColor','w');
    uicontrol(SF_WW,'Style','pushbutton','String','OK','Units','normalized', ...
        'Position',[0.415 0.06 0.180 0.120],'FontUnits','normalized','FontSize',0.30, ...
        'Callback',@(s,e) close(SF_WW));
    movegui(SF_WW,'center');
    uiwait(SF_WW);
end

% Help window
function open_help(~,~)
   page1 = {
        ''
        '   First, select the parent folder that contains all subject folders' 
        '   with FUNC subfolders (if necessary).'
        ''
        '   Then select the FUNC subfolder for the first subject and apply text'
        '   filter (e.g., *war*.nii, *wr*.nii, or *preproc*.nii.gz) to match all fMRI images.'
        ''
        ''
        '   ================= EXAMPLE #1 (SPM-like folder structure) =================='
        ''
        '   There is no need to change the parent folder to select functional files.'
        ''
        '   project/'
        '   ├─ rawdata/      # DICOM'
        '   └─ derivatives/  <-------------------- [Parent folder with FUNC subfolders (BY DEFAULT)]'
        '      ├─ sub-01/    <-------------------------------------------- [Selected subject folder]'
        '      │  ├─ anat/'
        '      │  │  ├─ *T1*.nii'
        '      │  │  └─ *T1*.nii derivatives'
        '      │  ├─ func/  <--------------------- [Select the FUNC subfolder for the first subject] (1)'
        '      │  │  ├─ sess-01/'
        '      │  │  │  ├─ Unprocessed functional files (*.nii)'
        '      │  │  │  └─ Preprocessed functional files:'
        '      │  │  │       ◦ smoothed + norm. + real. (e.g., swar*.nii)'
        '      │  │  │       ◦ unsmoothed + norm. + real. (e.g., war*.nii) <---- [Apply text filter] (2)' 
        '      │  │  └─ sess-02/ ...'
        '      │  └─ stat/              # first-level models (one folder per GLM)'
        '      │     ├─ GLM-01/'
        '      │     │  ├─ SPM.mat  <--------------------------------------- [Selected SPM.mat file]'
        '      │     │  └─ TMFC_denoise/ <------------------------------------------ [Output folder]'
        '      │     └─ GLM-02/ ...'
        '      └─ sub-02/ ...'
    };

    page2 = {
        ''
        '   ================ EXAMPLE #2 (BIDS-like folder structure) ================='
        ''
        '   "project/derivatives/firstlevel-spm" parent folder with STAT subfolders needs to be'
        '   changed to "project/derivatives/fmriprep" — parent folder with FUNC subfolders.'
        ''
        '   project/'
        '   ├── sub-01/'
        '   │   ├── ses-01/'
        '   │   │   ├── anat/'
        '   │   │   │   └── *T1*.nii'
        '   │   │   └── func/     # Unprocessed functional files'
        '   │   └── ses-02/ ...'
        '   ├── sub-02/ ...'
        '   └── derivatives'
        '       ├── fmriprep/<------------------ [Select parent folder (contains sub-*/ses-*/func)] (1)'
        '       │   ├── sub-01/'
        '       │   │   ├── ses-01/'
        '       │   │   │   └── func/<----------- [Select the FUNC subfolder for the first subject] (2)'
        '       │   │   │       └── Preprocessed functional files:'
        '       │   │   │           ◦ smoothed + normalized + realigned'
        '       │   │   │           ◦ unsmoothed + normalized + realigned <---- [Apply text filter] (3)'
        '       │   │   └── ses-02/ ...'
        '       │   └── sub-02/ ...'
        '       └── firstlevel-spm/  <---- [Parent folder with FUNC (BY DEFAULT)](Needs to be changed!)'
        '           ├── sub-01/     <------------------------------------ [Selected subject folder]'
        '           │   ├── GLM-01/'
        '           │   │   ├── SPM.mat        <--------------------------- [Selected SPM.mat file]'
        '           │   │   └── TMFC_denoise/  <----------------------------------- [Output folder]'
        '           │   └── GLM-02/ ...'
        '           └── sub-02/ ...'
    };

    page3 = {
        ''
        '   ================ EXAMPLE #3 (Other non-BIDS folder structure) ================'
        ''
        '   "project/firstlevel-spm" parent folder with STAT subfolders needs'
        '   to be changed to "project/nifti"  — parent folder with FUNC subfolders.'
        ''
        '   project/'
        '   ├─ rawdata/ # DICOM'
        '   ├─ nifti/   <-------------- [Select parent folder (contains sub-*/ses-*/func)]  (1)'
        '   │  ├─ sub-01/'
        '   │  │  ├─ anat/'
        '   │  │  │  ├─ *T1*.nii  '
        '   │  │  │  └─ *T1*.nii derivatives'
        '   │  │  └─ func/   <---------- [Select the FUNC subfolder for the first subject]  (2)'
        '   │  │     ├─ sess-01/'
        '   │  │     │  ├─ Unprocessed functional files (*.nii)'
        '   │  │     │  └─ Preprocessed functional files (*.nii):'
        '   │  │     │     ◦ smoothed + normalized + realigned'
        '   │  │     │     ◦ unsmoothed + normalized + realigned <---- [Apply text filter] (3)'
        '   │  │     └─ sess-02/ ...'
        '   │  └─ sub-02/ ...'
        '   └─ firstlevel-spm/  <-- [Parent folder with FUNC subfolders (BY DEFAULT)](Needs to be changed!)'
        '      ├─ sub-01/   <----------------------------------- [Selected subject folder]'
        '      │  ├─ GLM-01/'
        '      │  │  ├─ SPM.mat    <------------------------------ [Selected SPM.mat file]'
        '      │  │  └─ TMFC_denoise/   <--------------------------------- [Output folder]'
        '      │  └─ GLM-02/ ...'
        '      └─ sub-02/ ...'
    };

    pages = {page1,page2,page3};
    cur=1; total=numel(pages);

    ST_HELP = figure('Name','Help','NumberTitle','off','Units','normalized', ...
        'Position',[0.22 0.10 0.60 0.75], 'MenuBar','none','ToolBar','none', ...
        'Color','w','WindowStyle','Modal');

    txtArea = uicontrol(ST_HELP,'Style','edit', ...
        'Units','normalized','Position',[0.05 0.14 0.90 0.78], ...
        'BackgroundColor','w','Enable','inactive', ...
        'Max', 2, 'Min', 0, 'HorizontalAlignment','left', ...
        'FontName','Courier New','FontUnits','normalized','FontSize',0.025);

    btnPrev = uicontrol(ST_HELP,'Style','pushbutton','String','Previous', ...
        'Units','normalized','Position',[0.05 0.05 0.14 0.06], ...
        'FontUnits','normalized','FontSize',0.35, 'Callback',@go_prev);

    pageLbl = uicontrol(ST_HELP,'Style','text','String','', ...
        'Units','normalized','Position',[0.205 0.035 0.19 0.06], ...
        'BackgroundColor','w','HorizontalAlignment','center', ...
        'FontUnits','normalized','FontSize',0.40);

    btnNext = uicontrol(ST_HELP,'Style','pushbutton','String','Next', ...
        'Units','normalized','Position',[0.405 0.05 0.14 0.06], ...
        'FontUnits','normalized','FontSize',0.35, 'Callback',@go_next);

    btnOK = uicontrol(ST_HELP,'Style','pushbutton','String','OK', ...
        'Units','normalized','Position',[0.81 0.05 0.14 0.06], ...
        'FontUnits','normalized','FontSize',0.35, ...
        'Callback',@(s,e) close(ST_HELP));

    render_page(); movegui(ST_HELP,'center'); uiwait(ST_HELP);

    function render_page()
        set(txtArea,'String', pages{cur});
        set(pageLbl,'String', sprintf('Page %d of %d', cur, total));
        set(btnPrev,'Enable', tern(cur>1,'on','off'));
        set(btnNext,'Enable', tern(cur<total,'on','off'));
    end
    function go_prev(~,~), if cur>1, cur=cur-1; render_page(); end, end
    function go_next(~,~), if cur<total, cur=cur+1; render_page(); end, end
    function out = tern(cond,a,b), if cond, out=a; else, out=b; end, end
end

% -------------------------------------------------------------------------
uiwait(SF_MW);
delete(SF_MW);
end
