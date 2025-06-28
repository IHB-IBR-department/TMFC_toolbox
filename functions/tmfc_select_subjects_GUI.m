 function [SPM_paths, subject_paths, SPM_subfolder] = tmfc_select_subjects_GUI(SPM_check)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Opens a GUI window for selecting individual subject SPM.mat files
% created by SPM12 after 1-st level GLM estimation. Optionally checks
% SPM.mat files: 
% (1) checks if all SPM.mat files are present in the specified paths
% (2) checks if the same conditions are specified in all SPM.mat files
% (3) checks if output folders specified in SPM.mat files exist
% (4) checks if functional files specified in SPM.mat files exist
%
% FORMAT [paths] = tmfc_select_subjects_GUI(SPM_check)
%
%   Inputs: 
%   SPM_check         - 0 or 1 (don't check or check SPM.mat files)
%
%   Outputs:
%   SPM_paths         - Full paths to selected SPM.mat files
%   subject_paths     - Paths to selected subjects
%   SPM_subfolder     - Subfolder for SPM.mat files
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

if nargin == 0
	SPM_check = 1;
end
                    
% SS = select subjects, MW = main window 
SS_MW = figure('Name', 'Select subjects', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.36 0.25 0.35 0.575],'MenuBar', 'none','ToolBar', 'none','color','w','CloseRequestFcn',@SS_MW_exit);
SS_MW_S1 = uicontrol(SS_MW,'Style','text','String', 'Not Selected','ForegroundColor','red','Units', 'normalized', 'Position',[0.500 0.820 0.450 0.095],'backgroundcolor','w','FontUnits','normalized','FontSize',0.25);
SS_MW_S2 = uicontrol(SS_MW,'Style','text','String', 'Not Selected','ForegroundColor','red','Units', 'normalized', 'Position',[0.500 0.720 0.450 0.095],'backgroundcolor','w','FontUnits','normalized','FontSize',0.25);
SS_MW_LB1 = uicontrol(SS_MW, 'Style', 'listbox', 'String', '','Max', 100000,'Units', 'normalized', 'Position',[0.033 0.250 0.920 0.490],'FontUnits','points','FontSize',10,'Value', [],'callback', @SS_LB_select);
SS_MW_sel_sub = uicontrol(SS_MW,'Style','pushbutton', 'String', 'Select subject folders','Units', 'normalized', 'Position',[0.033 0.850 0.455 0.095],'FontUnits','normalized','FontSize',0.25,'callback', @select_sub);
SS_MW_sel_mat = uicontrol(SS_MW,'Style','pushbutton', 'String', 'Select SPM.mat file for Subject #1','Units', 'normalized', 'Position',[0.033 0.750 0.455 0.095],'FontUnits','normalized','FontSize',0.25,'callback', @select_SPM_mat);
SS_MW_add_new = uicontrol(SS_MW,'Style','pushbutton', 'String', 'Add new subject','Units', 'normalized', 'Position',[0.033 0.14 0.300 0.095],'FontUnits','normalized','FontSize',0.25,'callback', @add_new);
SS_MW_rem = uicontrol(SS_MW,'Style','pushbutton', 'String', 'Remove selected subject','Units', 'normalized', 'Position',[0.346 0.14 0.300 0.095],'FontUnits','normalized','FontSize',0.25,'callback', @remove_sub);
SS_MW_rem_all = uicontrol(SS_MW,'Style','pushbutton', 'String', 'Clear all subjects','Units', 'normalized', 'Position',[0.660 0.14 0.300 0.095],'FontUnits','normalized','FontSize',0.25,'callback', @remove_all);
SS_MW_conf = uicontrol(SS_MW,'Style','pushbutton', 'String', 'OK','Units', 'normalized', 'Position',[0.390 0.04 0.200 0.080],'FontUnits','normalized','FontSize',0.28,'callback', @confirm_paths);
movegui(SS_MW,'center');

% Temporary variables
subject_paths_tmp = {};       % Variable to store subject paths
subject_full_path = {};       % Variable to store full paths
SPM_mat_path = {};            % Varaible to store subfolder for SPM.mat file
selected_sub = {};            % Variable to store the selected list of paths (as INDEX)
add_new_subs = {};            % Variable used to add new subjects

%--------------------------------------------------------------------------
% Select subjects from the list
function SS_LB_select(~,~)
    index = get(SS_MW_LB1, 'Value');     
	selected_sub = index;                
end

%--------------------------------------------------------------------------
% Select subjects
function select_sub(~,~)
    
	set(SS_MW_LB1, 'String', '');              
    subject_paths_tmp = add_subjects();             

    if isempty(subject_paths_tmp)
    	disp('TMFC Subjects: 0 Subjects selected');
        set(SS_MW_S1,'String','Not selected','ForegroundColor','red');
        set(SS_MW_S2,'String','Not selected','ForegroundColor','red');
        SPM_mat_path = {};
        subject_full_path = {};
    else
        fprintf('TMFC Subjects: Subjects selected are: %d \n', size(subject_paths_tmp,1));
        disp('TMFC Subjects: Select SPM.mat file for the first subject.');
        set(SS_MW_S1,'String', strcat(num2str(size(subject_paths_tmp,1)),' selected'),'ForegroundColor',[0.219, 0.341, 0.137]);    
        set(SS_MW_S2,'String','Not selected','ForegroundColor','red');
        SPM_mat_path = {};
        subject_full_path = {};
    end
    
    if isempty(subject_full_path) && isempty(SPM_mat_path) && isempty(subject_paths_tmp)
        fprintf(2,'No subjects selected.\n');
        SPM_mat_path = {}; 
        subject_full_path = {}; 
        subject_paths_tmp = {}; 
        set(SS_MW_S1,'String','Not selected','ForegroundColor','red');
        set(SS_MW_S2,'String','Not selected','ForegroundColor','red');
    end 
end

%--------------------------------------------------------------------------
% Select SPM.mat file
function select_SPM_mat(~,~)
    if isempty(subject_paths_tmp)
        fprintf(2,'TMFC Subjects: Please select subject folders.\n');        
        
    elseif isempty(subject_full_path) && isempty(subject_paths_tmp)
        fprintf(2,'TMFC Subjects: Please select subject folders.\n');
        set(SS_MW_S2,'String','Not selected','ForegroundColor','red');
        set(SS_MW_LB1, 'String', '');
        
    else
        [subject_full_path, SPM_mat_path] = add_mat_file(subject_paths_tmp);            
        if ~isempty(SPM_mat_path)
            set(SS_MW_LB1, 'String', subject_full_path);                          
            disp('TMFC Subjects: The SPM.mat file has been succesfully selected.');
            set(SS_MW_S2,'String','Selected','ForegroundColor',[0.219, 0.341, 0.137]);
            SPM_subfolder = strrep(SPM_mat_path, '\SPM.mat', '');
        else
            fprintf(2,'TMFC Subjects: The SPM.mat file has not been selected.\n');
            set(SS_MW_S2,'String','Not selected','ForegroundColor','red');
            set(SS_MW_LB1,'String','');
        end 
    end
end

%--------------------------------------------------------------------------
% Add new subjects to the list
function add_new(~,~)   
    if isempty(subject_paths_tmp) || isempty(subject_full_path)
        fprintf(2,'TMFC Subjects: No existing list of subjects present. Please select subjects via ''Select subject folders'' button. \n');
        
    elseif isempty(SPM_mat_path)
        fprintf(2,'TMFC Subjects: Cannot add new subjects without SPM.mat file. Please select subjects via ''Select subject folders'' button and proceed to Select SPM.mat file. \n');

    else
        add_subs_full_path = {};                    
        add_new_subs = add_subjects();             
        
        if isempty(add_new_subs)
            fprintf(2,'TMFC Subjects: No newly selected subjects. \n');
        else
            subs_exist = size(subject_full_path,1); 
            
            for iSub = 1:size(add_new_subs,1)
               add_subs_full_path =  vertcat(add_subs_full_path,strcat(char(add_new_subs(iSub,:)),char(SPM_mat_path)));
            end
                        
            subject_full_path = vertcat(subject_full_path, add_subs_full_path);   % Joining exisiting list of subjects with new subjects
            new_subs_count = size(unique(subject_full_path)) - subs_exist;        % Removing duplicates
            subject_full_path = unique(subject_full_path);
            
            subject_paths_tmp = vertcat(subject_paths_tmp, add_new_subs);
            subject_paths_tmp = unique(subject_paths_tmp);

            if new_subs_count(1) == 0
                fprintf(2,'TMFC Subjects: Newly selected subjects are already present in the list, no new subjects added. \n');
            else
                fprintf('TMFC Subjects: New subjects selected: %d. \n', new_subs_count(1)); 
            end   
        end 
        
        set(SS_MW_LB1, 'String', subject_full_path);                              % Updating display with new subjects
        set(SS_MW_S1,'String', strcat(num2str(size(subject_full_path,1)),' selected'),'ForegroundColor',[0.219, 0.341, 0.137]);
        clear add_subs_full_path new_subs_count add_new_subs
    end
end

%--------------------------------------------------------------------------
% Remove subjects from the list
function remove_sub(~,~)    
    if isempty(selected_sub)
        fprintf(2,'TMFC Subjects: There are no selected subjects to remove from the list. Please select subjects to remove.\n');
    else
        subject_full_path(selected_sub,:) = [];   
        subject_paths_tmp(selected_sub,:) = [];   
        fprintf('TMFC Subjects: Number of subjects removed: %d. \n', size(selected_sub,2));
        
        set(SS_MW_LB1,'Value',[]);                                             
        set(SS_MW_LB1, 'String', subject_full_path);                              
        selected_sub ={};
        
        if size(subject_full_path,1) < 1
            set(SS_MW_S1,'String', 'Not selected','ForegroundColor','red');
        else
            set(SS_MW_S1,'String', strcat(num2str(size(subject_full_path,1)),' selected'),'ForegroundColor',[0.219, 0.341, 0.137])
        end
    end
end 

%--------------------------------------------------------------------------
% Remove all subjects
function remove_all(~,~)
    if isempty(subject_paths_tmp) || isempty(subject_full_path)
        fprintf(2,'TMFC Subjects: No subjects present to remove.\n');
    else
        subject_paths_tmp = {};
        subject_full_path = {};
        SPM_mat_path = {};
        add_new_subs = {};        
        disp('TMFC Subjects: All subjects have been removed.');
        set(SS_MW_LB1, 'String', '');                
        set(SS_MW_S1,'String', 'None selected','ForegroundColor','red');
        set(SS_MW_S2,'String', 'None selected','ForegroundColor','red');
    end 
end

%--------------------------------------------------------------------------
% Check SPM.mat files and export paths
function confirm_paths(~,~)
	file_correct = {};
    file_exist = {};
    file_dir = {};
    file_func = {};
         
    % Initial checks
    if isempty(subject_paths_tmp)
        fprintf(2,'TMFC Subjects: There are no selected subjects. Please select subjects and SPM.mat files.\n');
    elseif (isempty(subject_full_path) && isempty(SPM_mat_path)) || (~isempty(subject_paths_tmp) && isempty(SPM_mat_path))
        fprintf(2,'TMFC Subjects: Please select SPM.mat file for the first subject.\n');
    elseif (isempty(subject_full_path) && ~isempty(SPM_mat_path))
        fprintf(2,'TMFC Subjects: Please re-select subjects and SPM.mat file if required.\n');
    else
        SS_MW_exit(SS_MW);   
        
        % Check SPM.mat files
        if SPM_check == 1              
            % Stage 1 - Check SPM.mat files existence
            [file_exist,subject_file_exist] = check_file_exist(subject_full_path,subject_paths_tmp);        
            if size(file_exist,1) == 0
            	fprintf(2,'TMFC Subjects: (Stage 1 Check Failed) - Selected SPM.mat files are missing from the directories. Please try again.\n');
                reset_paths();
            else                
                % Stage 2 - Check task conditions
                [file_correct,subject_file_correct] = check_file_cond(file_exist,subject_file_exist);          
                if size(file_correct,1) == 0
                    fprintf(2,'TMFC Subjects: (Stage 2 Check Failed) - Selected SPM.mat files have different task conditions and/or number of sessions. Please select SPM.mat files with the same task conditions and number of sessions.\n');
                    reset_paths();
                else                    
                    % Stage 3 - Check output directories 
                    [file_dir,subject_file_dir] = check_file_dir(file_correct,subject_file_correct);
                    if size(file_dir,1) == 0
                        fprintf(2,'TMFC Subjects: (Stage 3 Check Failed) - Directory where the output files will be saved are missing (check SPM.swd). Please select correct SPM.mat files or change paths in SPM.mat files.\n');
                        reset_paths();
                    else                        
                        % Stage 4 - Check functional files
                        [file_func,subject_file_func] = check_file_func(file_dir,subject_file_dir);                
                        if size(file_func,1) == 0
                            fprintf(2,'TMFC Subjects: (Stage 4 Check Failed) - Functional files specified in SPM.mat file are missing. Please select correct SPM.mat files or change paths in SPM.mat files.\n');
                            reset_paths();
                        else
                            SPM_paths = file_func; 
                            subject_paths = subject_file_func;
                        end
                    end 
                end
            end
        else
            SPM_paths = subject_full_path; 
            subject_paths = subject_paths_tmp;
        end
    end                                                                 
end      

%--------------------------------------------------------------------------
% Clear temporary variables
function reset_paths(~,~)
	subject_paths_tmp = {};
    subject_full_path = {};
    SPM_mat_path = {};
    add_new_subs = {};
    selected_sub = {};
    SPM_subfolder = {};
end

%--------------------------------------------------------------------------
% Close select subjects GUI window
function SS_MW_exit(~,~) 
    if exist('paths', 'var') == 0
        SPM_paths = [];
        subject_paths = [];
        SPM_subfolder = [];
    end
    uiresume(SS_MW);
end

uiwait(SS_MW);
delete(SS_MW)
return;

end

%% Select subjects
function subject_dir = add_subjects(~,~)
    subjects = spm_select(inf,'dir','Select subject folders',{},pwd,'..');
    subject_dir = {};                % Cell to store subjects    
    % Updating list of Subjects
    for iSub = 1:size(subjects,1)
    	subject_dir = vertcat(subject_dir, deblank(subjects(iSub,:)));
    end
    subject_dir = unique(subject_dir);
end              

%% Select SPM.mat file
function [subject_full_path, SPM_mat_path] = add_mat_file(subject_dir)
    subject_full_path = {};  SPM_mat_path = {};
    [mat_file_path] = spm_select( 1,'any','Select SPM.mat file for the first subject',{}, strtrim(subject_dir(1,:)), 'SPM.*');
    if ~isempty(mat_file_path)
        [SPM_mat_path] = strrep(mat_file_path, strtrim(subject_dir(1,:)),'');   
    end
        
    for iSub = 1:size(subject_dir,1)
    	subject_full_path =  vertcat(subject_full_path,strcat(char(subject_dir(iSub,:)),char(SPM_mat_path)));
    end
end 

%% Check SPM.mat files existence
function [file_exist,subject_file_exist] = check_file_exist(subject_full_path,subject_paths_tmp)

    file_exist = {};
    file_not_exist = {};   
    subject_file_exist = {};
    
    % Check if SPM.mat files exist
    for iSub = 1:length(subject_full_path) 
        if exist(subject_full_path{iSub}, 'file')
            file_exist{iSub,1} = subject_full_path{iSub};
            subject_file_exist{iSub,1} = subject_paths_tmp{iSub};
        else
            file_not_exist{iSub,1} = subject_full_path{iSub};
        end                
    end 
    
    % Check if the variables storing the existing files are empty
    try
        file_exist = file_exist(~cellfun('isempty', file_exist)); 
        subject_file_exist = subject_file_exist(~cellfun('isempty', subject_file_exist));
    end

    try 
        file_not_exist = file_not_exist(~cellfun('isempty',file_not_exist)); 
    end
    
    % Show missing SPM.mat files
    if length(file_exist) ~= length(subject_full_path)      
    	% SS_WW = select subjects warning window
        SS_WW1 = figure('Name', 'Select subjects', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.32 0.26 0.35 0.28], 'color', 'w','MenuBar', 'none','ToolBar', 'none');
        SS_WW_LB = uicontrol(SS_WW1, 'Style', 'listbox', 'String', file_not_exist,'Max',inf,'Units', 'normalized', 'Position',[0.032 0.250 0.940 0.520],'FontUnits','points','FontSize',10);
        SS_WW_S1 = uicontrol(SS_WW1,'Style','text','String', 'Warning, the following SPM.mat files are missing:','Units', 'normalized', 'Position',[0.15 0.820 0.720 0.095], 'FontUnits','normalized','FontSize',0.5,'backgroundcolor', 'w');
        SS_WW_close = uicontrol(SS_WW1,'Style','pushbutton', 'String', 'OK','Units', 'normalized',  'Position',[0.415 0.06 0.180 0.120] ,'FontUnits','normalized','FontSize',0.30,'callback', @close_SS_WW);
        movegui(SS_WW1,'center');
        uiwait(SS_WW1);
    end
    
    function close_SS_WW(~,~)
        close(SS_WW1);
    end   
end

%% Check conditions specified in SPM.mat files
function [file_correct,subject_file_correct] = check_file_cond(file_exist,subject_file_exist)
    
    file_incorrect = {};
    file_correct = {};
    subject_file_correct = {};
    
    if length(file_exist) > 1

        w = waitbar(0,'Check conditions','Name','Check SPM.mat files');

        % Reference SPM.mat file
        file_correct{1,1} = file_exist{1};
        subject_file_correct{1,1} = subject_file_exist{1};
        SPM_ref = load(file_exist{1});

        % Reference structure for conditions
        for iSess = 1:length(SPM_ref.SPM.Sess)
            cond_ref(iSess).sess = struct('name', {SPM_ref.SPM.Sess(iSess).U(:).name});
        end

        % Start check
        for iSub = 2:length(file_exist)
            
            % SPM.mat file to check
            SPM = load(file_exist{iSub}).SPM;

            % Structure for conditions to check
            for jSess = 1:length(SPM.Sess)
                cond(jSess).sess = struct('name', {SPM.Sess(jSess).U(:).name});
            end 

            if ~isequaln(cond_ref, cond)
                file_incorrect{iSub,1} = file_exist{iSub};
            else
                file_correct{iSub,1} = file_exist{iSub};
                subject_file_correct{iSub,1} = subject_file_exist{iSub};
            end

            try
                waitbar(iSub/length(file_exist),w);
            end

            clear SPM cond
        end
        
    else
        file_correct = file_exist;
        subject_file_correct = subject_file_exist;
    end

    try
        close(w)
    end

    try
        file_correct = file_correct(~cellfun('isempty', file_correct));
        subject_file_correct = subject_file_correct(~cellfun('isempty', subject_file_correct));
    end

    try
        file_incorrect = file_incorrect(~cellfun('isempty', file_incorrect));
    end

    % Show incorrect SPM.mat files
    if length(file_correct) ~= length(file_exist)
        if isunix; fontscale = 0.85; else; fontscale = 1; end
        SS_WW2 = figure('Name', 'Select subjects', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.32 0.30 0.35 0.28], 'color', 'w','MenuBar', 'none','ToolBar', 'none');
        SS_WW_LB = uicontrol(SS_WW2, 'Style', 'listbox', 'String', file_incorrect,'Max',inf,'Units', 'normalized', 'Position',[0.032 0.250 0.940 0.520],'FontUnits','points','FontSize',10);
        SS_WW_S1 = uicontrol(SS_WW2,'Style','text','String', 'Warning, the following SPM.mat files have different conditions specified:','Units', 'normalized', 'Position',[0.15 0.820 0.720 0.095], 'FontUnits','normalized','FontSize',0.5*fontscale,'backgroundcolor', 'w');
        SS_WW_close = uicontrol(SS_WW2,'Style','pushbutton', 'String', 'OK','Units', 'normalized',  'Position',[0.415 0.06 0.180 0.120] ,'FontUnits','normalized','FontSize',0.30,'callback', @close_SS_WW);
        movegui(SS_WW2,'center');
        uiwait(SS_WW2);
    end
    
    function close_SS_WW(~,~)
        close(SS_WW2);
    end
end


%% Check output directories (SPM.swd) specified in SPM.mat files
function [file_dir,subject_file_dir] = check_file_dir(file_correct,subject_file_correct)

	file_dir = {};
    file_no_dir = {};
    subject_file_dir = {};

    w = waitbar(0,'Check directories','Name','Check SPM.mat files');

    for iSub = 1:length(file_correct)
        
        SPM = load(file_correct{iSub}).SPM;

        if ~isfield(SPM, 'swd')
            error('SPM.swd field does not exist. Check SPM.mat files. Try to estimate GLM.');
        end

        if exist(SPM.swd, 'dir') 
            file_dir{iSub,1} = file_correct{iSub};
            subject_file_dir{iSub,1} = subject_file_correct{iSub};
        else
            file_no_dir{iSub,1} = file_correct{iSub};
        end

        
        clear SPM
        
        try
            waitbar(iSub/length(file_correct),w);
        end
    end
  
    try
        close(w)
    end

    try
        file_dir = file_dir(~cellfun('isempty', file_dir));
        subject_file_dir = subject_file_dir(~cellfun('isempty', subject_file_dir));
    end
    
    try 
        file_no_dir = file_no_dir(~cellfun('isempty',file_no_dir)); 
    end
    
    % Show incorrect SPM.mat files
    if length(file_dir) ~= length(file_correct)
        if isunix; fontscale = 0.9; else; fontscale = 1; end
        SS_WW3 = figure('Name', 'Select subjects', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.32 0.26 0.35 0.28], 'color', 'w','MenuBar', 'none','ToolBar', 'none');
        SS_WW_LB = uicontrol(SS_WW3, 'Style', 'listbox', 'String', file_no_dir,'Max',inf,'Units', 'normalized', 'Position',[0.032 0.250 0.940 0.520],'FontUnits','points','FontSize',10);
        SS_WW_S1 = uicontrol(SS_WW3,'Style','text','String', 'Warning, the output folder (SPM.swd) specified in the following SPM.mat files do not exist:','Units', 'normalized', 'Position',[0.08 0.820 0.840 0.080], 'FontUnits','normalized','FontSize',0.56*fontscale,'backgroundcolor', 'w');
        SS_WW_close = uicontrol(SS_WW3,'Style','pushbutton', 'String', 'OK','Units', 'normalized',  'Position',[0.415 0.06 0.180 0.120] ,'FontUnits','normalized','FontSize',0.30,'callback', @close_SS_WW);
        movegui(SS_WW3,'center');
    	uiwait(SS_WW3);
    end
    
    function close_SS_WW(~,~)
    	close(SS_WW3);
    end   
end
            
%% Check functional files specified in SPM.mat files 
function [file_func,subject_file_func] = check_file_func(file_dir,subject_file_dir) 

    file_func = {};
    file_no_func = {};
    subject_file_func = {};

    w = waitbar(0,'Check functional files','Name','Check SPM.mat files');

    for iSub = 1:length(file_dir) 
        
        SPM = load(file_dir{iSub}).SPM;
        
        for jImage = 1:length(SPM.xY.VY)
            funct_check(jImage) = exist(SPM.xY.VY(jImage).fname, 'file');
        end
        
        if nnz(funct_check) == length(SPM.xY.VY)
            file_func{iSub,1} = file_dir{iSub};
            subject_file_func{iSub,1} = subject_file_dir{iSub};
        else
            file_no_func{iSub,1} = file_dir{iSub};
        end
        clear SPM funct_check      
        
        try
            waitbar(iSub/length(file_dir),w);
        end
        
    end

    try
        close(w)
    end
    
    try
        file_func = file_func(~cellfun('isempty', file_func));
        subject_file_func = subject_file_func(~cellfun('isempty', subject_file_func));
    end
    
    try 
        file_no_func = file_no_func(~cellfun('isempty',file_no_func)); 
    end
    
    if length(file_func) ~= length(file_dir)
        if isunix; fontscale = 0.85; else; fontscale = 1; end
        SS_WW4 = figure('Name', 'Select subjects', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.32 0.26 0.35 0.28], 'color', 'w','MenuBar', 'none','ToolBar', 'none');
        SS_WW_LB = uicontrol(SS_WW4, 'Style', 'listbox', 'String', file_no_func,'Max',inf,'Units', 'normalized', 'Position',[0.032 0.250 0.940 0.520],'FontUnits','points','FontSize',10);
        SS_WW_S1 = uicontrol(SS_WW4,'Style','text','String', 'Warning, the functional files specified in the following SPM.mat files do not exist:','Units', 'normalized', 'Position',[0.15 0.820 0.750 0.095], 'FontUnits','normalized','FontSize',0.5*fontscale,'backgroundcolor', 'w');
        SS_WW_close = uicontrol(SS_WW4,'Style','pushbutton', 'String', 'OK','Units', 'normalized',  'Position',[0.415 0.06 0.180 0.120] ,'FontUnits','normalized','FontSize',0.30,'callback', @close_SS_WW);
        movegui(SS_WW4,'center');
        uiwait(SS_WW4);
    end
    
    function close_SS_WW(~,~)
        close(SS_WW4);
    end
end

