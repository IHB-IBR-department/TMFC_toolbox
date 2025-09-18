function tmfc_statistics_GUI()

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Configure and run group-level tests on TMFC connectivity matrices:
%   • One-sample t-test
%   • Paired t-test
%   • Two-sample t-test 
%
% Expected input: user selects *.mat files; each file must contain exactly
% one variable that is either:
%   • ROI×ROI (one subject per file), or
%   • ROI×ROI×Subjects (multiple subjects in one file).
%
% Note: Permutation tests and Network-Based Statistics (NBS) 
% are planned for a future release.
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

% GUI Initialization
stats_GUI = figure('Name', 'TMFC statistics', 'NumberTitle', 'off', 'Units', 'normalized', ...
    'Position', [0.38 0.21 0.22 0.56],'MenuBar', 'none','ToolBar', 'none','color','w', 'Tag','TMFC_STATS_GUI', ...
    'CloseRequestFcn', @on_close);

movegui(stats_GUI, 'center');
RES_T1  = uicontrol(stats_GUI,'Style','text','String','TMFC statistics','Units', 'normalized', 'Position',[0.270 0.93 0.460 0.05],'fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w', 'FontWeight', 'bold');

% Pop-up menu to select type of test
ST_test_type  = uicontrol(stats_GUI,'Style','popupmenu','String', {'One-sample t-test', 'Paired t-test', 'Two-sample t-test'},'Units', 'normalized', 'Position',[0.045 0.87 0.91 0.05],'fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w');
 
% Listboxes to display *.mat file selection
ST_lst_0 = uicontrol(stats_GUI , 'Style', 'listbox', 'String', '','Max', 100,'value',[],'Units', 'normalized', 'Position',[0.045 0.56 0.91 0.300],'fontunits','points', 'fontSize',11);
ST_lst_1 = uicontrol(stats_GUI , 'Style', 'listbox', 'String', '','Max', 100,'value',[],'Units', 'normalized', 'Position',[0.045 0.56 0.440 0.300],'fontunits','points', 'fontSize', 11,'visible','off');
ST_lst_2 = uicontrol(stats_GUI , 'Style', 'listbox', 'String', '','Max', 100,'value',[],'Units', 'normalized', 'Position',[0.52 0.56 0.440 0.300],'fontunits','points', 'fontSize',11,'visible','off');

% Counter of ROIs × subjects for selected files
ST_L0_CTR = uicontrol(stats_GUI, 'Style', 'text','String', '0 ROIs x 0 subjects','Units', 'normalized', 'Position',[0.295 0.51 0.44 0.04],'fontunits','normalized', 'fontSize', 0.57, 'HorizontalAlignment','center','backgroundcolor','w','ForegroundColor',[0.773, 0.353, 0.067]);
ST_L1_CTR = uicontrol(stats_GUI, 'Style', 'text','String', '0 ROIs x 0 subjects','Units', 'normalized', 'Position',[0.045 0.51 0.44 0.04],'fontunits','normalized', 'fontSize', 0.57, 'HorizontalAlignment','center','backgroundcolor','w','ForegroundColor',[0.773, 0.353, 0.067],'visible', 'off');
ST_L2_CTR = uicontrol(stats_GUI, 'Style', 'text','String', '0 ROIs x 0 subjects','Units', 'normalized', 'Position',[0.52 0.51 0.44 0.04],'fontunits','normalized', 'fontSize', 0.57, 'HorizontalAlignment','center','backgroundcolor','w','ForegroundColor',[0.773, 0.353, 0.067],'visible', 'off');

% "Select & Remove" file buttons for each case
ST_L0_SEL = uicontrol(stats_GUI,'Style','pushbutton','String', 'Select','Units', 'normalized','Position',[0.045 0.45 0.445 0.054],'fontunits','normalized', 'fontSize', 0.36, 'UserData', struct('select','one_samp_sel'));
ST_L0_REM = uicontrol(stats_GUI,'Style','pushbutton','String', 'Remove','Units', 'normalized','Position',[0.52 0.45 0.445 0.054],'fontunits','normalized', 'fontSize', 0.36, 'UserData', struct('remove','one_samp_rem'));
ST_L1_SEL = uicontrol(stats_GUI,'Style','pushbutton','String', 'Select','Units', 'normalized','Position',[0.045 0.45 0.210 0.054],'fontunits','normalized', 'fontSize', 0.36, 'visible', 'off','UserData',struct('select', 'left_samp_sel'));
ST_L1_REM = uicontrol(stats_GUI,'Style','pushbutton','String', 'Remove','Units', 'normalized','Position',[0.275 0.45 0.210 0.054],'fontunits','normalized', 'fontSize', 0.36, 'visible', 'off', 'UserData', struct('remove','left_samp_rem'));
ST_L2_SEL = uicontrol(stats_GUI,'Style','pushbutton','String', 'Select','Units', 'normalized','Position',[0.52 0.45 0.210 0.054],'fontunits','normalized', 'fontSize', 0.36, 'visible', 'off','UserData', struct('select','right_samp_sel'));
ST_L2_REM = uicontrol(stats_GUI,'Style','pushbutton','String', 'Remove','Units', 'normalized','Position',[0.75 0.45 0.210 0.054],'fontunits','normalized', 'fontSize', 0.36, 'visible', 'off', 'UserData', struct('remove','right_samp_rem'));

% Boxes & layout for alpha & threshold values
ST_CONT = uipanel(stats_GUI,'Units', 'normalized','Position',[0.046 0.37 0.44 0.07],'HighLightColor',[0.78 0.78 0.78],'BorderType', 'line','BackgroundColor','w');
ST_CONT_txt  = uicontrol(stats_GUI,'Style','text','String', 'Contrast:','Units', 'normalized', 'Position',[0.095 0.38 0.38 0.04],'fontunits','normalized', 'fontSize', 0.55, 'HorizontalAlignment','Left','backgroundcolor','w');
ST_CONT_val  = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.278 0.382 0.18 0.045],'fontunits','normalized', 'fontSize', 0.50);
ST_ALP = uipanel(stats_GUI,'Units', 'normalized','Position',[0.52 0.37 0.44 0.07],'HighLightColor',[0.78 0.78 0.78],'BorderType', 'line','BackgroundColor','w');
ST_ALP_txt  = uicontrol(stats_GUI,'Style','text','String', 'Alpha:','Units', 'normalized', 'Position',[0.583 0.38 0.35 0.04],'fontunits','normalized', 'fontSize', 0.55, 'HorizontalAlignment','Left','backgroundcolor','w');
ST_ALP_val  = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.755 0.382 0.18 0.045],'fontunits','normalized', 'fontSize', 0.50);

% Type of threshold selection pop-up menu and conditional value
ST_THRES_TXT = uicontrol(stats_GUI,'Style','text','String', 'Threshold type: ','Units', 'normalized', 'Position',[0.098 0.298 0.38 0.04],'fontunits','normalized', 'fontSize', 0.58, 'HorizontalAlignment','Left','backgroundcolor','w');
ST_THRES_POP = uicontrol(stats_GUI,'Style','popupmenu','String', {'Uncorrected (Parametric)', 'FDR (Parametric)', 'Bonferroni (Parametric)', 'Uncorrected (Non-parametric)','FDR (Non-parametric)','NBS FWE (Non-parametric)','NBS TFCE (Non-parametric)'},'Units', 'normalized', 'Position',[0.358 0.295 0.6 0.05],'fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w');
ST_THRES_VAL_TXT = uicontrol(stats_GUI,'Style','text','String', 'Primary threshold (p-value): ','Units', 'normalized', 'Position',[0.098 0.23 0.5 0.04],'fontunits','normalized', 'fontSize', 0.58, 'HorizontalAlignment','Left','backgroundcolor','w', 'enable', 'off');
ST_THRES_VAL_UNI = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.76 0.234 0.2 0.04],'fontunits','normalized', 'fontSize', 0.58,'backgroundcolor','w', 'enable', 'off');
ST_PERM_TXT = uicontrol(stats_GUI,'Style','text','String', 'Number of permutations: ','Units', 'normalized', 'Position',[0.098 0.165 0.38 0.04],'fontunits','normalized', 'fontSize', 0.58, 'HorizontalAlignment','Left','backgroundcolor','w', 'enable', 'off');
ST_PERM_VAL = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.76 0.169 0.2 0.04],'fontunits','normalized', 'fontSize', 0.58,'backgroundcolor','w','enable', 'off');

% Run button
RES_RUN = uicontrol(stats_GUI, 'Style', 'pushbutton', 'String', 'Run','Units', 'normalized','Position',[0.4 0.05 0.210 0.054],'fontunits','normalized', 'fontSize', 0.36);

% Callback actions
set(ST_test_type, 'callback', @test_type);
set(ST_THRES_POP, 'callback', @threshold_type);
set(ST_lst_0, 'callback', @live_select_0)
set(ST_lst_1, 'callback', @live_select_1)
set(ST_lst_2, 'callback', @live_select_2)
set(ST_L0_SEL, 'callback', @(src, event) select_caller(get(src, 'UserData')));
set(ST_L1_SEL, 'callback', @(src, event) select_caller(get(src, 'UserData')));
set(ST_L2_SEL, 'callback', @(src, event) select_caller(get(src, 'UserData')));
set(ST_L0_REM, 'callback', @(src, event) remove_caller(get(src, 'UserData')));
set(ST_L1_REM, 'callback', @(src, event) remove_caller(get(src, 'UserData')));
set(ST_L2_REM, 'callback', @(src, event) remove_caller(get(src, 'UserData')));
set(RES_RUN, 'callback', @run);
% warning('off','backtrace')

M0 = {}; % variable to store the matrices for One-sample t-test
M1 = {}; % variable to store the matrices set 1 Paired & Two-sample t-test
M2 = {}; % variable to store the matrices set 2 Paired & Two-sample t-test

% Variables to store present selection of matrices from listboxes
selection_0 = '';
selection_1 = '';
selection_2 = '';
matrices_0 = [];
matrices_1 = [];
matrices_2 = [];

% Function to select respective call button
function select_caller(data)
    switch (data.select)
        case 'one_samp_sel'
            file_selector(M0, matrices_0, ST_lst_0, ST_L0_CTR,'one_samp_sel');
            
        case 'left_samp_sel'
            file_selector(M1, matrices_1, ST_lst_1, ST_L1_CTR,'left_samp_sel');
            
        case 'right_samp_sel'
            file_selector(M2, matrices_2, ST_lst_2, ST_L2_CTR,'right_samp_sel');
    end
end

% Function to select respective remove button
function remove_caller(data)
    switch (data.remove)
        case 'one_samp_rem'
            file_remove(selection_0, M0, ST_lst_0, ST_L0_CTR,'one_samp_sel');
            
        case 'left_samp_rem'
            file_remove(selection_1, M1, ST_lst_1, ST_L1_CTR,'left_samp_sel');
            
        case 'right_samp_rem'
            file_remove(selection_2, M2, ST_lst_2, ST_L2_CTR,'right_samp_sel');
    end
end

% Function to select and add files to listboxes
function file_selector(M_VAR, matrix, disp_box, disp_str, case_maker)

    % First case: Initial selection of data--------------------------------
    if isempty(M_VAR)
        
        % Checking if there exist pre-selected *.mat files
        M_VAR = select_mat_file();        % Select *.mat files
        M_VAR = unique(M_VAR);            % Remove duplicates

        % If *.mat files have been selected, perform multiple variable and dimension checks
        if ~isempty(M_VAR)          

            % Check if *.mat file consists of multiple variables
            mv_flag = multi_var_check(M_VAR);
            if mv_flag == 0   

                % Continue if the selected files do not contain multiple variables
                for i = 1:size(M_VAR,1)
                    % M(i).m = struct2array(load(M_VAR{i,:}));
                    tmp = load(M_VAR{i,:}); fn  = fieldnames(tmp); M(i).m = tmp.(fn{1});
                end
    
                try
                    matrix = cat(3,M(:).m);
                    if size(matrix,1) ~= size(matrix,2)
                        warning('Matrices are not square.');
                        clear M matrix   
                        M_VAR = {};
                    end
                catch
                    warning('Matrices have different dimensions.');
                    clear M  
                    M_VAR = {};
                end
                
            elseif mv_flag == 1
                M_VAR = {};
                warning('Each selected *.mat file must contain exactly one variable.');
            end 
        end
               
        % Updating the GUI 
        if ~exist('M_VAR', 'var') || isempty(M_VAR)
            % If all files selection was rejected during checks, reset GUI
            disp('No *.mat file(s) selected');
            set(disp_str, 'String', '0 ROIs x 0 subjects');
            set(disp_str, 'ForegroundColor',[0.773, 0.353, 0.067]);     
            M_VAR = {};
        elseif isempty(M_VAR{1}) 
            % If all files selection was rejected during checks, reset GUI
            disp('No *.mat file(s) selected');
            set(disp_str, 'String', '0 ROIs x 0 subjects');
            set(disp_str, 'ForegroundColor',[0.773, 0.353, 0.067]);     
            M_VAR = {};
        else
            % Show the number of *.mat files selected & update GUI
            fprintf('Number of .mat files selected: %d\n', size(M_VAR,1));
            set(disp_box,'String', M_VAR);
            set(disp_box,'Value', []);

            % Update the ROI x ROI x Subjects number in GUI
            set(disp_str, 'String', strcat(num2str(size(matrix,2)), ' ROIs x',32, num2str(size(matrix,3)),' subjects'));
            set(disp_str, 'ForegroundColor',[0.219, 0.341, 0.137]); 
        end
        
    % Second case: Add new matrices----------------------------------------
    else
        % Select new files via function
        new_M_VAR = select_mat_file();
        
        % If new files are selected then proceed 
        if ~isempty(new_M_VAR) && numel(new_M_VAR) >= 1 && ~isempty(new_M_VAR{1})
               
            % Check for multiple variables within selected files   
            if multi_var_check(new_M_VAR) ~= 1        
                
                % Continue if the selected files do not contain multiple variables
                for i = 1:size(new_M_VAR,1)
                    % M(i).m = struct2array(load(new_M_VAR{i,:}));
                    tmp = load(new_M_VAR{i,:}); fn  = fieldnames(tmp); M(i).m = tmp.(fn{1});
                end
    
                try
                    new_matrices = cat(3,M(:).m);
                    if size(new_matrices,1) ~= size(new_matrices,2)
                        warning('Matrices are not square.');
                        clear M 
                        new_M_VAR = {};
                    end
                catch
                    warning('Matrices have different dimensions.');
                    clear M  
                    new_M_VAR = {};
                end
                
                % Concatenate old and new matrices
                try
                    matrix = cat(3,matrix,new_matrices);
                    M_VAR = vertcat(M_VAR, new_M_VAR);
                    % Updating the GUI 
                    fprintf('Number of .mat files selected: %d\n', size(new_M_VAR,1));
                    set(disp_box,'String', M_VAR);
                    set(disp_box,'Value', []);
        
                    % Update the ROI x ROI x Subjects number
                    set(disp_str, 'String', strcat(num2str(size(matrix,2)), ' ROIs x',32, num2str(size(matrix,3)),' subjects'));
                    set(disp_str, 'ForegroundColor',[0.219, 0.341, 0.137]); 
                    clear M new_M_VAR
                catch
                    warning('Matrices have different numbers of ROIs.');
                    clear new_matrices new_M_VAR M
                end
            else
                warning('Selected *.mat file(s) consist(s) of multiple variables, please select *.mat files each containing only one variable.');
            end
        else          
            disp('No files added.');          
        end
    end
    
    switch (case_maker)    
        case 'one_samp_sel'
            M0 = M_VAR;
            matrices_0 = matrix;
        case 'left_samp_sel'
            M1 = M_VAR;
            matrices_1 = matrix;
        case 'right_samp_sel'
            M2 = M_VAR;
            matrices_2 = matrix;      
    end
end

% Function to perform removal of files from listboxes
function file_remove(sel_var, M_VAR, disp_box,disp_str,case_maker)
    if isempty(sel_var) && isempty(M_VAR)
        warning('There are no files present to remove, please select *.mat files to perform statistical analysis.');
    elseif isempty(sel_var) && ~isempty(M_VAR)
        warning('There are no selected matrices to remove from the list, please select matrices once again.');
    else
        M_VAR(sel_var,:) = [];
        holder = size(sel_var);
        fprintf('Number of .mat files removed: %d\n', holder(2));
    
        % Rebuild matrix from remaining files
        matrix = [];
        if ~isempty(M_VAR)
           try
               for i = 1:size(M_VAR,1)
                   tmp = load(M_VAR{i,:}); fn = fieldnames(tmp); M(i).m = tmp.(fn{1});
               end
               matrix = cat(3, M(:).m);
               if size(matrix,1) ~= size(matrix,2)
                   warning('Matrices are not square.');
                   matrix = [];
               end
               clear M
           catch
               warning('Matrices have different dimensions.');
               matrix = [];
               clear M
           end
        end
        % --------------------------------------------------
              
        set(disp_box,'Value', []);
        set(disp_box,'String', M_VAR);
        sel_var = {};
        
        if ~isempty(M_VAR)
            S  = whos('-file', M_VAR{1,:});  dims = S(1).size;
            subs = 0;
            for i = 1:length(M_VAR)
                temp = whos('-file', M_VAR{i,:}); temp_dim = temp(1).size;
                if numel(temp_dim) == 2, subs = subs + 1;
                elseif numel(temp_dim) == 3, subs = subs + temp_dim(3);
                else, warning('Files must be ROI×ROI or ROI×ROI×Subjects.');
                end
            end
            set(disp_str, 'String', sprintf('%d ROIs x %d subjects', dims(1), subs));
            set(disp_str, 'ForegroundColor',[0.219, 0.341, 0.137]);
        end
    end

    if isempty(M_VAR)
        set(disp_str, 'String', '0 ROIs x 0 subjects');
        set(disp_str, 'ForegroundColor',[0.773, 0.353, 0.067]);
    end
   
    selection_0 = '';  selection_1 = '';  selection_2 = '';
  
    switch (case_maker)
        case 'one_samp_sel'
            M0 = M_VAR; if isempty(M0), M0 = {}; end
            matrices_0 = matrix;
        case 'left_samp_sel'
            M1 = M_VAR; if isempty(M1), M1 = {}; end
            matrices_1 = matrix;
        case 'right_samp_sel'
            M2 = M_VAR; if isempty(M2), M2 = {}; end
            matrices_2 = matrix;
    end
end


% Variables to store live selections from lists
function live_select_0(~,~)
    index = get(ST_lst_0, 'Value');% Retrieves the users selection LIVE
    selection_0 = index;                % Variable for full selection
end
function live_select_1(~,~)
    index = get(ST_lst_1, 'Value');% Retrieves the users selection LIVE
    selection_1 = index;                % Variable for full selection
end
function live_select_2(~,~)
    index = get(ST_lst_2, 'Value');% Retrieves the users selection LIVE
    selection_2 = index;                % Variable for full selection
end


% Function to configure GUI based on selected test-type 
function test_type(~,~)
    
    % Extract the current Test mode selected by user
    contender = (ST_test_type.String{ST_test_type.Value});

    % Action relative to test type
    if strcmpi(contender, 'Paired t-test')
        
        % If Paired T Test is selected
        disp('Paired t-test selected.');
        
        % Reset GUI 
        set([ST_lst_0,ST_L0_CTR,ST_L0_SEL,ST_L0_REM],'visible', 'off');        
        set([ST_lst_1,ST_lst_2,ST_L1_CTR,ST_L2_CTR,...
            ST_L1_SEL,ST_L1_REM,ST_L2_SEL,ST_L2_REM],'visible', 'on','enable', 'on'); 
        set([ST_THRES_POP,ST_THRES_TXT,ST_CONT_txt,...
            ST_CONT_val,ST_ALP_txt,ST_ALP_val,RES_RUN],'enable', 'on');
        set([ST_L0_CTR, ST_L1_CTR,ST_L2_CTR], 'String', '0 ROIs x 0 subjects',...
            'ForegroundColor',[0.773, 0.353, 0.067]);
        set([ST_CONT_val, ST_ALP_val,ST_PERM_VAL, ST_THRES_VAL_UNI]...
            , 'String', []);
        
        % Reset Variables 
        M0 = {};
        M1 = {};
        M2 = {};
        selection_0 = '';
        selection_1 = '';
        selection_2 = '';
        
        % Reset GUI elements
        set(ST_lst_0,'String', M0,'Value', []);
        set(ST_lst_1,'String', M1,'Value', []);
        set(ST_lst_2,'String', M2,'Value', []);
       
    elseif strcmpi(contender, 'One-sample t-test')
        
        % If One-sample T Test is selected
        disp('One-sample t-test selected.');
        
        % Reset GUI
        set([ST_lst_0,ST_L0_CTR,ST_L0_SEL,ST_L0_REM],'visible', 'on');        
        set([ST_lst_1,ST_lst_2,ST_L1_CTR,ST_L2_CTR,...
            ST_L1_SEL,ST_L1_REM,ST_L2_SEL,ST_L2_REM],'visible', 'off');
        set([ST_L0_CTR, ST_L1_CTR,ST_L2_CTR], 'String', '0 ROIs x 0 subjects',...
            'ForegroundColor',[0.773, 0.353, 0.067]);
        set([ST_CONT_val, ST_ALP_val,ST_PERM_VAL, ST_THRES_VAL_UNI],...
            'String', []);
        set([ST_THRES_POP,ST_THRES_TXT,ST_CONT_txt,ST_CONT_val,...
            ST_ALP_txt,ST_ALP_val,RES_RUN],'enable', 'on');
        
        % Reset variables
        M0 = {};
        M1 = {};
        M2 = {};
        selection_0 = '';
        selection_1 = '';
        selection_2 = '';
        
        % Reset GUI elements
        set(ST_lst_0,'String', M0,'Value', []);
        set(ST_lst_1,'String', M1,'Value', []);
        set(ST_lst_2,'String', M2,'Value', []);
                
    elseif strcmpi(contender, 'Two-sample t-test')
        
        % If Two-sample T Test is selected
        disp('Two-sample t-test selected.');
        
        % Reset GUI 
        set([ST_lst_0,ST_L0_CTR,ST_L0_SEL,ST_L0_REM],'visible', 'off');    
        set([ST_lst_1,ST_lst_2,ST_L1_CTR,ST_L2_CTR,...
            ST_L1_SEL,ST_L1_REM,ST_L2_SEL,ST_L2_REM],'visible', 'on','enable', 'on');        
        set([ST_L0_CTR, ST_L1_CTR,ST_L2_CTR], 'String', '0 ROIs x 0 subjects',...
            'ForegroundColor',[0.773, 0.353, 0.067]);       
        set([ST_CONT_val, ST_ALP_val,ST_PERM_VAL, ST_THRES_VAL_UNI], 'String', []);        
        set([ST_THRES_POP,ST_THRES_TXT,ST_CONT_txt,ST_CONT_val,ST_ALP_txt,ST_ALP_val,RES_RUN],'enable', 'on');

        % Reset Variables
        M0 = {};
        M1 = {};
        M2 = {};
        selection_0 = '';
        selection_1 = '';
        selection_2 = '';
        
        % Reset GUI elements
        set(ST_lst_0,'String', M0,'Value', []);
        set(ST_lst_1,'String', M1,'Value', []);
        set(ST_lst_2,'String', M2,'Value', []);
    end    
end

% Type of threshold
function threshold_type(~,~)
    
    approach = (ST_THRES_POP.String{ST_THRES_POP.Value});
    
    if strcmpi(approach, 'Uncorrected (Parametric)') || strcmpi(approach, 'FDR (Parametric)') || strcmpi(approach, 'Bonferroni (Parametric)') || strcmpi(approach, 'NBS FWE (Non-parametric)') || strcmpi(approach, 'NBS TFCE (Non-parametric)') 
        set(ST_PERM_TXT, 'enable', 'off');
        set(ST_PERM_VAL, 'enable', 'off', 'String', []);
        
    elseif strcmpi(approach, 'Uncorrected (Non-parametric)') || strcmpi(approach, 'FDR (Non-parametric)') 
        set(ST_PERM_TXT, 'enable', 'on');
        set(ST_PERM_VAL, 'enable', 'on', 'String', []);

    end
    
    if strcmpi(approach, 'Uncorrected (Parametric)') || strcmpi(approach, 'FDR (Parametric)') || strcmpi(approach, 'Bonferroni (Parametric)') || strcmpi(approach, 'Uncorrected (Non-parametric)') || strcmpi(approach, 'FDR (Non-parametric)')
        set(ST_THRES_VAL_TXT, 'enable', 'off');
        set(ST_THRES_VAL_UNI, 'enable', 'off', 'String', []);

    elseif strcmpi(approach, 'NBS FWE (Non-parametric)') 
        set([ST_THRES_VAL_TXT,ST_PERM_TXT], 'enable', 'on');
        set([ST_THRES_VAL_UNI,ST_PERM_VAL], 'enable', 'on', 'String', []);

    elseif strcmpi(approach, 'NBS TFCE (Non-parametric)') 
        set([ST_THRES_VAL_TXT,ST_PERM_TXT], 'enable', 'off');
        set([ST_THRES_VAL_UNI,ST_PERM_VAL], 'enable', 'off', 'String', []);
    end
    
    if strcmpi(approach, 'Uncorrected (Non-parametric)') || strcmpi(approach, 'FDR (Non-parametric)') || strcmpi(approach, 'NBS FWE (Non-parametric)') || strcmpi(approach, 'NBS TFCE (Non-parametric)')
        set(RES_RUN,'enable', 'off');
        warning('Non-parametric analysis is under development. Please wait for a future update.');
    else
        set(RES_RUN,'enable', 'on');
    end
    
end

% Run Test
function run(~,~)
    test_type = (ST_test_type.String{ST_test_type.Value}); % Type of test - paired, one , two 
    
    %----------------------------------------------------------------------
    % PAIRED T-TEST
    %----------------------------------------------------------------------
    if strcmpi(test_type, 'Paired t-test')          
        if ~isempty(M1) && ~isempty(M2)
            flag_parameters = parameter_check(); % alpha & contrast
            if flag_parameters == 1
                flag_thres_parameters = perm_thres_check(); % permutations / primary threshold
                if flag_thres_parameters == 1
                    % Validate ROI dimensions
                    if isempty(matrices_1) || isempty(matrices_2)
                        warning('Matrices are empty.'); 
                        return;
                    end
                    if size(matrices_1,1) ~= size(matrices_2,1) || size(matrices_1,2) ~= size(matrices_2,2)
                        warning('The ROI×ROI dimensions differ between the two lists. Please select matrices with the same number of ROIs.');
                        set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]);
                        return;
                    end
            
                    % Validate subject counts for pairing
                    n1 = size(matrices_1,3);
                    n2 = size(matrices_2,3);
                    if n1 ~= n2
                        warning('Paired t-test requires the same number of subjects. Left: %d, Right: %d.', n1, n2);
                        set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]);
                        return;
                    end
            
                    % Map UI label -> threshold key 
                    thr_key = threshold_key(ST_THRES_POP.String{ST_THRES_POP.Value});
                    if isempty(thr_key)
                        warning('Selected threshold type is not supported for this analysis.');
                        return;
                    end
            
                    % Run
                    set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]); 
                    groups = {matrices_1, matrices_2};
                    alpha = safe_str2num_scalar(ST_ALP_val.String);
                    [thresholded,pval,tval,conval] = tmfc_ttest(groups, ...
                        str2num(ST_CONT_val.String), alpha, thr_key);
                    fprintf('Generating results...\n');
                    tmfc_results_GUI(thresholded,pval,tval,conval,alpha,thr_key);
                    clear thresholded pval tval conval groups alpha thr_key
                end
            end
        elseif ~isempty(M1) && isempty(M2)
            warning('Please select the second set of matrix files to run a paired t-test.');
        elseif isempty(M1) && ~isempty(M2)
            warning('Please select the first set of matrix files to run a paired t-test.');
        else
            warning('Please select matrix files to run a paired t-test.');
        end

    %----------------------------------------------------------------------
    % ONE-SAMPLE T-TEST
    %----------------------------------------------------------------------
    elseif strcmpi(test_type, 'One-sample t-test')
        
        % Check for matrix
        if ~isempty(M0)
            flag_parameters = parameter_check(); % Alpha, contrast check

            if flag_parameters == 1
                flag_thres_parameters = perm_thres_check(); % permutations, threshold check

                if flag_thres_parameters == 1
                    set(ST_L0_CTR, 'ForegroundColor',[0.219, 0.341, 0.137]); 
                    thr_key = threshold_key(ST_THRES_POP.String{ST_THRES_POP.Value});
                    if isempty(thr_key)
                        warning('Selected threshold type is not supported for this analysis.');
                        return;
                    else
                        alpha = safe_str2num_scalar(ST_ALP_val.String);
                        [thresholded,pval,tval,conval] = tmfc_ttest(matrices_0,str2num(ST_CONT_val.String),alpha,thr_key);
                        fprintf('Generating results...\n');
                        tmfc_results_GUI(thresholded,pval,tval,conval,alpha,thr_key);
                    end
                    clear thresholded pval tval conval alpha thr_key
                end
            end
            
        else
            warning('Please select matrix files to run a one-sample t-test.');
        end

    %----------------------------------------------------------------------
    % TWO-SAMPLE T-TEST
    %----------------------------------------------------------------------
    elseif strcmpi(test_type, 'Two-sample t-test')

        % Check for matrix
        if ~isempty(M1) && ~isempty(M2)
            flag_parameters = parameter_check();
            if flag_parameters == 1
                flag_thres_parameters = perm_thres_check();
                if flag_thres_parameters == 1
            
                    if isempty(matrices_1) || isempty(matrices_2)
                        warning('Matrices are empty.');
                        return;
                    end
                    if size(matrices_1,1) ~= size(matrices_2,1) || size(matrices_1,2) ~= size(matrices_2,2)
                        warning('The ROI×ROI dimensions differ between the two lists. Please select matrices with the same number of ROIs.');
                        set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]);
                        return;
                    end
            
                    set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]);
                    groups = {matrices_1, matrices_2};
            
                    thr_key = threshold_key(ST_THRES_POP.String{ST_THRES_POP.Value});
                    if isempty(thr_key)
                        warning('Selected threshold type is not supported for this analysis.');
                        return;
                    end
                    alpha = safe_str2num_scalar(ST_ALP_val.String);
                    [thresholded,pval,tval,conval] = tmfc_ttest2(groups, ...
                        str2num(ST_CONT_val.String),alpha,thr_key);
                    fprintf('Generating results...\n');
                    tmfc_results_GUI(thresholded,pval,tval,conval,alpha,thr_key);
                    clear thresholded pval tval conval groups alpha thr_key
                end
            end          
        elseif ~isempty(M1) && isempty(M2)
            warning('Please select the second set of matrix files to run a two-sample t-test.');
        elseif isempty(M1) && ~isempty(M2)
            warning('Please select the first set of matrix files to run a two-sample t-test.'); 
        else
            warning('Please select matrix files to run a two-sample t-test.');
        end 
    else
        warning('Files must be ROI×ROI or ROI×ROI×Subjects.');
    end 
end

% Check contrast and alpha
function flag = parameter_check(~,~)
    flag = 0;
    test_type = (ST_test_type.String{ST_test_type.Value});
    val_contrast = str2num(ST_CONT_val.String);   
    val_alpha = safe_str2num_scalar(ST_ALP_val.String);

    if isempty(val_contrast)
        warning('Please enter numeric values for the contrast.');
        return;
    end
    if isempty(val_alpha) || ~isscalar(val_alpha) || ~isfinite(val_alpha) || val_alpha < 0 || val_alpha > 1
        warning('Please enter an alpha value between 0 and 1.');
        return;
    end
    if any(strcmpi(test_type, {'Paired t-test','Two-sample t-test'}))
        if numel(val_contrast) ~= 2
            warning('Please enter two contrast values.');
            return;
        end
    elseif strcmpi(test_type, 'One-sample t-test')
        if numel(val_contrast) > 1
            warning('One-sample t-test accepts only a single contrast value.');
            return;
        end
    end
    flag = 1;
end 

% Handle alpha input
function x = safe_str2num_scalar(s)
    s = strtrim(s);
    if isempty(s), x = NaN; return; end
    % allow digits, decimal, + - * / ^, parentheses, whitespace, and e/E for exponents
    if isempty(regexp(s,'^[\d\.\+\-\*\/\^\(\)\seE]+$','once'))
        x = NaN; return;
    end
    tmp = str2num(s); 
    if isempty(tmp) || ~isscalar(tmp) || ~isfinite(tmp)
        x = NaN; return;
    end
    x = tmp;
end

% Validate permutation / primary threshold inputs 
function flag = perm_thres_check(~,~)
    flag = 0;
    thresh_type = (ST_THRES_POP.String{ST_THRES_POP.Value});

    switch thresh_type
        case {'Uncorrected (Non-parametric)', 'FDR (Non-parametric)'}
            permutation_value = str2double(ST_PERM_VAL.String);
            is_posint = ~isnan(permutation_value) && permutation_value > 0 && mod(permutation_value,1)==0;
            if is_posint
                warning('Non-parametric analysis is under development. Please wait for a future update.');
                flag = 0;
            else
                warning('Please enter a positive integer for the number of permutations.');
            end

        case 'NBS FWE (Non-parametric)'
            permutation_value = str2double(ST_PERM_VAL.String);
            is_posint = ~isnan(permutation_value) && permutation_value > 0 && mod(permutation_value,1)==0;

            threshold_value = safe_str2num_scalar(ST_THRES_VAL_UNI.String);
            is_valid_thresh = ~isnan(threshold_value) && threshold_value > 0 && threshold_value <= 1.0;

            if ~is_posint
                warning('Please enter a positive integer for the number of permutations.');
            elseif ~is_valid_thresh
                warning('Please enter a primary threshold (p-value) between 0 and 1.');
            else
                warning('Non-parametric analysis is under development. Please wait for a future update.');
                flag = 0;
            end

        case 'NBS TFCE (Non-parametric)'
            warning('Non-parametric analysis is under development. Please wait for a future update.');
            flag = 0;

        case {'Uncorrected (Parametric)', 'FDR (Parametric)', 'Bonferroni (Parametric)'}
            flag = 1;

        otherwise
            warning('Unsupported threshold type.');
    end
end

uiwait(stats_GUI);

function on_close(h, ~)
    try, uiresume(h); end
    delete(h);
end

end

%% Select *.mat files from the user via spm_select
function list_sel = select_mat_file(~,~)  
    files = spm_select(inf,'.mat','Select matrix files (.mat) for computation',{},pwd,'.');
    list_sel = {};
    list_sel = cellstr(files);
end

%% Check if the selected *.mat files contain multiple variables
%  0 = single variable in all files
%  1 = at least one file has 0 or >1 variables
% -1 = empty input (no files)
function flag = multi_var_check(input_file)
    if isempty(input_file) || isempty(input_file{1})
        flag = -1;  % nothing to check
        return;
    end
    flag = 0;
    for i = 1:size(input_file,1)
        varlist = who('-file', input_file{i,:});
        if numel(varlist) ~= 1
            if isempty(varlist)
                fprintf('No variables found in file: %s\n', input_file{i,:});
            else
                fprintf('Multiple variables found in file: %s\n', input_file{i,:});
            end
            flag = 1;
            break;
        end
    end
end

%% Convert internal labeling format to tmfc_ttest() labeling format
function small_string = threshold_key(input_string)
    small_string = '';    
    switch input_string 
        case 'Uncorrected (Parametric)'
            small_string = 'uncorr';
        
        case 'FDR (Parametric)'
            small_string = 'FDR';
        
        case 'Bonferroni (Parametric)'
            small_string = 'Bonf';            
    end
end