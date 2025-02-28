function tmfc_statistics_GUI()

% GUI Initialization
stats_GUI = figure('Name', 'TMFC: Results', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.45 0.25 0.22 0.56],'MenuBar', 'none','ToolBar', 'none','color','w'); movegui(stats_GUI, 'center');
RES_T1  = uicontrol(stats_GUI,'Style','text','String', 'TMFC results','Units', 'normalized', 'Position',[0.270 0.93 0.460 0.05],'fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w', 'FontWeight', 'bold');

% Pop up menu to select type of Test
ST_test_type  = uicontrol(stats_GUI,'Style','popupmenu','String', {'One-sample t-test', 'Paired t-test', 'Two-sample t-test'},'Units', 'normalized', 'Position',[0.045 0.87 0.91 0.05],'fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w');
 
% List boxes to show (.mat) file selection
ST_lst_0 = uicontrol(stats_GUI , 'Style', 'listbox', 'String', '','Max', 100,'value',[],'Units', 'normalized', 'Position',[0.045 0.56 0.91 0.300],'fontunits','points', 'fontSize',11);
ST_lst_1 = uicontrol(stats_GUI , 'Style', 'listbox', 'String', '','Max', 100,'value',[],'Units', 'normalized', 'Position',[0.045 0.56 0.440 0.300],'fontunits','points', 'fontSize', 11,'visible','off');
ST_lst_2 = uicontrol(stats_GUI , 'Style', 'listbox', 'String', '','Max', 100,'value',[],'Units', 'normalized', 'Position',[0.52 0.56 0.440 0.300],'fontunits','points', 'fontSize',11,'visible','off');

% Counter of subjects selected
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

% Boxes & Layout for Alpha & threshold values
ST_CONT = uipanel(stats_GUI,'Units', 'normalized','Position',[0.046 0.37 0.44 0.07],'HighLightColor',[0.78 0.78 0.78],'BorderType', 'line','BackgroundColor','w');
ST_CONT_txt  = uicontrol(stats_GUI,'Style','text','String', 'Contrast: ','Units', 'normalized', 'Position',[0.095 0.38 0.38 0.04],'fontunits','normalized', 'fontSize', 0.55, 'HorizontalAlignment','Left','backgroundcolor','w');
ST_CONT_val  = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.278 0.382 0.18 0.045],'fontunits','normalized', 'fontSize', 0.50);
ST_ALP = uipanel(stats_GUI,'Units', 'normalized','Position',[0.52 0.37 0.44 0.07],'HighLightColor',[0.78 0.78 0.78],'BorderType', 'line','BackgroundColor','w');
ST_ALP_txt  = uicontrol(stats_GUI,'Style','text','String', 'Alpha: ','Units', 'normalized', 'Position',[0.583 0.38 0.35 0.04],'fontunits','normalized', 'fontSize', 0.55, 'HorizontalAlignment','Left','backgroundcolor','w');
ST_ALP_val  = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.755 0.382 0.18 0.045],'fontunits','normalized', 'fontSize', 0.50);

% Type of Threshold selection Pop Up menu and conditional value
ST_THRES_TXT = uicontrol(stats_GUI,'Style','text','String', 'Threshold type: ','Units', 'normalized', 'Position',[0.098 0.298 0.38 0.04],'fontunits','normalized', 'fontSize', 0.58, 'HorizontalAlignment','Left','backgroundcolor','w');
ST_THRES_POP = uicontrol(stats_GUI,'Style','popupmenu','String', {'Uncorrected (Parametric)', 'FDR (Parametric)', 'Bonferroni (Parametric)', 'Uncorrected (Non-Parametric)','FDR (Non-Parametric)','NBS FWE (Non-Parametric)','NBS TFCE (Non-Parametric)'},'Units', 'normalized', 'Position',[0.358 0.295 0.6 0.05],'fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w');
ST_THRES_VAL_TXT = uicontrol(stats_GUI,'Style','text','String', 'Primary Threshold Value (Pval): ','Units', 'normalized', 'Position',[0.098 0.23 0.5 0.04],'fontunits','normalized', 'fontSize', 0.58, 'HorizontalAlignment','Left','backgroundcolor','w', 'enable', 'off');
ST_THRES_VAL_UNI = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.76 0.234 0.2 0.04],'fontunits','normalized', 'fontSize', 0.58,'backgroundcolor','w', 'enable', 'off');
ST_PERM_TXT = uicontrol(stats_GUI,'Style','text','String', 'Permutations: ','Units', 'normalized', 'Position',[0.098 0.165 0.38 0.04],'fontunits','normalized', 'fontSize', 0.58, 'HorizontalAlignment','Left','backgroundcolor','w', 'enable', 'off');
ST_PERM_VAL = uicontrol(stats_GUI,'Style','edit','String', '','Units', 'normalized', 'Position',[0.76 0.169 0.2 0.04],'fontunits','normalized', 'fontSize', 0.58,'backgroundcolor','w','enable', 'off');

% Run Results
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
warning('off','backtrace')

M0 = {}; % variable to store the matrices for One-sample t-test
M1 = {}; % variable to store the matrices set 1 Paired & Two-sample t-test
M2 = {}; % variable to store the matrices set 2 Paired & Two-sample t-test

% Variables to store present selection of matrices from list
selection_0 = '';
selection_1 = '';
selection_2 = '';
matrices_0 = 0;
matrices_1 = 0;
matrices_2 = 0;

% Function to select respective call button
function select_caller(data)
% Format: file_selector(M_VAR, matrix, disp_box, disp_str)
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

% Function to perform selection and addition of files to list boxes
function file_selector(M_VAR, matrix, disp_box, disp_str, case_maker)
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   % First case: Initial selection of data
    if isempty(M_VAR)
        
        % Checking if there exist pre-selected (.mat) files
        M_VAR = select_mat_file();         % Select (.mat) files
        M_VAR = unique(M_VAR);      % Remove duplicates

        
        % If (.mat) files have been selected, perform multiple variable and dimension checks
        if ~isempty(M_VAR)          

            % Check if (.mat) file consists of multiple variables
            if multi_var_check(M_VAR) == 0   

                % Continue if the selected files do not contain multiple variables
                for i = 1:size(M_VAR,1)
                    M(i).m = struct2array(load(M_VAR{i,:}));
                end
    
                try
                    matrix = cat(3,M(:).m);
                    if size(matrix,1) ~= size(matrix,2)
                        warning('Matrices are not square')
                        clear M matrices   
                        M_VAR = {};
                    end
                catch
                    warning('Matrices have different dimensions')
                    clear M  
                    M_VAR = {};
                end
                
            elseif multi_var_check(M_VAR) == 1
                % Warning if file has MULTIPLE VARIABLES within 
                M_VAR = {};
                warning('Selected *.mat file(s) consist(s) of multiple variables, please select *.mat files each containing only one variable.');
            end 
            
        end
               
        % Updating the GUI 
        if ~exist('M_VAR', 'var') || isempty(M_VAR)
            % If all files selection was rejected during checks, reset GUI
            disp('No (.mat) file(s) selected');
            set(disp_str, 'String', '0 ROIs x 0 subjects');
            set(disp_str, 'ForegroundColor',[0.773, 0.353, 0.067]);     
            M_VAR = {};
        elseif isempty(M_VAR{1}) 
            % If all files selection was rejected during checks, reset GUI
            disp('No (.mat) file(s) selected');
            set(disp_str, 'String', '0 ROIs x 0 subjects');
            set(disp_str, 'ForegroundColor',[0.773, 0.353, 0.067]);     
            M_VAR = {};
        else
            % Show the number of (.mat) files selected & update GUI
            fprintf('Number of (.mat) files selected are: %d \n', size(M_VAR,1));
            set(disp_box,'String', M_VAR);
            set(disp_box,'Value', []);

            % Update the ROI x ROI x Subjects number in GUI
            set(disp_str, 'String', strcat(num2str(size(matrix,2)), ' ROIs x',32, num2str(size(matrix,3)),' subjects'));
            set(disp_str, 'ForegroundColor',[0.219, 0.341, 0.137]); 
        end
        
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    % Second case: Add new matrices
    else
        % Select new files via function
        new_M_VAR = select_mat_file();
        
        % If new files are selected then proceed 
        if ~isempty(new_M_VAR{1})
               
            % Check for multiple variables within selected files   
            if multi_var_check(new_M_VAR) ~= 1        
                
                % Continue if the selected files do not contain multiple variables
                for i = 1:size(new_M_VAR,1)
                    M(i).m = struct2array(load(new_M_VAR{i,:}));
                end
    
                try
                    new_matrices = cat(3,M(:).m);
                    if size(new_matrices,1) ~= size(new_matrices,2)
                        warning('Matrices are not square.')
                        clear M 
                        new_M_VAR = {};
                    end
                catch
                    warning('Matrices have different dimensions.')
                    clear M  
                    new_M_VAR = {};
                end
                
                % Concatenate old and new matrices
                try
                    matrix = cat(3,matrix,new_matrices);
                    M_VAR = vertcat(M_VAR, new_M_VAR);
                    %Updating the GUI 
                    fprintf('Number of (.mat) files selected are: %d \n', size(new_M_VAR,1));
                    set(disp_box,'String', M_VAR);
                    set(disp_box,'Value', []);
        
                    % Update the ROI x ROI x Subjects number
                    set(disp_str, 'String', strcat(num2str(size(matrix,2)), ' ROIs x',32, num2str(size(matrix,3)),' subjects'));
                    set(disp_str, 'ForegroundColor',[0.219, 0.341, 0.137]); 
                    clear M new_M_VAR
                catch
                    warning('Matrices have different number of ROIs.');
                    clear new_matrices new_M_VAR M
                end
            
            else
                warning('Selected *.mat file(s) consist(s) of multiple variables, please select *.mat files each containing only one variable.');
            end
            
        % If no files are selected
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

% Function to perform removal of files from list boxes
function file_remove(sel_var, M_VAR, disp_box,disp_str,case_maker)

   if isempty(sel_var) && isempty(M_VAR)
       warning('There are no files present to remove, please select .mat files to perform statistical analysis.');
       
   elseif isempty(sel_var) && ~isempty(M_VAR)
       warning('There are no selected matrices to remove from the list, please select matrices once again.');
   else
       M_VAR(sel_var,:) = [];
       holder = size(sel_var);
       fprintf('Number of (.mat) files removed are: %d \n', holder(2));
              
       set(disp_box,'Value', []);
       set(disp_box,'String', M_VAR);
       sel_var = {};
       
       if ~isempty(M_VAR)           
            % Update the ROI x ROI x Subjects counter under each case
            % Partial load the first file to update 
            matObj = matfile(M_VAR{1,:});
            S = whos(matObj);
            dims = S.size;
            roi_sub = [];
            subs = 0;
            
            % Update ROI x Subjects
            for i = 1:length(M_VAR)
                matObj = matfile(M_VAR{i,:});% Extract size of each variable per iteration
                temp = whos(matObj);
                temp_dim = temp.size;

                if length(temp_dim) == 2
                    subs = subs + 1;
                elseif length(temp_dim) == 3
                    subs = subs + temp.size(3);
                else
                    warning('Error, the files must have ROI x ROI x Subjects dimensions.');
                end

            end
            
            % Update GUI 
            roi_sub = [dims(1), subs];
            set(disp_str, 'String', strcat(num2str(roi_sub(1)), ' ROIs x',32, num2str(roi_sub(2)),' subjects'));
            set(disp_str, 'ForegroundColor',[0.219, 0.341, 0.137]); 

       end
   end
  if isempty(M_VAR)
        set(disp_str, 'String', '0 ROIs x 0 subjects');
        set(disp_str, 'ForegroundColor',[0.773, 0.353, 0.067]);
  end
   
  selection_0 = '';
  selection_1 = '';
  selection_2 = '';
  
  switch (case_maker)
      
      case 'one_samp_sel'
          M0 = M_VAR;
          if isempty(M0)
              M0 = {};
          end
          
      case 'left_samp_sel'
          M1 = M_VAR;
          if isempty(M1)
              M1 = {};
          end
          
      case 'right_samp_sel'
          M2 = M_VAR;
          if isempty(M2)
              M2 = {};
          end
  end
     
end


% Variable to store Live selection from lists
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
    if strcmp(contender, 'Paired t-test')
        
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
       
    elseif strcmp(contender, 'One-sample t-test')
        
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
        
        % Reset Varaibles
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
                
    elseif strcmp(contender, 'Two-sample t-test')
        
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
        clear matrices 
    end    
end

% Type of Parameter
function threshold_type(~,~)
    
    approach = (ST_THRES_POP.String{ST_THRES_POP.Value});
    
    if strcmp(approach, 'Uncorrected (Parametric)') || strcmp(approach, 'FDR (Parametric)') || strcmp(approach, 'Bonferroni (Parametric)') || strcmp(approach, 'NBS FWE (Non-Parametric)') || strcmp(approach, 'NBS TFCE (Non-Parametric)') 
        set(ST_PERM_TXT, 'enable', 'off');
        set(ST_PERM_VAL, 'enable', 'off', 'String', []);
        
    elseif strcmp(approach, 'Uncorrected (Non-Parametric)') || strcmp(approach, 'FDR (Non-Parametric)') 
        set(ST_PERM_TXT, 'enable', 'on');
        set(ST_PERM_VAL, 'enable', 'on', 'String', []);

    end
    
    if strcmp(approach, 'Uncorrected (Parametric)') || strcmp(approach, 'FDR (Parametric)') || strcmp(approach, 'Bonferroni (Parametric)') || strcmp(approach, 'Uncorrected (Non-Parametric)') || strcmp(approach, 'FDR (Non-Parametric)')
        set(ST_THRES_VAL_TXT, 'enable', 'off');
        set(ST_THRES_VAL_UNI, 'enable', 'off', 'String', []);

    elseif strcmp(approach, 'NBS FWE (Non-Parametric)') 
        set([ST_THRES_VAL_TXT,ST_PERM_TXT], 'enable', 'on');
        set([ST_THRES_VAL_UNI,ST_PERM_VAL], 'enable', 'on', 'String', []);

    elseif strcmp(approach, 'NBS TFCE (Non-Parametric)') 
        set([ST_THRES_VAL_TXT,ST_PERM_TXT], 'enable', 'off');
        set([ST_THRES_VAL_UNI,ST_PERM_VAL], 'enable', 'off', 'String', []);
    end
    
    if strcmp(approach, 'Uncorrected (Non-Parametric)') || strcmp(approach, 'FDR (Non-Parametric)') || strcmp(approach, 'NBS FWE (Non-Parametric)') || strcmp(approach, 'NBS TFCE (Non-Parametric)')
        set(RES_RUN,'enable', 'off');
        warning('Non-parametric analysis under development. Please wait for future updates.')
    else
        set(RES_RUN,'enable', 'on');
    end
    
end

% Run Test
function run(~,~)

    test_type = (ST_test_type.String{ST_test_type.Value}); % Type of Test - paried, one , two 
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if strcmp(test_type, 'Paired t-test')          
    
        if ~isempty(M1) && ~isempty(M2)
        
             % ROI size calculation for Set 1 
            matObj = matfile(M1{1,:});
            S1 = whos(matObj);
            dims_file_1 = S1.size;

            % ROI size calculation for Set 1 
            matObj = matfile(M2{1,:});
            S2 = whos(matObj);
            dims_file_2 = S2.size;

            % Checking if ROI x ROIs of either files are consistent
            if dims_file_1(1) == dims_file_2(1) && dims_file_1(2) == dims_file_2(2)

                % Compare the no of subjects across values (ROI x ROI x Subjects)
                if length(dims_file_1) == 2 && length(dims_file_2) == 2
                    
                    % Check if dimensions are equal
                    if length(M1) == length(M2)
                        flag_parameters = parameter_check(); % alpha, contrast

                        if flag_parameters == 1
                            flag_thres_parameters = perm_thres_check(); % permutations, threshold

                            if flag_thres_parameters == 1

                                % Update GUI 
                                set([ST_L0_CTR,ST_L1_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]); 
                                max{1} = matrices_1; max{2} = matrices_2;

                                % Calling test function from main TMFC dir
                                [thresholded,pval,tval,conval] = tmfc_ttest(max,str2num(ST_CONT_val.String),eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                if ~isempty(test_key(ST_THRES_POP.String{ST_THRES_POP.Value}))
                                    fprintf('Generating results...\n');
                                    tmfc_results_GUI(thresholded,pval,tval,conval,eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                end
                                clear thresholded pval tval conval max;
                            end
                        end
                    else
                        warning('The number of matrices are inconsistent to perform Paired t-test. \nMatrice(s) for RM1: %d \nMatrice(s) for RM2: %d \nThe number of matrices for two repeated measurements must be the same.', length(M1), length(M2));
                        set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]); % Update GUI 
                    end

                % If file consists of multiple subjects - check dimensions
                % One matrix versus multiple subjects
                elseif length(dims_file_1) == 2 && length(dims_file_2) == 3
                    
                    % Check dimensions of files
                    if length(M1) == dims_file_2(3)
                        flag_parameters = parameter_check(); % alpha, contrast check
                        
                        if flag_parameters == 1
                            flag_thres_parameters = perm_thres_check(); % permutations, threshold check
                            
                            if flag_thres_parameters == 1
                                
                                % Update GUI
                                set([ST_L0_CTR,ST_L1_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]);
                                max{1} = matrices_1; max{2} = matrices_2;
                                
                                % Calling test function from main TMFC dir
                                [thresholded,pval,tval,conval] = tmfc_ttest(max,str2num(ST_CONT_val.String),eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                if ~isempty(test_key(ST_THRES_POP.String{ST_THRES_POP.Value}))
                                    fprintf('Generating results...\n');
                                    tmfc_results_GUI(thresholded,pval,tval,conval,eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                end
                                clear thresholded pval tval conval max;
                            end
                        end
                    else
                        warning('The number of matrices are inconsistent to perform Paired t-test. \nMatrice(s) for RM1: %d \nMatrice(s) for RM2: %d \nThe number of matrices for two repeated measurements must be the same.', length(M1), dims_file_2(3));
                        set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]); % Update GUI 
                    end

                % If file consists of multiple subjects - check dimensions
                % Multiple subjects versus one matrix
                elseif length(dims_file_1) == 3 && length(dims_file_2) == 2
                    
                    % Compare the no of subjects across values (ROI x ROI x Subjects)
                    if dims_file_1(3) == length(M2)
                        flag_parameters = parameter_check(); % alpha, contrast check
                        
                        if flag_parameters == 1
                            flag_thres_parameters = perm_thres_check(); % permutations, threshold check
                            
                            if flag_thres_parameters == 1
                                
                                % Update GUI 
                                set([ST_L0_CTR,ST_L1_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]); 
                                max{1} = matrices_1; max{2} = matrices_2;
                                
                                [thresholded,pval,tval,conval] = tmfc_ttest(max,str2num(ST_CONT_val.String),eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                if ~isempty(test_key(ST_THRES_POP.String{ST_THRES_POP.Value}))
                                    fprintf('Generating results...\n');
                                    tmfc_results_GUI(thresholded,pval,tval,conval,eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                end
                                clear thresholded pval tval conval max;
                            end
                        end
                    else
                        warning('The number of matrices are inconsistent to perform Paired t-test. \nMatrice(s) for RM1: %d \nMatrice(s) for RM2: %d \nThe number of matrices for two repeated measurements must be the same.', dims_file_1(3), length(M2));
                        set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]); % Update GUI 
                    end

                % if file consists of multiple subjects - check dimensions
                % check multiple subjects versus multiple subjets
                elseif length(dims_file_1) == 3 && length(dims_file_2) == 3
                    
                    if dims_file_1(3) == dims_file_2(3)
                        flag_parameters = parameter_check(); % alpha, contrast check
                        
                        if flag_parameters == 1
                            flag_thres_parameters = perm_thres_check(); % permutations, threshold check
                            
                            % Update GUI 
                            if flag_thres_parameters == 1
                                set([ST_L0_CTR,ST_L1_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]); 
                                max{1} = matrices_1; max{2} = matrices_2;
                                
                                % Calling test function from main TMFC dir
                                [thresholded,pval,tval,conval] = tmfc_ttest(max,str2num(ST_CONT_val.String),eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                if ~isempty(test_key(ST_THRES_POP.String{ST_THRES_POP.Value}))
                                    fprintf('Generating results...\n');
                                    tmfc_results_GUI(thresholded,pval,tval,conval,eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                                end
                                clear thresholded pval tval conval max;
                            end
                        end
                    else
                        warning('The number of matrices are inconsistent to perform Paired t-test. \nMatrice(s) for RM1: %d \nMatrice(s) for RM2: %d \nThe number of matrices for two repeated measurements must be the same.', dims_file_1(3), dims_file_2(3));
                        set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]); % Update GUI 
                    end

                else 
                    warning('Error, the files must have ROI x ROI x Subjects dimensions.');
                end
                
            else
               warning('The number of ROI x ROIs between the two lists are inconsistent, please select matrices with consistent number of ROIs.');
               set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]); % Update GUI 

            end

        elseif ~isempty(M1) && isempty(M2)
            warning('Please select SECOND set of Matrices to perform Paired t-test.');

        elseif isempty(M1) && ~isempty(M2)
            warning('Please select FIRST set of Matrices to perform Paired t-test.');
            
        else
            warning('Please select matrices files to perform Paired t-test.');
        end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    elseif strcmp(test_type, 'One-sample t-test')
        
        % Check for matrix
        if ~isempty(M0)
            flag_parameters = parameter_check(); % Alpha, contrast check

            if flag_parameters == 1
                flag_thres_parameters = perm_thres_check(); % permutations, threshold check

                if flag_thres_parameters == 1
                    % Update GUI 
                    set([ST_L0_CTR,ST_L1_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]); 

                    % Calling test function from main TMFC dir
                    [thresholded,pval,tval,conval] = tmfc_ttest(matrices_0, str2num(ST_CONT_val.String),eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));                   
                    if ~isempty(test_key(ST_THRES_POP.String{ST_THRES_POP.Value}))
                        fprintf('Generating results...\n');
                        tmfc_results_GUI(thresholded,pval,tval,conval,eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                    end
                    clear thresholded pval tval conval;
                end
            end
            
        else
            warning('Please select matrices files to perform One-sample t-test.');
        end
        
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    elseif strcmp(test_type, 'Two-sample t-test')

        % Check for matrix
        if ~isempty(M1) && ~isempty(M2)
        
            % ROI size calculation for Set 1 
            matObj = matfile(M1{1,:});
            S1 = whos(matObj);
            dims_file_1 = S1.size;
                        
            % ROI size calculation for Set 1 
            matObj = matfile(M2{1,:});
            S2 = whos(matObj);
            dims_file_2 = S2.size;

            % Checking Dimensions of Files
            if dims_file_1(1) == dims_file_2(1) && dims_file_1(2) == dims_file_2(2)
                flag_parameters = parameter_check(); % alpha, contrast check
                
                if flag_parameters == 1
                    flag_thres_parameters = perm_thres_check(); % permutations, threshold check
                    
                    if flag_thres_parameters == 1
                        % Update GUI 
                        set([ST_L0_CTR,ST_L1_CTR], 'ForegroundColor',[0.219, 0.341, 0.137]); 
                        max{1} = matrices_1; max{2} = matrices_2;
                        
                        % Calling test function from main TMFC dir
                        [thresholded,pval,tval,conval] = tmfc_ttest2(max,str2num(ST_CONT_val.String),eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                        if ~isempty(test_key(ST_THRES_POP.String{ST_THRES_POP.Value}))
                            fprintf('Generating results...\n');
                            tmfc_results_GUI(thresholded,pval,tval,conval,eval(get(ST_ALP_val, 'String')),test_key(ST_THRES_POP.String{ST_THRES_POP.Value}));
                        end
                        clear thresholded pval tval conval max;
                    end
                end
            else
               warning('The number of ROI x ROIs between the two lists are inconsistent, please select matrices with consistent number of ROIs.');
               set([ST_L1_CTR,ST_L2_CTR], 'ForegroundColor',[0.773, 0.353, 0.067]); % Update GUI 
            end

        elseif ~isempty(M1) && isempty(M2)
            warning('Please select SECOND set of Matrices to perform Paired t-test.');
            
        elseif isempty(M1) && ~isempty(M2)
            warning('Please select FIRST set of Matrices to perform Paired t-test.'); 
            
        else
            warning('Please select matrices files to perform Paired t-test.');
            
        end 
        
    else
        warning('Error, the files must have ROI x ROI x Subjects dimensions.');
    end 
    
end

% Function to check Contrast, alpha logical numeric bounds 
function flag = parameter_check(~,~)
        
    flag = 0;
    test_type = (ST_test_type.String{ST_test_type.Value}); % Type of Test - paried, one , two 
    val_contrast = str2num(ST_CONT_val.String);        % Contrast
    try
        val_alpha = eval(get(ST_ALP_val, 'String'));    % Alpha
    catch
        val_alpha = NaN;
    end
    
     if isempty(val_contrast)
         warning('Please enter numeric values for the contrast.');
     else
         if strcmp(test_type, 'Paired t-test') || strcmp(test_type, 'Two-sample t-test')
             if length(val_contrast) < 2
                 warning('Please enter TWO contrast values for computation.');
             elseif length(val_contrast) > 2
                 warning('Number of contrast values cannot be greater than TWO, please re-enter contrast values.');
             else
                  if isnan(val_alpha)
                     warning('Please enter a numeric alpha value for computation.');
                  else
                    if val_alpha > 1 || val_alpha < 0
                        warning('Please re-enter alpha value between [0.0, 1.0].');
                    else
                        flag = 1;
                    end
                 end
            end
        elseif strcmp(test_type, 'One-sample t-test')
            if length(val_contrast) >=2
                warning('Number of contrast values cannot exceed ONE, please re-enter contrast value for computation.');
            else
                 if isnan(val_alpha)
                    warning('Please enter alpha value for computation.');
                 else
                    if (val_alpha > 1) || (val_alpha < 0)
                        warning('Please re-enter alpha value between [0.0, 1.0].');
                    else
                        flag = 1;
                    end
                 end
                
            end
        end
    end

end 

% Function to check permutation & threshold logical numeric bounds
function flag = perm_thres_check(~,~)

    flag = 0;
    thresh_type = (ST_THRES_POP.String{ST_THRES_POP.Value});

     if strcmp(thresh_type, 'Uncorrected (Non-Parametric)') || strcmp(thresh_type, 'FDR (Non-Parametric)') 
         
         % Check if permutations is a number and not a floating point value
         permutation_value = str2num(ST_PERM_VAL.String);
         if ~isempty(permutation_value) && permutation_value > 0
             warning('Non-parametric analysis under development.');
             flag = 1;
         elseif permutation_value <= 0
             warning('Please enter a postive numeric value for Number of Permutations.');
         else
             warning('Please enter a numeric value for Number of Permutations.');
         end         
         
     elseif strcmp(thresh_type, 'NBS FWE (Non-Parametric)') || strcmp(thresh_type, 'NBS TFCE (Non-Parametric)') 
         permutation_value = str2num(ST_PERM_VAL.String);
         threshold_value = str2double(ST_THRES_VAL_UNI.String);
         if ~isnan(threshold_value) && threshold_value > 0 && threshold_value <=1.0
             warning('Non-parametric analysis under development.');
             flag = 1;
         elseif threshold_value <= 0 || threshold_value > 1.0
             warning('Please enter a Primary Threshold value between [0.0, 1.0] for computation.');
         else
             warning('Please enter a Primary Threshold value for computation.');
         end       
     elseif strcmp(thresh_type, 'Uncorrected (Parametric)') || strcmp(thresh_type, 'FDR (Parametric)') || strcmp(thresh_type, 'Bonferroni (Parametric)' )
         flag = 1;
     end
end

% Wait for GUI to complete/close to continue using TMFC
if ~isempty(findobj('Tag', 'TMFC_GUI'))
    uiwait();
end

end
%%
% Function to select (.mat) files from the user via spm_select
function list_sel = select_mat_file(~,~)  
    files = spm_select(inf,'.mat','Select matrices for computation',{},pwd,'.');
    list_sel = {};
    list_sel = cellstr(files);
end

%%
% Function to check if the selecte (.mat) files consists of multiple
% variables - Returns Binary Flag where, 
% 0 = no multiple variables in selected (.mat) files 
% 1 = Multiple variables EXIST in selected (.mat) files 
function flag = multi_var_check(input_file)
    
    main_size = size(input_file);    % Store size
    list_var = [];                   % variable to store list of multiple vars
    ctr = 1;                         % Counter the number of files present
    flag = 0;                        % Binary Flag to indicate status of multiple vars
    
    if ~isempty(input_file{1})
        % Loop to iterate through all possible (.mat) files 
        for i = 1:main_size(1)

            var = who('-file', input_file{i,:});   % Listing the variable into temp Workspace - Cell Datatype
            %var = who('-file', D(i,:));  % Listing the variable into temp Workspace - Standalone Datatype
            A = size(var);

            % If there exists files with multiple variables within, then disp
            if A(1) > 1
                fprintf('Multiple variables present in the file :%s \n', input_file{i,:})
                list_var(ctr) = i;
                ctr = ctr+1;
                flag = 1;
            end        
        end
    else
        flag = -1;
    end
    
end

%%
% Function to check selected files are of same dimension or not
% i.e. 2D = ROI x ROI 
% i.e. 3D = ROI x ROI x Subjects

% The selection of 2D or 3D is based on the first selected file i.e. if the 
% first file is of 2D dimensions then the code will check if all remaining 
% files are of 2D dimensions else it will not select the files.

% Similarly if a 3D file is selected as the first file, then it will check 
% if all files are in 3D format. 

% Function returns Binary Flag 
% if flag = 1, selected files HAVE INCONSISTENT dimensions
% if flag = 0, selected files have CONSISTENT dimensions
function flag = dimension_check(input_files)

    ctr = 1;                          % Counter to store number of files 
    file_address = [];                % variable to store file addresses
    dimension = size(input_files);    % Size of input list of files
    flag = 0;                         % Final Flag to indicate accpet or reject dimension check
    
    % Run the loop for given list of selected files
    for i = 1:dimension(1)
        
        % Intializing the comparison format (2D or 3D)
        % The file to compare against the rest of the file
        if i == 1
            
            % Format to load OBJECT parts of the file without data
            matObj = matfile(input_files{1,:});
            S = whos(matObj);                    % Get characteristics of the data
            file_1_size = S.size;                % Extract Size
        else
            
        % For all other iterations, directly extract the dimension of files
            
            % Extraction of dimension of the variables from (.mat) files
            matObj = matfile(input_files{i,:});
            S = whos(matObj);
            file_2_size = S.size;
            
            
        % Checking if the dimensions of FIRST file is same against the rest
            if length(file_1_size) ~= length(file_2_size)
                
                % Print the Dimension & Path to the files that are inconsistent
                fprintf('\n The Dimensions of following files are not equal, please re-select files with consistent dimensions:')
                fprintf('\nFile 1: %s', input_files{1,:});
                fprintf('\nDimensions: %s ', num2str(file_1_size));
                
                fprintf('\nFile %d: %s has been excluded from the selection',i, input_files{i,:});
                fprintf('\nDimensions: %s \n', num2str(file_2_size));
                file_address(ctr) = i;
                ctr = ctr+1;   
                flag = 1;                
            end            
        end
    end
end

%% Function to convert internal labelling format to tmfc_ttest() labelling format
function small_string = test_key(input_string)

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