function [conditions] = tmfc_conditions_GUI(SPM_path,input_case)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Opens a GUI window for gPPI analysis. Allows to choose conditions of
% interest for gPPI analysis.
% 
% FORMAT [conditions] = tmfc_conditions_GUI(SPM_path, input_case)
%   SPM_path          - Path to individual subject SPM.mat file
%   input_case        - Select conditions of interest, three cases:
%                       [1 = Specify ROI set, 2 = gPPI, 3 = LSS]
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

% Get all conditions from the SPM.mat file
all_cond = generate_conditions(SPM_path);

% Check if SPM.mat is not empty
if isempty(all_cond)
    error('Selected SPM.mat file is empty or invalid.');
else
    % Select conditions via GUI
    conditions = conditions_GUI(all_cond, input_case);
end
end


%% =========[ Select conditions of interest for analysis via GUI ]=========
function [conditions] = conditions_GUI(all_cond, input_case)

    cond_L1 = {};       % Variable to store all conditions in GUI 
    cond_L2 = {};       % Variable to store selected conditions in GUI 
    conditions_MW_SE1 = {};   % Variable to store the selected list of conditions in BOX 1(as INDEX)
    conditions_MW_SE2 = {};   % Variable to store the selected list of conditions in BOX 2(as INDEX)

    switch(input_case)
        case 1
            MW_string = 'Select ROIs: Select conditions';
            HW_string = 'Select ROIs: Help';
        case 2
            MW_string = 'gPPI: Select conditions';
            HW_string = 'gPPI: Help';  
        case 3
            MW_string = 'LSS: Select conditions';
            HW_string = 'LSS: Help';
            all_cond([all_cond(:).pmod]~=1)=[];
    end


    for iCond = 1:length(all_cond)
        cond_L1 = vertcat(cond_L1, all_cond(iCond).list_name);        
    end
    
    conditions_MW = figure('Name', MW_string, 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.45 0.25 0.22 0.56],'MenuBar', 'none','ToolBar', 'none','color','w','Resize','on','WindowStyle','modal','CloseRequestFcn', @MW_exit);

    cond_MW_S1  = uicontrol(conditions_MW,'Style','text','String', 'Select conditions of interest','Units', 'normalized', 'Position',[0.270 0.93 0.460 0.05],'fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w');
    cond_MW_S2  = uicontrol(conditions_MW,'Style','text','String', 'All conditions:','Units', 'normalized', 'Position',[0.045 0.88 0.450 0.05],'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w');
    cond_MW_S3  = uicontrol(conditions_MW,'Style','text','String', 'Conditions of interest:','Units', 'normalized', 'Position',[0.045 0.425 0.450 0.05],'HorizontalAlignment', 'left','fontunits','normalized', 'fontSize', 0.50,'backgroundcolor','w');
    cond_MW_LB1 = uicontrol(conditions_MW , 'Style', 'listbox', 'String', cond_L1,'Max', 1000,'Units', 'normalized', 'Position',[0.045 0.59 0.900 0.300],'fontunits','points', 'fontSize', 11,'Value', [],'callback', @MW_LB1_select);
    cond_MW_LB2 = uicontrol(conditions_MW , 'Style', 'listbox', 'String', cond_L2,'Max', 1000,'Units', 'normalized', 'Position',[0.045 0.135 0.900 0.300],'fontunits','points', 'fontSize', 11,'Value', [],'callback', @MW_LB2_select);

    cond_MW_add = uicontrol(conditions_MW,'Style','pushbutton','String', 'Add selected','Units', 'normalized','Position',[0.045 0.50 0.270 0.065],'fontunits','normalized', 'fontSize', 0.32,'callback', @MW_add);
    cond_MW_add_all = uicontrol(conditions_MW,'Style','pushbutton','String', 'Add all','Units', 'normalized','Position',[0.360 0.50 0.270 0.065],'fontunits','normalized', 'fontSize', 0.32,'callback', @MW_add_all); 
    cond_MW_help = uicontrol(conditions_MW,'Style','pushbutton','String', 'Help','Units', 'normalized','Position',[0.680 0.50 0.270 0.065],'fontunits','normalized', 'fontSize', 0.32,'callback', @MW_help); 

    cond_MW_confirm = uicontrol(conditions_MW,'Style','pushbutton','String', 'OK','Units', 'normalized','Position',[0.045 0.05 0.270 0.065],'fontunits','normalized', 'fontSize', 0.32,'callback', @MW_confirm); 
    cond_MW_remove = uicontrol(conditions_MW,'Style','pushbutton','String', 'Remove selected','Units', 'normalized','Position',[0.360 0.05 0.270 0.065],'fontunits','normalized', 'fontSize', 0.32,'callback', @MW_remove); 
    cond_MW_remove_all = uicontrol(conditions_MW,'Style','pushbutton','String', 'Remove all','Units', 'normalized','Position',[0.680 0.05 0.270 0.065],'fontunits','normalized', 'fontSize', 0.32,'callback', @MW_remove_all);
    movegui(conditions_MW,'center');

    %----------------------------------------------------------------------
    % Function to close conditions GUI without selection of conditions
    %----------------------------------------------------------------------
    function MW_exit(~,~)
        delete(conditions_MW);
        conditions = NaN;
    end

    %----------------------------------------------------------------------
    function MW_LB1_select(~,~)
        index = get(cond_MW_LB1, 'Value'); 
        conditions_MW_SE1 = index;      
    end

    %----------------------------------------------------------------------
    function MW_LB2_select(~,~)
        index = get(cond_MW_LB2, 'Value');  
        conditions_MW_SE2 = index;             
    end

    %----------------------------------------------------------------------
    % Function to add single condition
    %----------------------------------------------------------------------
    function MW_add(~,~)
        
        if isempty(conditions_MW_SE1)
            warning('No conditions selected.');
        else
            len_exist = length(cond_L2);     
            new_conds = {};                  

            % Based on the selection add variables to a selected list
            new_conds = vertcat(new_conds,cond_L1(conditions_MW_SE1)); 

            % Addition & extraction of unique selected conditions
            cond_L2 = vertcat(cond_L2,new_conds);
            new_cond_count = length(unique(cond_L2)) - len_exist;
            cond_L2 = unique(cond_L2);

            % Check if newly selected conditions have been added
            if new_cond_count == 0
                warning('Newly selected conditions are already present in the list, no new conditions added.');
            else
                fprintf('Conditions selected: %d. \n', new_cond_count(1)); 
                cond_L2 = sort_selected_conditions(cond_L2,all_cond);
            end 

            % Set sorted list of conditions into GUI
            set(cond_MW_LB2,'String',cond_L2);
        end
    end

    %----------------------------------------------------------------------
    % Function to add all conditions 
    %----------------------------------------------------------------------
    function MW_add_all(~,~) 

        if length(cond_L2) == length(cond_L1)
            warning('All conditions are already selected.');
        else
        	len_exist = length(cond_L2);
            new_conds = {};                                             
            new_conds = vertcat(new_conds,cond_L1);             

            % Addition & extraction of unique selected conditions
            cond_L2 = vertcat(cond_L2,new_conds);
            new_cond_count = length(unique(cond_L2)) - len_exist;
            cond_L2 = unique(cond_L2);

            % Check if newly selected conditions have been added
            if new_cond_count == 0
                warning('Newly selected conditions are already present in the list, no new conditions added.');
            else
                fprintf('New conditions selected: %d. \n', new_cond_count(1)); 
                cond_L2 = sort_selected_conditions(cond_L2,all_cond);
            end 

            % Set sorted list of conditions into GUI
            set(cond_MW_LB2,'String',cond_L2);
        end
    end

    %----------------------------------------------------------------------
    % Function to export list of conditions
    %----------------------------------------------------------------------
    function MW_confirm(~,~)

    	if isempty(cond_L2)
        	warning('Please select conditions.');
        else
            cond = struct;
            n_cond = 1;
            for jCond = 1:length(all_cond)
               for kCond = 1:length(cond_L2)
                   match = strcmp(all_cond(jCond).list_name,cond_L2(kCond));
                   if match == 1
                       cond(n_cond).sess = all_cond(jCond).sess;
                       cond(n_cond).number = all_cond(jCond).number;
                       cond(n_cond).pmod = all_cond(jCond).pmod;
                       cond(n_cond).name = all_cond(jCond).name;
                       cond(n_cond).list_name = all_cond(jCond).list_name;
                       cond(n_cond).file_name = all_cond(jCond).file_name;
                       n_cond = n_cond + 1;
                   end
               end
            end
            delete(conditions_MW);
            disp(strcat(num2str(length(cond_L2)),' conditions successfully selected.'));
            conditions = cond;
            clear cond n_cond jCond kCond match
        end
    end

    %----------------------------------------------------------------------
    % Function to remove single condition
    %----------------------------------------------------------------------
    function MW_remove(~,~)
        if isempty(cond_L2)
        	warning('No conditions present to remove.');
        elseif isempty(conditions_MW_SE2)
        	warning('No conditions selected to remove.');
        else
            cond_L2(conditions_MW_SE2,:) = [];
            fprintf('Number of conditions removed: %d. \n', length(conditions_MW_SE2));
            set(cond_MW_LB2, 'Value', []);
            set(cond_MW_LB2, 'String', cond_L2);
            conditions_MW_SE2 = {};
        end
    end

    %----------------------------------------------------------------------
    % Function to remove all conditions
    %----------------------------------------------------------------------
    function MW_remove_all(~,~) 
        if isempty(cond_L2)
            warning('No conditions present to remove.');
        else
            cond_L2 = {};                                             
            set(cond_MW_LB2, 'String', []);
            conditions_MW_SE2 = {};
            warning('All selected conditions have been removed.');
        end
    end

    %----------------------------------------------------------------------
    % Help window for conditions selection
    %----------------------------------------------------------------------
    function MW_help(~,~)

        cond_HW = figure('Name', HW_string, 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.67 0.31 0.22 0.50],'MenuBar', 'none','ToolBar', 'none','color','w','Resize','off', 'WindowStyle', 'Modal');

        if strcmp(HW_string, 'Select ROIs: Help')
            string_info = {'Suppose you have two separate sessions.','Both sessions contains task regressors for "Cond A", "Cond B" and "Errors"','',...
            'If you are only interested in "Cond A" and "Cond B" comparison, the following conditions should be selected:','','Cond A (Sess1)',...
            'Cond B (Sess1)','Cond A (Sess2)','Cond B (Sess2)','','For all selected conditions of interest, the TMFC toolbox will calculate the omnibus F-contrast.',...
            '','The center of the moving sphere will be moved to the local maximum inside a fixed sphere of larger radius.','',...
            'Local maxima are determined using the omnibus F-contrast for the selected conditions of interest.'};
            
        elseif strcmp(HW_string, 'gPPI: Help')
            string_info = {'Suppose you have two separate sessions.','','Both sessions contains task regressors for', '"Cond A", "Cond B" and "Errors"', '','If you are only interested in "Cond A" and "Cond B" comparison, the following conditions should be selected:',...
            '','1)  Cond A (Sess1)','2)  Cond B (Sess1)','3)  Cond A (Sess2)','4)  Cond B (Sess2)','','For all selected conditions of interest, the TMFC toolbox will create psycho-physiological (PPI) regressors. Thus, for each condition of interest, the generalized PPI (gPPI) model will contain two regressors: (1) psychological regressor and (2) PPI regressor.'...
            '','For trials of no interest (here, "Errors"), the gPPI model will contain only the psychological regressor.'}; 
        
        else
            string_info = {'Suppose you have two separate sessions.','','Both sessions contains task regressors for', '"Cond A", "Cond B" and "Errors"', '','If you are only interested in "Cond A" and "Cond B" comparison, the following conditions should be selected:',...
            '','1)  Cond A (Sess1)','2)  Cond B (Sess1)','3)  Cond A (Sess2)','4)  Cond B (Sess2)','','For all selected conditions of interest, the TMFC toolbox will calculate individual trial beta-images using Least-Squares Separate (LSS) approach.',...
            '','For each individual trial (event), the LSS approach estimates a separate general linear model (GLM) with two regressors. The first regressor models the expected BOLD response to the current trial of interest, and the second (nuisance) regressor models the BOLD response to all other trials (of interest and no interest).',...
            '','For trials of no interest (here, "Errors"), separate GLMs will not be estimated. Trials of no interest will be used only for the second (nuisance) regressor.'};
        end
    
        % Create if-else chain, compare if HW_string == ROI, gPPI, LSS and
        % then assign the repsective help string, change window size if
        % required. 

        cond_HW_S1 = uicontrol(cond_HW,'Style','text','String',string_info ,'Units', 'normalized', 'Position', [0.05 0.12 0.89 0.85], 'HorizontalAlignment', 'left','backgroundcolor','w','fontunits','normalized', 'fontSize', 0.0301);
        cond_HW_OK = uicontrol(cond_HW,'Style','pushbutton','String', 'OK','Units', 'normalized', 'Position', [0.34 0.06 0.30 0.06],'callback', @cond_HW_close,'fontunits','normalized', 'fontSize', 0.35);
        movegui(cond_HW,'center');

        function cond_HW_close(~,~)
        	close(cond_HW);
        end
    end

    uiwait(conditions_MW);
end


%% ===========[ Function to get information about conditions ]=============
function [cond_list] = generate_conditions(SPM_path)
    cond_list = {}; 
    try
        load(SPM_path);
        kCond = 1;
        for iSess = 1:length(SPM.Sess)
            for jCond = 1:length(SPM.Sess(iSess).U)
                for kPmod = 1:length(SPM.Sess(iSess).U(jCond).name)
                    cond_list(kCond).sess = iSess;
                    cond_list(kCond).number = jCond;
                    cond_list(kCond).pmod = kPmod;
                    cond_list(kCond).name = char(SPM.Sess(iSess).U(jCond).name(kPmod));
                    cond_list(kCond).list_name = [char(SPM.Sess(iSess).U(jCond).name(kPmod)) ' (Sess' num2str(iSess) ', Cond' num2str(jCond) ')'];
                    cond_list(kCond).file_name = ['[Sess_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
                                                 regexprep(char(SPM.Sess(iSess).U(jCond).name(kPmod)),' ','_') ']'];
                    kCond = kCond + 1;
                end
            end 
        end
    catch 
        disp('Selected SPM.mat file does not exist or invalid.');
        cond_list = {};
    end
end

%% ========[ Function to perform sorting of selected conditions ]==========
function [sorted_list] = sort_selected_conditions(selected_cond,all_cond)

    temp = {};
    sort_index = 1;
    for iCond = 1:length(selected_cond)
        for jCond = 1:length(all_cond)
            if strcmp(selected_cond(iCond),all_cond(jCond).list_name)
                if sort_index == 1
                    temp = all_cond(jCond);
                    sort_index = sort_index + 1;
                else 
                    temp(sort_index) = all_cond(jCond);
                    sort_index = sort_index + 1;
                end
            end
        end
    end

    [~,index] = sortrows([temp.sess; temp.number]');
    reindexed_list = temp(index); 

    sorted_list = {};
    for iCond = 1:length(reindexed_list) 
        sorted_list = vertcat(sorted_list, reindexed_list(iCond).list_name);
    end

    clear index temp sort_index iCond jCond reindexed_list
end