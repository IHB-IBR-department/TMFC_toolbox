function [ROI_set] = tmfc_select_ROIs_GUI(tmfc)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Opens a GUI window for selecting ROI masks. Creates group mean binary 
% mask based on 1st-level masks (see SPM.VM) and applies it to all selected
% ROIs. Empty ROIs will be removed. Masked ROIs will be limited to only
% voxels which have data for all subjects. The dimensions, orientation, and
% voxel sizes of the masked ROI images will be adjusted according to the
% group mean binary mask.
%
% FORMAT [ROI_set] = tmfc_select_ROIs_GUI(tmfc)
%
% Input:
%   tmfc.subjects.path     - Paths to individual SPM.mat files
%   tmfc.project_path      - The path where all results will be saved
%
% Output:
%   ROI_Set                - Structure with information about selected ROIs
%
% =========================================================================
%
% Copyright (C) 2024 Ruslan Masharipov
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

% Initial checks
if ~isfield(tmfc,'subjects')
    error('Select subjects.');
elseif strcmp(tmfc.subjects(1).path, '')
    error('Select subjects.');
elseif ~exist(tmfc.subjects(1).path,'file')
    error('SPM.mat file for the first subject does not exist.')
elseif ~isfield(tmfc,'project_path')
    error('Select TMFC project folder.');
end

% Specify ROI set name
[ROI_set_name] = ROI_set_name_GUI();

% Specify ROI set structure   
if ~strcmp(ROI_set_name,'')
    ROI_set = ROI_set_generation(ROI_set_name);  
else
    warning('ROIs not selected.');
    ROI_set = [];
end

% -------------------------------------------------------------------------
% Select ROIs, create ROI masks and remove heavily cropped ROIs
function [ROI_set] = ROI_set_generation(ROI_set_name)
    
    ROI_set.set_name = ROI_set_name;
    
    % Select ROIs
    try
        [ROI_paths] = spm_select(inf,'any','Select ROI masks',{},pwd);
        for iROI = 1:size(ROI_paths,1)
            [~, ROI_set.ROIs(iROI).name, ~] = fileparts(deblank(ROI_paths(iROI,:)));
            ROI_set.ROIs(iROI).path = deblank(ROI_paths(iROI,:));
            ROI_set.ROIs(iROI).path_masked = fullfile(tmfc.project_path,'ROI_sets',ROI_set_name,'Masked_ROIs',[ROI_set.ROIs(iROI).name '_masked.nii']);
        end
    catch
        warning('ROIs not selected.');
    end
    
    % Calculate ROI size and remove heavily cropped ROIs
    if ~isempty(ROI_paths)
        
        % Clear & create 'Masked_ROIs' folder
        if isdir(fullfile(tmfc.project_path,'ROI_sets',ROI_set_name))
            rmdir(fullfile(tmfc.project_path,'ROI_sets',ROI_set_name),'s');
        end
        
        if ~isdir(fullfile(tmfc.project_path,'ROI_sets',ROI_set_name,'Masked_ROIs'))
            mkdir(fullfile(tmfc.project_path,'ROI_sets',ROI_set_name,'Masked_ROIs'));
        end
        
        % Create group mean binary mask
        for iSub = 1:length(tmfc.subjects)
            sub_mask{iSub,1} = [tmfc.subjects(iSub).path(1:end-7) 'mask.nii'];
        end
        group_mask = fullfile(tmfc.project_path,'ROI_sets',ROI_set_name,'Masked_ROIs','Group_mask.nii');
        
        if length(tmfc.subjects) == 1
            copyfile(sub_mask{1,1},group_mask);
        else
            spm_imcalc(sub_mask,group_mask,'prod(X)',{1,0,1,2});
        end
        
        % Calculate ROI size before masking
        w = waitbar(0,'Please wait...','Name','Calculating raw ROI sizes');
        group_mask = spm_vol(group_mask);
        nROI = numel(ROI_set.ROIs);
        for iROI = 1:nROI
            ROI_mask = spm_vol(ROI_set.ROIs(iROI).path);
            Y = zeros(group_mask.dim(1:3));
            % Loop through slices
            for p = 1:group_mask.dim(3)
                % Adjust dimensions, orientation, and voxel sizes to group mask
                B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
                X = zeros(1,prod(group_mask.dim(1:2)));
                M = inv(B * inv(group_mask.mat) * ROI_mask.mat);
                d = spm_slice_vol(ROI_mask, M, group_mask.dim(1:2), 1);
                d(isnan(d)) = 0;
                X(1,:) = d(:)';
                Y(:,:,p) = reshape(X,group_mask.dim(1:2));
            end
            % Raw ROI size (in voxels)
            ROI_set.ROIs(iROI).raw_size = nnz(Y);
            try
                waitbar(iROI/nROI,w,['ROI No ' num2str(iROI,'%.f')]);
            end
        end
        
        try
            close(w);
        end
        
        % Mask the ROI images by the goup mean binary mask
        w = waitbar(0,'Please wait...','Name','Masking ROIs by group mean mask');
        input_images{1,1} = group_mask.fname;
        for iROI = 1:nROI
            input_images{2,1} = ROI_set.ROIs(iROI).path;
            ROI_mask = ROI_set.ROIs(iROI).path_masked;
            spm_imcalc(input_images,ROI_mask,'(i1>0).*(i2>0)',{0,0,1,2});
            try
                waitbar(iROI/nROI,w,['ROI No ' num2str(iROI,'%.f')]);
            end
        end
        
        try
            close(w)
        end
        
        % Calculate ROI size after masking
        w = waitbar(0,'Please wait...','Name','Calculating masked ROI sizes');
        for iROI = 1:nROI
            ROI_set.ROIs(iROI).masked_size = nnz(spm_read_vols(spm_vol(ROI_set.ROIs(iROI).path_masked)));
            ROI_set.ROIs(iROI).masked_size_percents = 100*ROI_set.ROIs(iROI).masked_size/ROI_set.ROIs(iROI).raw_size;
            try
                waitbar(iROI/nROI,w,['ROI No ' num2str(iROI,'%.f')]);
            end
        end
        
        try
            close(w)
        end
        
        % Check for empty ROIs
        empty_ROI_list = {};
        empty_ROI_index = 1;
        for iROI = 1:length(ROI_set.ROIs)
            if ROI_set.ROIs(iROI).masked_size_percents == 0
                empty_ROI_list{empty_ROI_index,1} = iROI;
                empty_ROI_list{empty_ROI_index,2} = ROI_set.ROIs(iROI).name;
                empty_ROI_index = empty_ROI_index + 1;
            end
        end
        
        % GUI interface for removing heavily cropped ROIs
        if ~isempty(empty_ROI_list)
            disp_empty_ROI_list = {};
            for iROI = 1:size(empty_ROI_list,1)
                ROI_string = horzcat('No ',num2str(empty_ROI_list{iROI,1}),': ',empty_ROI_list{iROI,2});
                disp_empty_ROI_list = vertcat(disp_empty_ROI_list, ROI_string);
            end
            
            % GUI to show empty ROIs
            ROI_remove_empty_GUI(disp_empty_ROI_list);
            
            % Removing the empty ROIs
            ROI_index = 0;
            for iROI = 1:size(empty_ROI_list,1)
                ROI_set.ROIs(empty_ROI_list{iROI,1}-ROI_index) = [];
                ROI_index = ROI_index +1;
            end
            
            % Ask user to remove heavily cropped ROIs
            if isempty(ROI_set.ROIs)
                warning('All ROIs are empty. Select different ROIs.');
            else
                ROI_set = ROI_remove_crop(ROI_set);
            end
        else
            ROI_set = ROI_remove_crop(ROI_set);
        end
    else
        warning('ROIs not selected.');
    end
    
    if ~isfield(ROI_set,'set_name') || ~isfield(ROI_set,'ROIs')
        ROI_set = [];
    end
end   
end

%% ====================[ Specify ROI set name GUI ]========================
function [ROI_set_name] = ROI_set_name_GUI(~,~)
    
    ROI_set_name = '';

    ROI_set_name_MW = figure('Name', 'Select ROIs', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.62 0.50 0.16 0.14],'Resize','off','color','w','MenuBar', 'none','ToolBar', 'none','WindowStyle', 'modal','CloseRequestFcn', @ROI_set_name_MW_EXIT);    
    ROI_set_name_MW_S = uicontrol(ROI_set_name_MW,'Style','text','String', 'Enter a name for the ROI set','Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.40,'backgroundcolor',get(ROI_set_name_MW,'color'),'Position',[0.14 0.60 0.700 0.230]);
    ROI_set_name_MW_E = uicontrol(ROI_set_name_MW,'Style','edit','String','','Units', 'normalized','fontunits','normalized', 'fontSize', 0.45,'HorizontalAlignment','left','Position',[0.10 0.44 0.800 0.190]);
    ROI_set_name_MW_OK = uicontrol(ROI_set_name_MW,'Style','pushbutton', 'String', 'OK','Units', 'normalized','fontunits','normalized', 'fontSize', 0.45,'Position',[0.10 0.16 0.310 0.180],'callback', @check_ROI_set_name);
    ROI_set_name_MW_HELP = uicontrol(ROI_set_name_MW,'Style','pushbutton', 'String', 'Help','Units', 'normalized','fontunits','normalized', 'fontSize', 0.45,'Position',[0.59 0.16 0.310 0.180], 'callback', @ROI_set_name_HW);    
    movegui(ROI_set_name_MW,'center');

    %----------------------------------------------------------------------
    % Close ROI set name GUI
    %----------------------------------------------------------------------
    function ROI_set_name_MW_EXIT(~,~)
    	delete(ROI_set_name_MW);
    	ROI_set_name = '';
    end
    
    %----------------------------------------------------------------------
    % Check ROI set name
    %----------------------------------------------------------------------
    function check_ROI_set_name(~,~)
        tmp_name = get(ROI_set_name_MW_E, 'String');
        if ~strcmp(tmp_name,'') && ~strcmp(tmp_name(1),' ')    
        	fprintf('Name of ROI set: %s.\n', tmp_name);
            delete(ROI_set_name_MW);      
            ROI_set_name = tmp_name;     
        else
            warning('Name not entered or invalid, please re-enter.');
        end
    end

    %----------------------------------------------------------------------
    % ROI set name: Help window
    %----------------------------------------------------------------------
    function ROI_set_name_HW(~,~)

        help_string = {'First, define a name for the set of ROIs. TMFC results for this ROI set will be stored in:','','"TMFC_project_name\ROI_sets\ROI_set_name"','',...
        'Second, select one or more ROI masks (*.nii files). TMFC toolbox will create a group mean binary mask based on individual subjects 1st-level masks (see SPM.VM) and apply it to all selected ROIs Empty ROIs will be excluded from further analysis. Masked ROIs will be limited to only voxels which have data for all subjects. The dimensions, orientation, and voxel sizes of the masked ROI images will be adjusted according to the group mean binary mask. These files will be stored in "Masked_ROIs"',...
        '','Third, exclude heavily cropped ROIs from further analysis, if necessary.','','Note: You can define several ROI sets and switch between them. Push the "ROI_set" button and then push "Add new ROI set". Each time you need to switch between ROI sets push the "ROI_set" button.'};

        ROI_set_name_HW_MW = figure('Name', 'Select ROIs', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.67 0.31 0.2 0.45],'Resize','off','color','w','MenuBar', 'none','ToolBar', 'none','Windowstyle', 'Modal');
        ROI_set_name_HW_MW_S = uicontrol(ROI_set_name_HW_MW,'Style','text','String', help_string,'Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.0356,'HorizontalAlignment', 'left', 'Position',[0.05 0.14 0.89 0.82],'backgroundcolor',get(ROI_set_name_HW_MW,'color'));
        ROI_set_name_HW_MW_OK = uicontrol(ROI_set_name_HW_MW,'Style','pushbutton', 'String', 'OK','Units', 'normalized','fontunits','normalized', 'fontSize', 0.45,'Position',[0.34 0.06 0.30 0.06],'callback', @ROI_set_name_HW_MW_close);
        movegui(ROI_set_name_HW_MW,'center');

        function ROI_set_name_HW_MW_close(~,~)
            close(ROI_set_name_HW_MW);
        end
    end

    uiwait();
end

%% ================ [ GUI window to remove empty ROIs ]====================
function ROI_remove_empty_GUI(empty_ROI_list)

    ROI_remove_string = {'Warning, the following ROIs do not',...
                         'contain data for at least one subject and',...
                         'will be excluded from the analysis:'};

    ROI_remove_MW = figure('Name', 'Select ROIs', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.35 0.40 0.28 0.35],'Resize','off','color','w','MenuBar', 'none','ToolBar', 'none');

    ROI_remove_MW_list = uicontrol(ROI_remove_MW , 'Style', 'listbox', 'String', empty_ROI_list,'Max', 100,'Units', 'normalized', 'Position',[0.048 0.22 0.91 0.40],'fontunits','normalized', 'fontSize', 0.105,'Value', []);
    ROI_remove_MW_S1 = uicontrol(ROI_remove_MW,'Style','text','String',ROI_remove_string,'Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.22,'backgroundcolor',get(ROI_remove_MW,'color'), 'Position',[0.20 0.73 0.600 0.2]);
    ROI_remove_MW_S2 = uicontrol(ROI_remove_MW,'Style','text','String', 'Empty ROIs:','Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.55,'backgroundcolor',get(ROI_remove_MW,'color'), 'Position',[0.04 0.62 0.200 0.08]);    
    ROI_remove_MW_OK = uicontrol(ROI_remove_MW,'Style','pushbutton', 'String', 'OK','Units', 'normalized','fontunits','normalized', 'fontSize', 0.4, 'Position',[0.38 0.07 0.28 0.10],'callback', @ROI_remove_MW_close);
    movegui(ROI_remove_MW,'center');

    %----------------------------------------------------------------------
    % Close GUI
    %----------------------------------------------------------------------
    function ROI_remove_MW_close(~,~)
        close(ROI_remove_MW);
    end

    disp(['Removed ' num2str(length(empty_ROI_list)) ' ROI(s) from the ROI set.']);
    uiwait();  
end


%% =============[ GUI window to remove heavily cropped ROIs ]==============
function [ROI_set_crop] = ROI_remove_crop(ROI_set)

    ROI_string = {};
    ROI_set_crop = []; 

    for iROI = 1:length(ROI_set.ROIs)
        full_string = {iROI,horzcat('No ',num2str(iROI),': ',ROI_set.ROIs(iROI).name, ' :: ', ...
            num2str(ROI_set.ROIs(iROI).raw_size),' voxels', ' :: ' , num2str(ROI_set.ROIs(iROI).masked_size), ...
            ' voxels ' , ':: ',num2str(ROI_set.ROIs(iROI).masked_size_percents), ' %'), ROI_set.ROIs(iROI).masked_size_percents};
        ROI_string = vertcat(ROI_string, full_string);
    end

    ROI_crop_MW_L1 = ROI_string;
    ROI_crop_MW_L2 = {};
    ROI_crop_MW_INFO1 = {'Remove heavily cropped ROIs with insufficient data, if necessary.'};
    ROI_crop_MW_INFO2 = {'No # :: ROI name :: Voxels before masking :: Voxels after masking :: Percent left'};    
    ROI_crop_MW_IND1 = {};          
    ROI_crop_MW_IND2 = {};         

    ROI_crop_MW = figure('Name', 'Select ROIs', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.35 0.40 0.32 0.55],'Resize','off','color','w','MenuBar', 'none','ToolBar', 'none','Windowstyle', 'Modal','CloseRequestFcn', @ROI_crop_MW_EXIT);
    ROI_crop_MW_LB1 = uicontrol(ROI_crop_MW , 'Style', 'listbox', 'String', ROI_crop_MW_L1(:,2,1),'Max', 100,'Units', 'normalized', 'Position',[0.048 0.565 0.91 0.30],'fontunits','normalized', 'fontSize', 0.098, 'Value', [], 'callback', @LB1_SEL);
    ROI_crop_MW_LB2 = uicontrol(ROI_crop_MW , 'Style', 'listbox', 'String', ROI_crop_MW_L2,'Max', 100,'Units', 'normalized', 'Position',[0.048 0.14 0.91 0.25],'fontunits','normalized', 'fontSize', 0.119, 'Value', [], 'callback', @LB2_SEL);

    ROI_crop_MW_S1 = uicontrol(ROI_crop_MW,'Style','text','String', ROI_crop_MW_INFO1,'Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.54,'Position',[0.10 0.92 0.8 0.05],'backgroundcolor',get(ROI_crop_MW,'color'));
    ROI_crop_MW_S2 = uicontrol(ROI_crop_MW,'Style','text','String', ROI_crop_MW_INFO2,'Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.64,'HorizontalAlignment', 'left','Position',[0.048 0.87 0.91 0.040],'backgroundcolor',get(ROI_crop_MW,'color'));
    ROI_crop_MW_S3 = uicontrol(ROI_crop_MW,'Style','text','String', '% threshold','Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.44,'HorizontalAlignment', 'left','Position',[0.84 0.475 0.13 0.055],'backgroundcolor',get(ROI_crop_MW,'color'));
    ROI_crop_MW_S4 = uicontrol(ROI_crop_MW,'Style','text','String', 'Removed ROIs:','Units', 'normalized', 'fontunits','normalized', 'fontSize', 0.50,'HorizontalAlignment', 'left','Position',[0.05 0.395 0.2 0.05],'backgroundcolor',get(ROI_crop_MW,'color'));

    ROI_crop_MW_remove_selected = uicontrol(ROI_crop_MW,'Style','pushbutton', 'String', 'Remove selected','Units', 'normalized','fontunits','normalized', 'fontSize', 0.4,'Position',[0.047 0.48 0.24 0.063], 'callback', @remove_selected);
    ROI_crop_MW_remove_thresholded = uicontrol(ROI_crop_MW,'Style','pushbutton', 'String', 'Remove ROIs under % threshold','Units', 'normalized','fontunits','normalized', 'fontSize', 0.4,'Position',[0.32 0.48 0.40 0.063], 'callback', @remove_thresholded);
    ROI_crop_MW_OK = uicontrol(ROI_crop_MW,'Style','pushbutton', 'String', 'OK','Units', 'normalized','fontunits','normalized', 'fontSize', 0.4,'Position',[0.047 0.056 0.24 0.063], 'callback', @confirm_selection);

    ROI_crop_MW_return_selected = uicontrol(ROI_crop_MW,'Style','pushbutton', 'String', 'Return selected','Units', 'normalized','fontunits','normalized', 'fontSize', 0.4,'Position',[0.39 0.056 0.24 0.063], 'callback', @return_selected);
    ROI_crop_MW_return_all = uicontrol(ROI_crop_MW,'Style','pushbutton', 'String', 'Return all','Units', 'normalized','fontunits','normalized', 'fontSize', 0.4,'Position',[0.72 0.056 0.24 0.063], 'callback', @retern_all);
    ROI_crop_MW_thr = uicontrol(ROI_crop_MW,'Style','edit','String',[],'Units', 'normalized','fontunits','normalized', 'fontSize', 0.42,'HorizontalAlignment','center','Position',[0.74 0.48 0.1 0.06]);
    movegui(ROI_crop_MW,'center');    

    %----------------------------------------------------------------------
    function ROI_crop_MW_EXIT(~,~)
    	delete(ROI_crop_MW);
    	disp('ROIs not selected.');
    end

    %----------------------------------------------------------------------
    function LB1_SEL(~,~)
    	index = get(ROI_crop_MW_LB1, 'Value');  
        ROI_crop_MW_IND1 = index;      
    end

    %----------------------------------------------------------------------
    function LB2_SEL(~,~)
        index = get(ROI_crop_MW_LB2, 'Value');  
        ROI_crop_MW_IND2 = index;             
    end

    %----------------------------------------------------------------------
    function remove_selected(~,~)
        if isempty(ROI_crop_MW_IND1)
            warning('No ROIs selected for removal.');
        else
            updated_ROIs = {};  
            new_ROIs_flag = 0;

            % Create list of selected ROIs
            updated_ROIs = vertcat(updated_ROIs, ROI_crop_MW_L1(ROI_crop_MW_IND1,:,:)); 

            if ~isempty(ROI_crop_MW_L2)
                % Condition 1: if ROIs have been previously selected, execute after checking for duplicates 
                if size(updated_ROIs,1) >= 2 
                    % Case, when more than 2 ROIs are selected
                    temp_ROI_list = []; 
                    ROI_index = 1; 
                    for iROI = 1:size(updated_ROIs,1) 
                        for jROI = 1:size(ROI_crop_MW_L2,1)
                            if strcmp(updated_ROIs(iROI,2,1), ROI_crop_MW_L2(jROI,2,1))
                               temp_ROI_list(ROI_index) = iROI;
                               ROI_index = ROI_index+1;
                            end
                        end
                    end            
                else
                    % Case, when only 1 ROI is selected
                    temp_ROI_list = [];
                    ROI_index = 1;
                    for iROI = 1:size(ROI_crop_MW_L2,1)
                        if strcmp(updated_ROIs(1,2,1),ROI_crop_MW_L2(iROI,2,1))
                            temp_ROI_list(ROI_index) = iROI;
                            ROI_index = ROI_index+1;
                        end
                    end                 
                end

                % Remove the respective ROIs in the main ROI list for GUI
                if length(temp_ROI_list)>=2
                    % Case, when more than 2 ROIs are selected
                    set_index = 0;                   
                    for iROI = 1:length(temp_ROI_list)
                        updated_ROIs(temp_ROI_list(iROI)-set_index,:,:) = [];
                        set_index = set_index + 1;
                    end
                    new_ROIs_flag = size(updated_ROIs,1);  
                else
                    % Case, when only 1 ROI is selected
                    for iROI = 1:size(updated_ROIs,1)
                        if updated_ROIs{iROI,1,1} == temp_ROI_list
                            updated_ROIs(iROI,:,:) = [];
                        end
                    end
                    new_ROIs_flag = size(updated_ROIs,1);
                end
                ROI_crop_MW_L2 = sortrows(vertcat(ROI_crop_MW_L2, updated_ROIs),1);              
            else
                % Condition 2: if ROIs have not been previously selected, directly add to list
                ROI_crop_MW_L2 = vertcat(ROI_crop_MW_L2, updated_ROIs); 
                new_ROIs_flag = 2;
            end

            % Check if newly selected ROIs for removal have been added
            if new_ROIs_flag == 2
                    fprintf('ROIs selected for removal: %d. \n', size(ROI_crop_MW_L2,1));
            elseif new_ROIs_flag == 0
                    warning('Selected ROIs are already present in the removal list, no new ROIs to remove.');
            else
                    fprintf('New selected ROIs for removal: %d. \n', new_ROIs_flag); 
            end 

            % Update sorted list of ROIs into GUI
            set(ROI_crop_MW_LB2, 'String', ROI_crop_MW_L2(:,2,1));
        end
    end

    %----------------------------------------------------------------------
    function remove_thresholded(~,~)

        updated_ROIs = {}; 
        new_ROIs_flag = 0;  

        ROI_crop_thr = get(ROI_crop_MW_thr, 'String'); % Get threshold

        if ~strcmp(ROI_crop_thr,'')

            threshold = str2double(ROI_crop_thr); 

            if isnan(threshold)
                warning('Entered threshold should be a natural number, please re-enter.');

            elseif (threshold<0) || (threshold>100)
                warning('Please enter a threshold between 0 and 100%.');

            else
                temp_ROI_list = [];
                ROI_index = 1;
                for iROI = 1:size(ROI_crop_MW_L1,1)
                    if ROI_crop_MW_L1{iROI,3,1} <= threshold
                        temp_ROI_list(ROI_index) = iROI;
                        ROI_index = ROI_index + 1;
                    end
                end
                updated_ROIs = ROI_crop_MW_L1(temp_ROI_list,:,:);

                % Compiling the list of ROIs to remove
                if isempty(ROI_crop_MW_L2)
                    % If List is exported for the firt time
                    ROI_crop_MW_L2 = vertcat(ROI_crop_MW_L2, updated_ROIs); 
                    new_ROIs_flag = 2;
                else
                    % Check for duplicates
                    if size(updated_ROIs,1) >= 2
                        % Case, when more than 2 ROIs are present under threshold
                        temp_ROI_list = []; 
                        ROI_index = 1;      
                        for iROI = 1:size(updated_ROIs,1)
                            for jSet = 1:size(ROI_crop_MW_L2,1)
                                if strcmp(updated_ROIs(iROI,2,1), ROI_crop_MW_L2(jSet,2,1))
                                   temp_ROI_list(ROI_index) = iROI;
                                   ROI_index = ROI_index+1;
                                end
                            end
                        end
                    else
                        % Case, when only 1 ROI is present under threshold
                        temp_ROI_list = []; 
                        ROI_index = 1;
                        for iROI = 1:size(ROI_crop_MW_L2,1)
                            if strcmp(updated_ROIs(1,2,1),ROI_crop_MW_L2(iROI,2,1))
                                temp_ROI_list(ROI_index) = iROI;
                                ROI_index = ROI_index+1;
                            end
                        end
                    end

                    % Remove the respective ROIs in the main ROI list for GUI
                    if length(temp_ROI_list)>=2
                        % Case, when more than 2 ROIs are selected
                        ROI_index = 0;      
                        for c = 1:length(temp_ROI_list)
                            updated_ROIs(temp_ROI_list(c)-ROI_index,:,:) = [];
                            ROI_index = ROI_index + 1;
                        end
                        new_ROIs_flag = size(updated_ROIs,1);
                    else
                        % Case, when only 1 ROI is selected
                        updated_ROIs(temp_ROI_list,:,:) = [];
                        new_ROIs_flag = size(updated_ROIs,1);
                    end
                    ROI_crop_MW_L2 = sortrows(vertcat(ROI_crop_MW_L2, updated_ROIs),1);
                end

                % Check if newly selected ROIs have been added
                if new_ROIs_flag == 2 && size(ROI_crop_MW_L2,1) ~= 0
                        fprintf('ROIs selected for removal: %d. \n', size(ROI_crop_MW_L2,1));

                elseif new_ROIs_flag(1) == 2 && size(ROI_crop_MW_L2,1) == 0
                        warning('ROIs below this threshold do not exist.');

                    elseif new_ROIs_flag(1) == 0
                        warning('All ROIs below this threshold have already been removed.');

                    else
                        fprintf('%d ROIs selected for removal at threshold %d percents. \n', new_ROIs_flag(1),threshold);     
                end    

                % Update sorted list of ROIs into GUI
                set(ROI_crop_MW_LB2, 'String', ROI_crop_MW_L2(:,2,1));

            end

        else            
            warning('The entered threshold is empty or invalid, please re-enter.');
        end
    end

    %----------------------------------------------------------------------
    function return_selected(~,~)
        if isempty(ROI_crop_MW_L2)
            warning('No ROIs present to return.');
        elseif isempty(ROI_crop_MW_IND2)
            warning('No ROIs selected to return.');
        else
            if length(ROI_crop_MW_IND2) >= 2
                set_index = 0;
                for c = 1:length(ROI_crop_MW_IND2)
                    ROI_crop_MW_L2(ROI_crop_MW_IND2(c)-set_index,:,:) = [];
                    set_index = set_index + 1;
                end
                fprintf('Number of ROIs removed are %d. \n', set_index);
            else
                ROI_crop_MW_L2(ROI_crop_MW_IND2,:,:) = [];
                disp('Selected ROI has been returned.');
            end

            if isempty(ROI_crop_MW_L2)
                ROI_crop_MW_L2 = {};
               set(ROI_crop_MW_LB2, 'String', ROI_crop_MW_L2); 
            else
                set(ROI_crop_MW_LB2, 'String', ROI_crop_MW_L2(:,2,1));
                set(ROI_crop_MW_LB2, 'Value', []);
            end
        end
    end

    %----------------------------------------------------------------------
    function retern_all(~,~)
        if isempty(ROI_crop_MW_L2)
            warning('No ROIs present to return.');
        else
            ROI_crop_MW_L2 = {};
            set(ROI_crop_MW_LB2, 'String', ROI_crop_MW_L2);
            set(ROI_crop_MW_LB2, 'Value', []);
            fprintf('%d ROIs have been returned. \n',size(ROI_crop_MW_L2,1));
        end
    end
    
    %----------------------------------------------------------------------
    function confirm_selection(~,~)
        if isempty(ROI_crop_MW_L2)
            disp('New ROI set has been defined. All selected ROIs have been saved.');
            delete(ROI_crop_MW);
        else
            if length(ROI_crop_MW_L1) == length(ROI_crop_MW_L2)
                warning('All ROIs have been removed, please try again.');
            else
                disp('New ROI set has been defined. Highly cropped ROIs have been removed.');
                ROI_index = 0;
                for jROI = 1:size(ROI_crop_MW_L2,1)
                    ROI_set.ROIs(ROI_crop_MW_L2{jROI,1,1} - ROI_index) = [];
                    ROI_index = ROI_index + 1;
                end
                delete(ROI_crop_MW);
            end
        end
        ROI_set_crop = ROI_set;    
    end

    uiwait();
end



