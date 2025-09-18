function [sub_check] = tmfc_ROI_to_ROI_contrast(tmfc,type,contrast_number,ROI_set_number)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Calculates linear contrasts of FC matrices.
%
% FORMAT [sub_check] = tmfc_ROI_to_ROI_contrast(tmfc,type,contrast_number)
%
%   type                   - TMFC analysis type
%                            1: gPPI
%                            2: gPPI-FIR
%                            3: BSC-LSS
%                            4: BSC-LSS after FIR
%
%   contrast_number        - Indices of contrasts to compute in tmfc struct
%    
%   tmfc.subjects.path            - Paths to individual SPM.mat files
%   tmfc.subjects.name            - Subject names within the TMFC project
%                           ('Subject_XXXX' naming will be used by default)
%   tmfc.project_path             - Path where all results will be saved
%   tmfc.defaults.parallel        - 0 or 1 (sequential/parallel computing)
%   tmfc.ROI_set                  - List of selected ROIs
%   tmfc.ROI_set.type             - Type of the ROI set
%   tmfc.ROI_set.set_name         - Name of the ROI set
%   tmfc.ROI_set.ROIs.name        - Name of the selected ROI
%   tmfc.ROI_set.ROIs.path_masked - Paths to the ROI images masked by group
%                                   mean binary mask 
%   
% -------------------------------------------------------------------------
% For gPPI and gPPI-FIR analyses:
% tmfc.ROI_set(ROI_set_number).gPPI.conditions - List of conditions of
%                                                interest for gPPI and
%                                                gPPI-FIR analyses
%                                                (see tmfc_PPI)
% For BSC analyses:
% tmfc.LSS.conditions           - List of conditions of interest
%                                 for BSC-LSS analysis
%                                 (see tmfc_LSS)
% tmfc.LSS_after_FIR.conditions - List of conditions of interest
%                                 for BSC-LSS (after FIR) analysis
%                                 (see tmfc_LSS_after_FIR)
% -------------------------------------------------------------------------
%
%
% Example of an ROI set (see tmfc_select_ROIs_GUI):
%
%   tmfc.ROI_set(1).set_name = 'two_ROIs';
%   tmfc.ROI_set(1).type = 'binary_images';
%   tmfc.ROI_set(1).ROIs(1).name = 'ROI_1';
%   tmfc.ROI_set(1).ROIs(2).name = 'ROI_2';
%   tmfc.ROI_set(1).ROIs(1).path_masked = 'C:\ROI_set\two_ROIs\ROI_1.nii';
%   tmfc.ROI_set(1).ROIs(2).path_masked = 'C:\ROI_set\two_ROIs\ROI_2.nii';                             
%
% FORMAT [sub_check] = tmfc_ROI_to_ROI_contrast(tmfc,type,contrast_number,ROI_set_number)
% Run the function for the selected ROI set
%
%   tmfc                   - As above
%   ROI_set_number         - Number of the ROI set in the tmfc structure
%   
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru


if nargin < 4
   ROI_set_number = 1;
end

% Check subject names
if ~isfield(tmfc.subjects,'name')
    for iSub = 1:length(tmfc.subjects)
        tmfc.subjects(iSub).name = ['Subject_' num2str(iSub,'%04.f')];
    end
end

w = waitbar(0,'Please wait...','Name','Compute ROI-to-ROI contrasts');
start_time = tic;
count_sub = 1;
cleanupObj = onCleanup(@unfreeze_after_ctrl_c);

nSub = length(tmfc.subjects);
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
nCon = length(contrast_number);

switch type
    %================================gPPI==================================
    case 1
        for iSub = 1:nSub
            % Load default contrasts for conditions of interest
            cond_list = tmfc.ROI_set(ROI_set_number).gPPI.conditions;
            for jCond = 1:length(cond_list)                             
                % If ROI-to-ROI analysis has not been performed for the default contrasts
                if ~exist(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI','ROI_to_ROI','symmetrical', ...
                         [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']),'file')
                    disp('ROI-to-ROI analysis has not been performed previously. Performing ROI-to-ROI gPPI analysis, please wait...');
                    analysis = tmfc.defaults.analysis;
                    tmfc.defaults.analysis = 2;
                    tmfc_gPPI(tmfc,ROI_set_number,1);
                    tmfc.defaults.analysis = analysis;
                end
                
                load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI','ROI_to_ROI','asymmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']));
                load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI','ROI_to_ROI','symmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']));
    
                matrices(jCond,:) = ppi_matrix(:)';
                symm_matrices(jCond,:) = symm_ppi_matrix(:)';
                clear ppi_matrix symm_ppi_matrix
            end
            % Calculate and save contrasts
            for jCon = 1:nCon
                ppi_matrix = reshape(tmfc.ROI_set(ROI_set_number).contrasts.gPPI(contrast_number(jCon)).weights*matrices,[nROI,nROI]);
                symm_ppi_matrix = reshape(tmfc.ROI_set(ROI_set_number).contrasts.gPPI(contrast_number(jCon)).weights*symm_matrices,[nROI,nROI]);
                save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI','ROI_to_ROI','asymmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(contrast_number(jCon),'%04.f') ...
                    '_[' regexprep(tmfc.ROI_set(ROI_set_number).contrasts.gPPI(contrast_number(jCon)).title,' ','_') '].mat']),'ppi_matrix');
                save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI','ROI_to_ROI','symmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(contrast_number(jCon),'%04.f') ...
                    '_[' regexprep(tmfc.ROI_set(ROI_set_number).contrasts.gPPI(contrast_number(jCon)).title,' ','_') '].mat']),'symm_ppi_matrix');
                clear ppi_matrix symm_ppi_matrix
            end
            % Update waitbar
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            count_sub = count_sub + 1;
            time_remaining = (nSub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            try
                waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end      
            sub_check(iSub) = 1;
            clear matrices
        end 

    %==============================gPPI-FIR================================
    case 2
        for iSub = 1:nSub          
            % Load default contrasts for conditions of interest
            cond_list = tmfc.ROI_set(ROI_set_number).gPPI.conditions;
            for jCond = 1:length(cond_list)                              
                % If ROI-to-ROI analysis has not been performed for the default contrasts
                if ~exist(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','symmetrical', ...
                                  [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']),'file')
                    disp('ROI-to-ROI analysis has not been performed previously. Performing ROI-to-ROI gPPI-FIR analysis, please wait...');
                    analysis = tmfc.defaults.analysis;
                    tmfc.defaults.analysis = 2;
                    tmfc_gPPI_FIR(tmfc,ROI_set_number,1);
                    tmfc.defaults.analysis = analysis;
                end
    
                load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','asymmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']));
                load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','symmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']));
    
                matrices(jCond,:) = ppi_matrix(:)';
                symm_matrices(jCond,:) = symm_ppi_matrix(:)';
                clear ppi_matrix symm_ppi_matrix
            end
            % Calculate and save contrasts
            for jCon = 1:nCon
                ppi_matrix = reshape(tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(contrast_number(jCon)).weights*matrices,[nROI,nROI]);
                symm_ppi_matrix = reshape(tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(contrast_number(jCon)).weights*symm_matrices,[nROI,nROI]);
                save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','asymmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(contrast_number(jCon),'%04.f') ...
                    '_[' regexprep(tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(contrast_number(jCon)).title,' ','_') '].mat']),'ppi_matrix');
                save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','symmetrical', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(contrast_number(jCon),'%04.f') ...
                    '_[' regexprep(tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(contrast_number(jCon)).title,' ','_') '].mat']),'symm_ppi_matrix');
                clear ppi_matrix symm_ppi_matrix
            end
            % Update waitbar
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            count_sub = count_sub + 1;
            time_remaining = (nSub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            try
                waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end 
            sub_check(iSub) = 1;
            clear matrices
        end 

    %===============================BSC-LSS================================
    case 3
        for iSub = 1:nSub
            % Load default contrasts for conditions of interest
            cond_list = tmfc.LSS.conditions;
            for jCond = 1:length(cond_list)                            
                % If ROI-to-ROI analysis has not been performed for the default contrasts
                if ~exist(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS','ROI_to_ROI', ...
                         [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']),'file')
                    disp('ROI-to-ROI analysis has not been performed previously. Performing ROI-to-ROI BSC-LSS analysis, please wait...');
                    analysis = tmfc.defaults.analysis;
                    tmfc.defaults.analysis = 2;
                    tmfc_BSC(tmfc,ROI_set_number,0);
                    tmfc.defaults.analysis = analysis;
                end
                
                load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS','ROI_to_ROI', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']));
    
                matrices(jCond,:) = z_matrix(:)';
                clear z_matrix
            end
            % Calculate and save contrasts
            for jCon = 1:nCon
                z_matrix = reshape(tmfc.ROI_set(ROI_set_number).contrasts.BSC(contrast_number(jCon)).weights*matrices,[nROI,nROI]);
                save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS','ROI_to_ROI', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(contrast_number(jCon),'%04.f') ...
                    '_[' regexprep(tmfc.ROI_set(ROI_set_number).contrasts.BSC(contrast_number(jCon)).title,' ','_') '].mat']),'z_matrix');
                clear z_matrix
            end
            % Update waitbar
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            count_sub = count_sub + 1;
            time_remaining = (nSub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            try
                waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end         
            sub_check(iSub) = 1;
            clear matrices
        end 
    
    %==========================BSC-LSS after FIR===========================
    case 4
        for iSub = 1:nSub
            % Load default contrasts for conditions of interest
            cond_list = tmfc.LSS_after_FIR.conditions;
            for jCond = 1:length(cond_list)               
                % If ROI-to-ROI analysis has not been performed for the default contrasts
                if ~exist(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','ROI_to_ROI', ...
                         [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']),'file')
                    disp('ROI-to-ROI analysis has not been performed previously. Performing ROI-to-ROI BSC-LSS (after FIR) analysis, please wait...');
                    analysis = tmfc.defaults.analysis;
                    tmfc.defaults.analysis = 2;
                    tmfc_BSC_after_FIR(tmfc,ROI_set_number,0);
                    tmfc.defaults.analysis = analysis;
                end
                
                load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','ROI_to_ROI', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']));

                matrices(jCond,:) = z_matrix(:)';
                clear z_matrix
            end
            % Calculate and save contrasts
            for jCon = 1:nCon
                z_matrix = reshape(tmfc.ROI_set(ROI_set_number).contrasts.BSC_after_FIR(contrast_number(jCon)).weights*matrices,[nROI,nROI]);
                save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','ROI_to_ROI', ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(contrast_number(jCon),'%04.f') ...
                    '_[' regexprep(tmfc.ROI_set(ROI_set_number).contrasts.BSC_after_FIR(contrast_number(jCon)).title,' ','_') '].mat']),'z_matrix');
                clear z_matrix
            end
            % Update waitbar
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            count_sub = count_sub + 1;
            time_remaining = (nSub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            try
                waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end       
            sub_check(iSub) = 1;
            clear matrices
        end    
end

% Close waitbar
try
    delete(w)
end

function unfreeze_after_ctrl_c()    
    try
        delete(findall(0,'type','figure','Tag', 'tmfc_waitbar'));
        GUI = guidata(findobj('Tag','TMFC_GUI')); 
        set([GUI.TMFC_GUI_B1, GUI.TMFC_GUI_B2, GUI.TMFC_GUI_B3, GUI.TMFC_GUI_B4,...
           GUI.TMFC_GUI_B5a, GUI.TMFC_GUI_B5b, GUI.TMFC_GUI_B6, GUI.TMFC_GUI_B7,...
           GUI.TMFC_GUI_B8, GUI.TMFC_GUI_B9, GUI.TMFC_GUI_B10, GUI.TMFC_GUI_B11,...
           GUI.TMFC_GUI_B12a,GUI.TMFC_GUI_B12b,GUI.TMFC_GUI_B13a,GUI.TMFC_GUI_B13b,...
           GUI.TMFC_GUI_B14a, GUI.TMFC_GUI_B14b], 'Enable', 'on');
    end 
end
end
