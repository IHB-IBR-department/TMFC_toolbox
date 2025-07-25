function [sub_check] = tmfc_BGFC(tmfc,ROI_set_number,start_sub)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Calculates background functional connectivity (BGFC).
% 
% Extracts residual time-series from volumes of interest (VOIs). Regresses
% out confounds and co-activations using FIR model. Applies whitening and
% high-pass filtering. Calculates Pearson's correlation between residual
% time-series. Converts Pearson's r to Fisher's Z.
%
% FORMAT [sub_check] = tmfc_BGFC(tmfc)
% Run a function starting from the first subject in the list.
%
%   tmfc.subjects.path            - Paths to individual SPM.mat files
%   tmfc.subjects.name            - Subject names within the TMFC project
%                           ('Subject_XXXX' naming will be used by default)
%   tmfc.project_path             - Path where all results will be saved
%   tmfc.defaults.parallel        - 0 or 1 (sequential/parallel computing)
%
%   tmfc.ROI_set                  - List of selected ROIs
%   tmfc.ROI_set.type             - Type of the ROI set
%   tmfc.ROI_set.set_name         - Name of the ROI set
%   tmfc.ROI_set.ROIs.name        - Name of the selected ROI
%   tmfc.ROI_set.ROIs.path_masked - Paths to the ROI images masked by group
%                                   mean binary mask 
%
% Example of the ROI set (see tmfc_select_ROIs_GUI):
%
%   tmfc.ROI_set(1).set_name = 'two_ROIs';
%   tmfc.ROI_set(1).type = 'binary_images';
%   tmfc.ROI_set(1).ROIs(1).name = 'ROI_1';
%   tmfc.ROI_set(1).ROIs(2).name = 'ROI_2';
%   tmfc.ROI_set(1).ROIs(1).path_masked = 'C:\ROI_set\two_ROIs\ROI_1.nii';
%   tmfc.ROI_set(1).ROIs(2).path_masked = 'C:\ROI_set\two_ROIs\ROI_2.nii';
%
% FORMAT [sub_check] = tmfc_BGFC(tmfc,ROI_set_number,start_sub)
% Run the function starting from a specific subject in the path list for
% the selected ROI set.
%
%   tmfc                   - As above
%   ROI_set_number         - Number of the ROI set in the tmfc structure
%   start_sub              - Subject number on the path list to start with
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

if nargin == 1
   ROI_set_number = 1;
   start_sub = 1;
elseif nargin == 2
   start_sub = 1;
end

% Check subject names
if ~isfield(tmfc.subjects,'name')
    for iSub = 1:length(tmfc.subjects)
        tmfc.subjects(iSub).name = ['Subject_' num2str(iSub,'%04.f')];
    end
end

nSub = length(tmfc.subjects);
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
SPM = load(tmfc.subjects(1).path).SPM;
nSess = length(SPM.Sess);
sub_check = zeros(1,nSub);
if start_sub > 1
    sub_check(1:start_sub) = 1;
end

% Initialize waitbar
w = waitbar(0,'Please wait...','Name','Calculating residuals and BGFC','Tag','tmfc_waitbar');
start_time = tic;
count_sub = 1;
cleanupObj = onCleanup(@unfreeze_after_ctrl_c);

spm('defaults','fmri');
spm_jobman('initcfg');

if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','ROI_to_ROI'))
    mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','ROI_to_ROI'));
end

for iSub = start_sub:nSub
    if isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',tmfc.subjects(iSub).name))
        rmdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',tmfc.subjects(iSub).name),'s');
    end

    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',tmfc.subjects(iSub).name))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',tmfc.subjects(iSub).name));
    end

    for jSess = 1:nSess
        for kROI = 1:nROI
            matlabbatch{1}.spm.util.voi.spmmat = {fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name,'SPM.mat')};
            matlabbatch{1}.spm.util.voi.adjust = NaN; % Adjust for everything 
            matlabbatch{1}.spm.util.voi.session = jSess;
            matlabbatch{1}.spm.util.voi.name = tmfc.ROI_set(ROI_set_number).ROIs(kROI).name;
            switch tmfc.ROI_set(ROI_set_number).type
                case {'binary_images','fixed_spheres'}
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {tmfc.ROI_set(ROI_set_number).ROIs(kROI).path_masked};
                otherwise
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {tmfc.ROI_set(ROI_set_number).ROIs(kROI).path_masked(iSub).subjects};
            end
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.1;
            matlabbatch{1}.spm.util.voi.expression = 'i1';           
            batch{kROI} = matlabbatch;
            clear matlabbatch
        end
        
        switch tmfc.defaults.parallel
            case 0  % Sequential
                for kROI = 1:nROI
                    spm('defaults','fmri');
                    spm_jobman('initcfg');
                    spm_get_defaults('cmdline',true);
                    spm_jobman('run',batch{kROI});
                    movefile(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs', ... 
                                      tmfc.subjects(iSub).name,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']));
                    if exist(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']),'file')
                        delete(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']));
                    else
                        delete(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
                end
                
            case 1  % Parallel
                try
                    parpool;
                    figure(findobj('Tag','TMFC_GUI'));
                end
                
                parfor kROI = 1:nROI
                    spm('defaults','fmri');
                    spm_jobman('initcfg');
                    spm_get_defaults('cmdline',true);
                    spm_jobman('run',batch{kROI});
                    movefile(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs', ... 
                                      tmfc.subjects(iSub).name,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']));
                    if exist(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']),'file')
                        delete(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']));
                    else
                        delete(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(tmfc.project_path,'FIR_regression',tmfc.subjects(iSub).name, ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
                end
        end

        clear batch
    end

    % Calculate BGFC matrix
    for jSess = 1:nSess
        for kROI = 1:nROI
            Y(:,kROI) = load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs', ... 
                tmfc.subjects(iSub).name,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat'])).Y;
        end
        z_matrix = atanh(corr(Y));
        z_matrix(1:size(z_matrix,1)+1:end) = nan;
        save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','ROI_to_ROI', ... 
                [tmfc.subjects(iSub).name '_Session_' num2str(jSess) '.mat']),'z_matrix');
        clear Y z_matrix  
    end
    
    sub_check(iSub) = 1;
    
    % Update waitbar
    elapsed_time = toc(start_time);
    time_per_sub = elapsed_time/count_sub;
    count_sub = count_sub + 1;
    time_remaining = (nSub-iSub)*time_per_sub;
    hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
    try
        waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
    end
end

% Close waitbar
try
    close(w)
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
