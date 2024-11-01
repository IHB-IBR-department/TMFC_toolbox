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
%   tmfc.project_path             - Path where all results will be saved
%   tmfc.defaults.parallel        - 0 or 1 (sequential/parallel computing)
%
%   tmfc.ROI_set                  - List of selected ROIs
%   tmfc.ROI_set.set_name         - Name of the ROI set
%   tmfc.ROI_set.ROIs.name        - Name of the selected ROI
%   tmfc.ROI_set.ROIs.path_masked - Paths to the ROI images masked by group
%                                   mean binary mask 
%
% Example of the ROI set:
%
%   tmfc.ROI_set(1).set_name = 'two_ROIs';
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

if nargin == 1
   ROI_set_number = 1;
   start_sub = 1;
elseif nargin == 2
   start_sub = 1;
end

nSub = length(tmfc.subjects);
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
SPM = load(tmfc.subjects(1).path);
nSess = length(SPM.SPM.Sess);

% Initialize waitbar for parallel or sequential computing
switch tmfc.defaults.parallel
    case 0    % Sequential
        w = waitbar(0,'Please wait...','Name','Calculating residuals and BGFC','Tag', 'tmfc_waitbar');
        cleanupObj = onCleanup(@unfreeze_after_ctrl_c);
    case 1    % Parallel
        try   % Waitbar for MATLAB R2017a and higher
        	D = parallel.pool.DataQueue;            
            w = waitbar(0,'Please wait...','Name','Calculating residuals and BGFC','Tag', 'tmfc_waitbar');
            afterEach(D, @tmfc_parfor_waitbar);    
            tmfc_parfor_waitbar(w,nSub,start_sub);     
        catch % No waitbar for MATLAB R2016b and earlier
            opts = struct('WindowStyle','non-modal','Interpreter','tex');
            w = warndlg({'\fontsize{12}Sorry, waitbar progress update is not available for parallel computations in MATLAB R2016b and earlier.',[],...
                'Please wait until all computations are completed.',[],...
                'If you want to interrupt computations:',...
                '   1) Do not close this window;',...
                '   2) Select MATLAB main window;',...
                '   3) Press Ctrl+C.'},'Please wait...',opts);
        end
        cleanupObj = onCleanup(@unfreeze_after_ctrl_c);
end

spm('defaults','fmri');
spm_jobman('initcfg');

if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','ROI_to_ROI'))
    mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','ROI_to_ROI'));
end

for iSub = start_sub:nSub
    tic
    if isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',['Subject_' num2str(iSub,'%04.f')]))
        rmdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',['Subject_' num2str(iSub,'%04.f')]),'s');
    end

    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',['Subject_' num2str(iSub,'%04.f')]))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs',['Subject_' num2str(iSub,'%04.f')]));
    end

    for jSess = 1:nSess
        for kROI = 1:nROI
            matlabbatch{1}.spm.util.voi.spmmat = {fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')],'SPM.mat')};
            matlabbatch{1}.spm.util.voi.adjust = NaN; % Adjust for everything 
            matlabbatch{1}.spm.util.voi.session = jSess;
            matlabbatch{1}.spm.util.voi.name = tmfc.ROI_set(ROI_set_number).ROIs(kROI).name;
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {tmfc.ROI_set(ROI_set_number).ROIs(kROI).path_masked};
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
                    movefile(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs', ... 
                                      ['Subject_' num2str(iSub,'%04.f')],['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']));
                    if exist(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']),'file')
                        delete(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']));
                    else
                        delete(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
                end
                
            case 1  % Parallel
                parfor kROI = 1:nROI
                    spm('defaults','fmri');
                    spm_jobman('initcfg');
                    spm_get_defaults('cmdline',true);
                    spm_jobman('run',batch{kROI});
                    movefile(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs', ... 
                                      ['Subject_' num2str(iSub,'%04.f')],['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']));
                    if exist(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']),'file')
                        delete(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '_eigen.nii']));
                    else
                        delete(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(tmfc.project_path,'FIR_regression',['Subject_' num2str(iSub,'%04.f')], ...
                                      ['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
                end
        end

        clear batch
    end

    % Calculate BGFC matrix
    for jSess = 1:nSess
        for kROI = 1:nROI
            VOI = load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','FIR_VOIs', ... 
                ['Subject_' num2str(iSub,'%04.f')],['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(jSess) '.mat']));
            Y(:,kROI) = VOI.Y; 
            clear VOI
        end
        z_matrix = atanh(corr(Y));
        z_matrix(1:size(z_matrix,1)+1:end) = nan;
        save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BGFC','ROI_to_ROI', ... 
                ['Subject_' num2str(iSub,'%04.f') '_Session_' num2str(jSess) '.mat']),'z_matrix');
        clear VOI z_matrix  
    end
    
    sub_check(iSub) = 1;
    pause(0.0001)
    
    % Update waitbar
    switch tmfc.defaults.parallel
        case 0  % Sequential
            hms = fix(mod(((nSub-iSub)*toc/iSub), [0, 3600, 60]) ./ [3600, 60, 1]);
            try
                waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1)) ':' num2str(hms(2)) ':' num2str(hms(3)) ' [hr:min:sec] remaining']);
            end
        case 1  % Parallel
            try 
                send(D,[]);
            end
    end
end

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
%% ========================================================================

% Waitbar for parallel mode
function tmfc_parfor_waitbar(waitbarHandle,iterations,start_sub)
    persistent count h nSub start
    if nargin == 3
        count = start_sub - 1;
        h = waitbarHandle;
        nSub = iterations;
        start = tic;        
    else
        if isvalid(h)         
            count = count + 1;
            time = toc(start);
            hms = fix(mod(((nSub-count)*time/count), [0, 3600, 60]) ./ [3600, 60, 1]);
            waitbar(count/nSub, h, [num2str(count/nSub*100,'%.f') '%, ' num2str(hms(1)) ':' num2str(hms(2)) ':' num2str(hms(3)) ' [hr:min:sec] remaining']);
        end
    end
end
