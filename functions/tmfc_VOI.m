function [sub_check] = tmfc_VOI(tmfc,ROI_set_number,start_sub)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Extracts time-series from volumes of interest (VOIs). Calculates 
% F-contrast for all conditions of interest. Regresses out conditions of no
% interest and confounds. Applies whitening and high-pass filtering.
%
% FORMAT [sub_check] = tmfc_VOI(tmfc)
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
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions        - List of conditions of interest
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.sess   - Session number (as specified in SPM.Sess)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.number - Condition number (as specified in SPM.Sess.U)
%
% Session number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two sessions. Both sessions 
% contains three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified:
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).sess   = 1;   
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).number = 1; - "Cond A", 1st session
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).sess   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).number = 2; - "Cond B", 1st session
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).number = 1; - "Cond A", 2nd session
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).number = 2; - "Cond B", 2nd session 
%
% Example of the ROI set:
%
%   tmfc.ROI_set(1).set_name = 'two_ROIs';
%   tmfc.ROI_set(1).ROIs(1).name = 'ROI_1';
%   tmfc.ROI_set(1).ROIs(2).name = 'ROI_2';
%   tmfc.ROI_set(1).ROIs(1).path_masked = 'C:\ROI_set\two_ROIs\ROI_1.nii';
%   tmfc.ROI_set(1).ROIs(2).path_masked = 'C:\ROI_set\two_ROIs\ROI_2.nii';
%
% FORMAT [sub_check] = tmfc_VOI(tmfc,ROI_set_number,start_sub)
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

try
    main_GUI = guidata(findobj('Tag','TMFC_GUI'));                           
    set(main_GUI.TMFC_GUI_S3,'String', 'Updating...','ForegroundColor',[0.772, 0.353, 0.067]);       
end

nSub = length(tmfc.subjects);
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
cond_list = tmfc.ROI_set(ROI_set_number).gPPI.conditions;
nCond = length(cond_list);
sess = []; sess_num = []; nSess = [];
for iSub = 1:nCond
    sess(iSub) = cond_list(iSub).sess;
end
sess_num = unique(sess);
nSess = length(sess_num);
sub_check = zeros(1,nSub);
if start_sub > 1
    sub_check(1:start_sub) = 1;
end

% Initialize waitbar for sequential or parallel computations
switch tmfc.defaults.parallel
    case 0    % Sequential
        w = waitbar(0,'Please wait...','Name','VOI time-series extraction','Tag','tmfc_waitbar');
        cleanupObj = onCleanup(@unfreeze_after_ctrl_c);
    case 1    % Parallel
        try   % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;             
            w = waitbar(0,'Please wait...','Name','VOI time-series extraction','Tag','tmfc_waitbar');
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
        
        try   % Bring TMFC main window to the front 
            figure(findobj('Tag','TMFC_GUI'));
        end
end

spm('defaults','fmri');
spm_jobman('initcfg');

for iSub = start_sub:nSub
    tic
    % Calculate F-contrast for all conditions of interest
    SPM = load(tmfc.subjects(iSub).path);
    matlabbatch{1}.spm.stats.con.spmmat = {tmfc.subjects(iSub).path};
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'F_conditions_of_interest';
    weights = zeros(nCond,size(SPM.SPM.xX.X,2));
    for jCond = 1:nCond
        weights(jCond,SPM.SPM.Sess(cond_list(jCond).sess).col(cond_list(jCond).number)) = 1;
    end
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = weights;
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    spm_get_defaults('cmdline',true);
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    SPM = load(tmfc.subjects(iSub).path);

    if isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',['Subject_' num2str(iSub,'%04.f')]))
        rmdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',['Subject_' num2str(iSub,'%04.f')]),'s');
    end

    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',['Subject_' num2str(iSub,'%04.f')]))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',['Subject_' num2str(iSub,'%04.f')]));
    end
    
    for jSess = 1:nSess
        for kROI = 1:nROI
            matlabbatch{1}.spm.util.voi.spmmat = {tmfc.subjects(iSub).path};
            matlabbatch{1}.spm.util.voi.adjust = length(SPM.SPM.xCon);
            matlabbatch{1}.spm.util.voi.session = sess_num(jSess);
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
                    movefile(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs', ... 
                                      ['Subject_' num2str(iSub,'%04.f')],['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']));
                    if exist(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']),'file')
                        delete(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']));
                    else
                        delete(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
                end
                
            case 1  % Parallel
                parfor kROI = 1:nROI
                    spm('defaults','fmri');
                    spm_jobman('initcfg');
                    spm_get_defaults('cmdline',true);
                    spm_jobman('run',batch{kROI});
                    movefile(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs', ... 
                                      ['Subject_' num2str(iSub,'%04.f')],['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']));
                    if exist(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']),'file')
                        delete(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']));
                    else
                        delete(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(SPM.SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
                end
        end

        clear batch
    end
    
    sub_check(iSub) = 1;
    pause(0.0001)

    % Update main TMFC GUI
    try  
        main_GUI = guidata(findobj('Tag','TMFC_GUI'));                        
        set(main_GUI.TMFC_GUI_S3,'String', strcat(num2str(iSub), '/', num2str(nSub), ' done'),'ForegroundColor',[0.219, 0.341, 0.137]);    
    end
    
    % Update waitbar
    switch tmfc.defaults.parallel
        case 0   % Sequential
            hms = fix(mod(((nSub-iSub)*toc), [0, 3600, 60]) ./ [3600, 60, 1]);
            try
                waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end   
        case 1   % Parallel
            try
                send(D,[]);
            end
    end

    clear SPM
end

% Close waitbar
try
    delete(w);
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
function tmfc_parfor_waitbar(waitbarHandle,iterations,firstsub)
    persistent count h N start t1 t2
    if nargin == 3
        count = firstsub - 1;
        h = waitbarHandle;
        N = iterations;
        start = tic;
        t1 = 0; t2 = 0;
    else
        if isvalid(h)         
            count = count + 1;
            t1 = toc(start);
            time = t1 - t2;
            hms = fix(mod(((N-count)*time), [0, 3600, 60]) ./ [3600, 60, 1]);
            waitbar(count/N, h, [num2str(count/N*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            t2 = toc(start);
        end
    end
end
