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
%   tmfc.subjects.name            - Subject names within the TMFC project
%                                   ('Subject_XXXX' naming will be used by default)
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
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions        - List of conditions of interest
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.sess   - Session number (as specified in SPM.Sess)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.number - Condition number (as specified in SPM.Sess.U)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.pmod   - Parametric/Time modulator number (see SPM.Sess.U.P)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.name   - Condition name (as specified in SPM.Sess.U.name(kPmod))
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.file_name - Condition-specific file names:
%   (['[Sess_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
%    regexprep(char(SPM.Sess(iSess).U(jCond).name(1)),' ','_') ']'];)
%
% Session number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two sessions. Both sessions 
% contains three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified (see tmfc_conditions_GUI, nested function:
% [cond_list] = generate_conditions(SPM_path)):
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).sess   = 1;   
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).number = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).name = 'Cond_A'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).file_name = '[Sess_1]_[Cond_1]_[Cond_A]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).sess   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).name = 'Cond_B'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).file_name = '[Sess_1]_[Cond_1]_[Cond_B]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).number = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).name = 'Cond_A'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).file_name = '[Sess_2]_[Cond_1]_[Cond_A]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).name = 'Cond_B'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).file_name = '[Sess_2]_[Cond_2]_[Cond_B]';
%
% If GLMs contain parametric or time modulators, add the following fields:
% e.g. first modulator for fourth condition:
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).sess   = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).pmod   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).name = 'Cond_BxModulator1^1';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).file_name = '[Sess_2]_[Cond_2]_[Cond_BxModulator1^1]'; 
% e.g. second modulator for second condition:
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).sess   = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).number = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).pmod = 3; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).name = 'Cond_BxModulator2^1'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).file_name = '[Sess_2]_[Cond_2]_[Cond_BxModulator2^1]'; 
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

try
    main_GUI = guidata(findobj('Tag','TMFC_GUI'));                           
    set(main_GUI.TMFC_GUI_S3,'String', 'Updating...','ForegroundColor',[0.772, 0.353, 0.067]);       
end

nSub = length(tmfc.subjects);
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
cond_list = tmfc.ROI_set(ROI_set_number).gPPI.conditions;
nCond = length(cond_list);
sess = []; sess_num = []; nSess = [];
for iCond = 1:nCond
    sess(iCond) = cond_list(iCond).sess;
end
sess_num = unique(sess);
nSess = length(sess_num);
sub_check = zeros(1,nSub);
if start_sub > 1
    sub_check(1:start_sub) = 1;
end

% Initialize waitbar
w = waitbar(0,'Please wait...','Name','VOI time-series extraction','Tag','tmfc_waitbar');
start_time = tic;
count_sub = 1;
cleanupObj = onCleanup(@unfreeze_after_ctrl_c);

spm('defaults','fmri');
spm_jobman('initcfg');

for iSub = start_sub:nSub
    % Calculate F-contrast for all conditions of interest
    SPM = load(tmfc.subjects(iSub).path).SPM;
    matlabbatch{1}.spm.stats.con.spmmat = {tmfc.subjects(iSub).path};
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'F_conditions_of_interest';
    cond_col = [];
    for iCond = 1:length(cond_list)
        FCi = [];
        FCi = SPM.Sess(cond_list(iCond).sess).Fc(cond_list(iCond).number).i; 
        try
            FCp = []; 
            FCp = SPM.Sess(cond_list(iCond).sess).Fc(cond_list(iCond).number).p; 
            FCi = FCi(FCp==cond_list(iCond).pmod);
        end
        cond_col = [cond_col SPM.Sess(cond_list(iCond).sess).col(FCi)];
    end 
    weights = zeros(length(cond_col),size(SPM.xX.X,2));
    for iCond = 1:length(cond_col)
        weights(iCond,cond_col(iCond)) = 1;
    end
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = weights;
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    spm_get_defaults('cmdline',true);
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    SPM = load(tmfc.subjects(iSub).path).SPM;

    if isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',tmfc.subjects(iSub).name))
        rmdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',tmfc.subjects(iSub).name),'s');
    end

    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',tmfc.subjects(iSub).name))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs',tmfc.subjects(iSub).name));
    end
    
    for jSess = 1:nSess
        for kROI = 1:nROI
            matlabbatch{1}.spm.util.voi.spmmat = {tmfc.subjects(iSub).path};
            matlabbatch{1}.spm.util.voi.adjust = length(SPM.xCon);
            matlabbatch{1}.spm.util.voi.session = sess_num(jSess);
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
                    movefile(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs', ... 
                                      tmfc.subjects(iSub).name,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']));
                    if exist(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']),'file')
                        delete(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']));
                    else
                        delete(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
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
                    movefile(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']), ...
                             fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs', ... 
                                      tmfc.subjects(iSub).name,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '.mat']));
                    if exist(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']),'file')
                        delete(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(sess_num(jSess)) '_eigen.nii']));
                    else
                        delete(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_eigen.nii']));
                    end
                    delete(fullfile(SPM.swd,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_mask.nii']));
                end
        end

        clear batch
    end
    
    sub_check(iSub) = 1;

    % Update main TMFC GUI
    try  
        main_GUI = guidata(findobj('Tag','TMFC_GUI'));                        
        set(main_GUI.TMFC_GUI_S3,'String', strcat(num2str(iSub), '/', num2str(nSub), ' done'),'ForegroundColor',[0.219, 0.341, 0.137]);    
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
