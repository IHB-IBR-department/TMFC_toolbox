function [sub_check] = tmfc_LSS(tmfc,start_sub)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% For each individual trial, the Least-Squares Separate (LSS) approach
% estimates a separate GLM with two regressors. The first regressor models
% the expected BOLD response to the current trial of interest, and the 
% second (nuisance) regressor models the BOLD response to all other trials
% (of interest and no interest). For trials of no interest (e.g., errors),
% individual GLMs are not estimated. Trials of no interest are used only
% for the second (nuisance) regressor.
%
% This function uses SPM.mat file (which contains the specification of the
% 1st-level GLM with canonical HRF) to specify and estimate 1st-level GLMs
% for each individual trial of interest (LSS approach).
%
% Note: If the original GLMs contain parametric or time modulators, they
% will be removed from the LSS GLMs. If the original GLMs contain time and
% dispersion derivatives, they will be removed from the LSS GLMs.
%
% FORMAT [sub_check] = tmfc_LSS(tmfc)
% Run a function starting from the first subject in the list.
%
%   tmfc.subjects.path     - Paths to individual SPM.mat files
%   tmfc.subjects.name     - Subject names within the TMFC project
%                           ('Subject_XXXX' naming will be used by default)
%   tmfc.project_path      - Path where all results will be saved
%   tmfc.defaults.parallel - 0 or 1 (sequential or parallel computing)
%   tmfc.defaults.maxmem   - e.g. 2^31 = 2GB (how much RAM can be used at
%                            the same time during GLM estimation)
%   tmfc.defaults.resmem   - true or false (store temporaty files during
%                            GLM estimation in RAM or on disk)
%
%   tmfc.LSS.conditions           - List of conditions of interest
%   tmfc.LSS.conditions.sess      - Session number 
%                                   (as specified in SPM.Sess)
%   tmfc.LSS.conditions.number    - Condition number
%                                   (as specified in SPM.Sess.U)
%   tmfc.LSS.conditions.name      - Condition name
%                                   (as specified in SPM.Sess.U.name)
%   tmfc.LSS.conditions.file_name - Condition-specific file names:
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
%   tmfc.LSS.conditions(1).sess   = 1;     
%   tmfc.LSS.conditions(1).number = 1; 
%   tmfc.LSS.conditions(1).name = 'Cond_A';
%   tmfc.LSS.conditions(1).file_name = '[Sess_1]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS.conditions(2).sess   = 1;
%   tmfc.LSS.conditions(2).number = 2;
%   tmfc.LSS.conditions(2).name = 'Cond_B';
%   tmfc.LSS.conditions(2).file_name = '[Sess_1]_[Cond_2]_[Cond_B]';  
%
%   tmfc.LSS.conditions(3).sess   = 2;
%   tmfc.LSS.conditions(3).number = 1;
%   tmfc.LSS.conditions(3).name = 'Cond_A';
%   tmfc.LSS.conditions(3).file_name = '[Sess_2]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS.conditions(4).sess   = 2;
%   tmfc.LSS.conditions(4).number = 2;
%   tmfc.LSS.conditions(4).name = 'Cond_B';
%   tmfc.LSS.conditions(4).file_name = '[Sess_2]_[Cond_2]_[Cond_B]'; 
%
% FORMAT [sub_check] = tmfc_LSS(tmfc, start_sub)
% Run the function starting from a specific subject in the path list.
%
%   tmfc                   - As above
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
    set(main_GUI.TMFC_GUI_S6,'String', 'Updating...','ForegroundColor',[0.772, 0.353, 0.067]);       
end

spm('defaults','fmri');
spm_jobman('initcfg');

nSub = length(tmfc.subjects);
cond_list = tmfc.LSS.conditions;
nCond = length(cond_list);
sess = []; sess_num = []; nSess = [];
for iCond = 1:nCond
    sess(iCond) = cond_list(iCond).sess;
end
sess_num = unique(sess);
nSess = length(sess_num);

% Initialize waitbar
w = waitbar(0,'Please wait...','Name','LSS regression','Tag','tmfc_waitbar');
start_time = tic;
count_sub = 1;
cleanupObj = onCleanup(@unfreeze_after_ctrl_c);

% Loop through subjects
for iSub = start_sub:nSub
    SPM = load(tmfc.subjects(iSub).path).SPM;

    % Check if SPM.mat has concatenated sessions 
    % (if spm_fmri_concatenate.m sript was used)
    if size(SPM.nscan,2) == size(SPM.Sess,2)
        SPM_concat(iSub) = 0;
    else
        SPM_concat(iSub) = 1;
    end
    concat(iSub).scans = SPM.nscan;
    
    if isdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name))
        rmdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name),'s');
    end

    if ~isdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name))
        mkdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,'Betas'));
        mkdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,'GLM_batches'));
    end

    % Loop through sessions
    for jSess = 1:nSess       
        
        % Trials of interest
        nTrial = 0;
        ons_of_int = [];
        dur_of_int = [];
        cond_of_int = [];
        trial.cond = [];
        trial.number = [];
        for kCond = 1:nCond
            if cond_list(kCond).sess == sess_num(jSess)
                nTrial = nTrial + length(SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons);
                ons_of_int = [ons_of_int; SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons];
                dur_of_int = [dur_of_int; SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).dur];
                cond_of_int = [cond_of_int cond_list(kCond).number];
                trial.cond = [trial.cond; repmat(cond_list(kCond).number,length(SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons),1)];
                trial.number = [trial.number; (1:length(SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons))'];
            end
        end

        all_trials_number = (1:nTrial)';  

        % Trials of no interest
        cond_of_no_int = setdiff((1:length(SPM.Sess(sess_num(jSess)).U)),cond_of_int);
        ons_of_no_int = [];
        dur_of_no_int = [];
        for kCondNoInt = 1:length(cond_of_no_int)
            ons_of_no_int = [ons_of_no_int; SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).ons];
            dur_of_no_int = [dur_of_no_int; SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).dur];
        end
        
        % Loop through trials of interest
        for kTrial = 1:nTrial

            if isdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]))
                rmdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]),'s');
            end
            
            mkdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]));
            matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)])};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = SPM.xBF.UNITS;
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = SPM.xY.RT;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = SPM.xBF.T;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SPM.xBF.T0;
                        
            % Functional images
            if SPM_concat(iSub) == 0
                for image = 1:SPM.nscan(sess_num(jSess))
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{image,1} = [SPM.xY.VY(SPM.Sess(sess_num(jSess)).row(image)).fname ',' ...
                                                                              num2str(SPM.xY.VY(SPM.Sess(sess_num(jSess)).row(image)).n(1))];       
                end
            else
                for image = 1:size(SPM.xY.VY,1)
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{image,1} = [SPM.xY.VY(SPM.Sess(jSess).row(image)).fname ',' ...
                                                                              num2str(SPM.xY.VY(SPM.Sess(jSess).row(image)).n(1))];
                end
            end
    
            % Current trial vs all other trials (of interest and no interrest)
            current_trial_ons = ons_of_int(kTrial);
            current_trial_dur = dur_of_int(kTrial);
            other_trials = all_trials_number(all_trials_number~=kTrial);
            other_trials_ons = [ons_of_int(other_trials); ons_of_no_int];
            other_trials_dur = [dur_of_int(other_trials); dur_of_no_int];
            
            % Conditions
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Current_trial';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = current_trial_ons;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = current_trial_dur;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Other_trials';
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = other_trials_ons;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = other_trials_dur;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;

            % Confounds       
            for conf = 1:length(SPM.Sess(sess_num(jSess)).C.name)
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(conf).name = SPM.Sess(sess_num(jSess)).C.name{1,conf};
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(conf).val = SPM.Sess(sess_num(jSess)).C.C(:,conf);
            end   

            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    
            % HPF, HRF, mask 
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = SPM.xX.K(sess_num(jSess)).HParam;    
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = SPM.xGX.iGXcalc;
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = SPM.xM.gMT;           
            
            try
                matlabbatch{1}.spm.stats.fmri_spec.mask = {SPM.xM.VM.fname};
            catch
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            end
            
            if strcmp(SPM.xVi.form,'i.i.d') || strcmp(SPM.xVi.form,'none')
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'None';
            elseif strcmp(SPM.xVi.form,'fast') || strcmp(SPM.xVi.form,'FAST')
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST';
            else
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            end
        
            if strcmp(SPM.xVi.form,'wls')
                rWLS(iSub) = 1;
            else
                rWLS(iSub) = 0;
            end

            matlabbatch_2{1}.spm.stats.fmri_est.spmmat(1) = {fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'SPM.mat')};
            matlabbatch_2{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch_2{1}.spm.stats.fmri_est.method.Classical = 1;

            batch{kTrial} = matlabbatch;
            batch_2{kTrial} = matlabbatch_2;
            clear matlabbatch matlabbatch_2 current* other*
        end
         
        % Sequential or parallel computing
        switch tmfc.defaults.parallel                                 
            % -------------------- Sequential computing -------------------
            case 0
                for kTrial = 1:nTrial                                          
                    % Specify LSS GLM
                    spm('defaults','fmri');
                    spm_jobman('initcfg');
                    spm_get_defaults('cmdline',true);
                    spm_get_defaults('stats.resmem',tmfc.defaults.resmem);
                    spm_get_defaults('stats.maxmem',tmfc.defaults.maxmem);
                    spm_get_defaults('stats.fmri.ufp',1);
                    spm_jobman('run',batch{kTrial});
                    % Concatenated sessions
                    if SPM_concat(iSub) == 1
                        spm_fmri_concatenate(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name, ...
                            ['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'SPM.mat'),concat(iSub).scans);
                    end
                    % Check for rWLS
                    if rWLS(iSub) == 0
                        spm_jobman('run', batch_2{kTrial});
                    else
                        tmfc_rwls_LSS(tmfc,sess_num,iSub,jSess,kTrial);
                    end

                    % Save individual trial beta image
                    copyfile(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'beta_0001.nii'),...
                        fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,'Betas', ...
                        ['Beta_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name(1)),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].nii']));

                    % Save GLM_batch.mat file
                    tmfc_parsave_batch(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,'GLM_batches',...
                        ['GLM_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name(1)),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].mat']),batch{kTrial});

                    % Remove temporary LSS directory
                    rmdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]),'s');
                    
                    pause(0.01)
                    condition(trial.cond(kTrial)).trials(trial.number(kTrial)) = 1;
                end
            % --------------------- Parallel computing --------------------      
            case 1
                try
                    parpool;
                    figure(findobj('Tag','TMFC_GUI'));
                end

                parfor kTrial = 1:nTrial                   
                    % Specify LSS GLM
                    spm('defaults','fmri');
                    spm_jobman('initcfg');
                    spm_get_defaults('cmdline',true);
                    spm_get_defaults('stats.resmem',tmfc.defaults.resmem);
                    spm_get_defaults('stats.maxmem',tmfc.defaults.maxmem);
                    spm_get_defaults('stats.fmri.ufp',1);
                    spm_jobman('run',batch{kTrial});
                    % Concatenated sessions
                    if SPM_concat(iSub) == 1
                        spm_fmri_concatenate(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name, ...
                            ['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'SPM.mat'),concat(iSub).scans);
                    end
                    % Check for rWLS
                    if rWLS(iSub) == 0
                        spm_jobman('run', batch_2{kTrial});
                    else
                        tmfc_rwls_LSS(tmfc,sess_num,iSub,jSess,kTrial);
                    end

                    % Save individual trial beta image
                    copyfile(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'beta_0001.nii'),...
                        fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,'Betas', ...
                        ['Beta_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name(1)),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].nii']));

                    % Save GLM_batch.mat file
                    tmfc_parsave_batch(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,'GLM_batches',...
                        ['GLM_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name(1)),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].mat']),batch{kTrial});

                    % Remove temporary LSS directory
                    rmdir(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name,['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]),'s');
                    
                    trials(kTrial) = 1;                   
                end

                for kTrial = 1:nTrial
                    condition(trial.cond(kTrial)).trials(trial.number(kTrial)) = trials(kTrial);
                end
                clear trials
        end

        sub_check(iSub).session(sess_num(jSess)).condition = condition;
        pause(0.0001)

        clear nTrial ons* dur* cond_of_int cond_of_no_int trial all_trials_number condition 
    end
    
    % Update main TMFC GUI
    try                                                             
        main_GUI = guidata(findobj('Tag','TMFC_GUI'));
        set(main_GUI.TMFC_GUI_S6,'String', strcat(num2str(iSub), '/', num2str(nSub), ' done'),'ForegroundColor',[0.219, 0.341, 0.137]);
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

    clear SPM batch

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

% Estimate rWLS LSS model
function tmfc_rwls_LSS(tmfc,sess_num,iSub,jSess,kTrial)
    SPM_LSS = load(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name, ...
          ['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'SPM.mat')).SPM;
    SPM_LSS.xVi.form = 'wls';
    nScan = sum(SPM_LSS.nscan);
    for jScan = 1:nScan
        SPM_LSS.xVi.Vi{jScan} = sparse(nScan,nScan);
        SPM_LSS.xVi.Vi{jScan}(jScan,jScan) = 1;
    end
    original_dir = pwd;
    cd(fullfile(tmfc.project_path,'LSS_regression',tmfc.subjects(iSub).name, ...
          ['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]));
    tmfc_spm_rwls_spm(SPM_LSS);
    cd(original_dir);
end

% Save batches in parallel mode
function tmfc_parsave_batch(fname,matlabbatch)
	save(fname, 'matlabbatch')
end
