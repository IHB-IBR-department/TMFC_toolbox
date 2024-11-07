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
% FORMAT [sub_check] = tmfc_LSS(tmfc)
% Run a function starting from the first subject in the list.
%
%   tmfc.subjects.path     - Paths to individual SPM.mat files
%   tmfc.project_path      - Path where all results will be saved
%   tmfc.defaults.parallel - 0 or 1 (sequential or parallel computing)
%   tmfc.defaults.maxmem   - e.g. 2^31 = 2GB (how much RAM can be used at
%                            the same time during GLM estimation)
%   tmfc.defaults.resmem   - true or false (store temporaty files during
%                            GLM estimation in RAM or on disk)
%
%   tmfc.LSS.conditions        - List of conditions of interest
%   tmfc.LSS.conditions.sess   - Session number (as specified in SPM.Sess)
%   tmfc.LSS.conditions.number - Condition number (as specified in SPM.Sess.U)
%
% Session number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two sessions. Both sessions 
% contains three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified:
%
%   tmfc.LSS.conditions(1).sess   = 1;   
%   tmfc.LSS.conditions(1).number = 1; - "Cond A", 1st session
%   tmfc.LSS.conditions(2).sess   = 1;
%   tmfc.LSS.conditions(2).number = 2; - "Cond B", 1st session
%   tmfc.LSS.conditions(3).sess   = 2;
%   tmfc.LSS.conditions(3).number = 1; - "Cond A", 2nd session
%   tmfc.LSS.conditions(4).sess   = 2;
%   tmfc.LSS.conditions(4).number = 2; - "Cond B", 2nd session
%
% FORMAT [sub_check] = tmfc_LSS(tmfc, start_sub)
% Run the function starting from a specific subject in the path list.
%
%   tmfc                   - As above
%   start_sub              - Subject number on the path list to start with
%
% =========================================================================
%
% Copyright (C) 2023 Ruslan Masharipov
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

% Initialize waitbar for sequential or parallel computations
switch tmfc.defaults.parallel
    case 0    % Sequential
        w = waitbar(0,'Please wait...','Name','LSS regression','Tag','tmfc_waitbar');
        cleanupObj = onCleanup(@unfreeze_after_ctrl_c);
    case 1    % Parallel       
        try   % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;            
            w = waitbar(0,'Please wait...','Name','LSS regression','Tag','tmfc_waitbar');
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

% Loop through subjects
for iSub = start_sub:nSub
    tic
    SPM = load(tmfc.subjects(iSub).path);
    
    if isdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')]))
        rmdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')]),'s');
    end

    if ~isdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')]))
        mkdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],'Betas'));
        mkdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],'GLM_batches'));
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
                nTrial = nTrial + length(SPM.SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons);
                ons_of_int = [ons_of_int; SPM.SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons];
                dur_of_int = [dur_of_int; SPM.SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).dur];
                cond_of_int = [cond_of_int cond_list(kCond).number];
                trial.cond = [trial.cond; repmat(cond_list(kCond).number,length(SPM.SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons),1)];
                trial.number = [trial.number; (1:length(SPM.SPM.Sess(sess_num(jSess)).U(cond_list(kCond).number).ons))'];
            end
        end

        all_trials_number = (1:nTrial)';  

        % Trials of no interest
        cond_of_no_int = setdiff((1:length(SPM.SPM.Sess(sess_num(jSess)).U)),cond_of_int);
        ons_of_no_int = [];
        dur_of_no_int = [];
        for kCondNoInt = 1:length(cond_of_no_int)
            ons_of_no_int = [ons_of_no_int; SPM.SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).ons];
            dur_of_no_int = [dur_of_no_int; SPM.SPM.Sess(sess_num(jSess)).U(cond_of_no_int(kCondNoInt)).dur];
        end
        
        % Loop through trials of interest
        for kTrial = 1:nTrial

            if isdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]))
                rmdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]),'s');
            end
            
            mkdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]));
            matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)])};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = SPM.SPM.xBF.UNITS;
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = SPM.SPM.xY.RT;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = SPM.SPM.xBF.T;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SPM.SPM.xBF.T0;
                        
            % Functional images
            for image = 1:SPM.SPM.nscan(sess_num(jSess))
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans{image,1} = SPM.SPM.xY.VY(SPM.SPM.Sess(sess_num(jSess)).row(image)).fname;
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
            for conf = 1:length(SPM.SPM.Sess(sess_num(jSess)).C.name)
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(conf).name = SPM.SPM.Sess(sess_num(jSess)).C.name{1,conf};
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(conf).val = SPM.SPM.Sess(sess_num(jSess)).C.C(:,conf);
            end   

            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    
            % HPF, HRF, mask 
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = SPM.SPM.xX.K(sess_num(jSess)).HParam;    
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = SPM.SPM.xGX.iGXcalc;
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = SPM.SPM.xM.gMT;           
            
            try
                matlabbatch{1}.spm.stats.fmri_spec.mask = {SPM.SPM.xM.VM.fname};
            catch
                matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            end
            
            if strcmp(SPM.SPM.xVi.form,'i.i.d')
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'None';
            elseif strcmp(SPM.SPM.xVi.form,'AR(0.2)')
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            else
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST';
            end

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'SPM.mat')};
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            batch{kTrial} = matlabbatch;
            clear matlabbatch current* other*
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

                    % Save individual trial beta image
                    copyfile(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'beta_0001.nii'),...
                        fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],'Betas', ...
                        ['Beta_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].nii']));

                    % Save GLM_batch.mat file
                    tmfc_parsave_batch(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],'GLM_batches',...
                        ['GLM_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].mat']),batch{kTrial});

                    % Remove temporal LSS directory
                    rmdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]),'s');
                    
                    pause(0.01)
                    condition(trial.cond(kTrial)).trials(trial.number(kTrial)) = 1;
                end
            % --------------------- Parallel computing --------------------      
            case 1
                parfor kTrial = 1:nTrial                   
                    % Specify LSS GLM
                    spm('defaults','fmri');
                    spm_jobman('initcfg');
                    spm_get_defaults('cmdline',true);
                    spm_get_defaults('stats.resmem',tmfc.defaults.resmem);
                    spm_get_defaults('stats.maxmem',tmfc.defaults.maxmem);
                    spm_get_defaults('stats.fmri.ufp',1);
                    spm_jobman('run',batch{kTrial});

                    % Save individual trial beta image
                    copyfile(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)],'beta_0001.nii'),...
                        fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],'Betas', ...
                        ['Beta_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].nii']));

                    % Save GLM_batch.mat file
                    tmfc_parsave_batch(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],'GLM_batches',...
                        ['GLM_[Sess_' num2str(sess_num(jSess)) ']_[Cond_' num2str(trial.cond(kTrial)) ']_[' regexprep(char(SPM.SPM.Sess(sess_num(jSess)).U(trial.cond(kTrial)).name),' ','_') ']_[Trial_' num2str(trial.number(kTrial)) '].mat']),batch{kTrial});

                    % Remove temporal LSS directory
                    rmdir(fullfile(tmfc.project_path,'LSS_regression',['Subject_' num2str(iSub,'%04.f')],['LSS_Sess_' num2str(sess_num(jSess)) '_Trial_' num2str(kTrial)]),'s');
                    
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

    % Update waitbar for sequential or parallel computing
    switch(tmfc.defaults.parallel)
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

% Save batches in parallel mode
function tmfc_parsave_batch(fname,matlabbatch)
	save(fname, 'matlabbatch')
end

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
