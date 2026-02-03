function [preDVARS,postDVARS] = tmfc_calculate_DVARS(FD,SPM_paths,options,masks,output_paths)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
%
% (1) Calculates Derivative of root mean square VARiance over voxelS (DVARS).
% DVARS is the root mean square (RMS) of the temporal derivative of
% BOLD time series within the GM mask (Muschelli et al., 2014).
%
% (2) Computes FD-DVARS correlations. 
%
% (3) Computes task-DVARS correlations (Pearson's r) per session using task
% regressors from the SPM design matrix, with summary statistics (mean / max |r|).
%
% DVARS is computed before and after noise regression.
% (original vs. updated GLM, respectively).
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

% Specify filename for DVARS after denoising
%--------------------------------------------------------------------------
if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 0
    aCompCor_fname = ['[aCompCor_' num2str(options.aCompCor(1)) 'WM_' num2str(options.aCompCor(2)) 'CSF]'];
elseif (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 1
    aCompCor_fname = ['[aCompCor_' num2str(options.aCompCor(1)) 'WM_' num2str(options.aCompCor(2)) 'CSF_Ort]'];
elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 0
    aCompCor_fname = '[aCompCor50]';
elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 1
    aCompCor_fname = ['[aCompCor50_Ort]'];
end

DVARS2_fname = ['DVARS_[' options.motion ']'];
if options.rWLS == 1; DVARS2_fname = strcat(DVARS2_fname,'_[rWLS]'); end
if options.spikereg == 1; DVARS2_fname = strcat(DVARS2_fname,['_[SpikeReg_' num2str(options.spikeregFDthr) 'mm]']); end
if ~strcmpi(options.GSR,'none'); DVARS2_fname = strcat(DVARS2_fname,['_[' options.GSR ']']); end
if ~strcmpi(options.WM_CSF,'none'); DVARS2_fname = strcat(DVARS2_fname,['_[' options.WM_CSF ']']); end
if sum(options.aCompCor)~=0; DVARS2_fname = strcat(DVARS2_fname,['_' aCompCor_fname]); end
DVARS2_fname = strcat(DVARS2_fname,'.mat');

% Calculate DVARS
%--------------------------------------------------------------------------
w = waitbar(0,'Please wait...','Name','Calculating DVARS');
for iSub = 1:length(SPM_paths)

    % GM mask path
    %----------------------------------------------------------------------
    GM_path = masks.GM{iSub};
    if isempty(GM_path) || ~exist(GM_path,'file')
        error('GM mask not found for subject %d.', iSub);
    end

    % DVARS (before denoising) (use original GLM)
    %----------------------------------------------------------------------
    if ~exist(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat'),'file')
        SPM = load(SPM_paths{iSub}).SPM;
        GM = tmfc_extract_GM_residuals(SPM,GM_path,0);
        FD_DVARS_corr = [];
        for jSess = 1:size(SPM.Sess,2)
            DVARS.Sess(jSess).DVARS_ts = sqrt(mean([zeros(1,size(GM(SPM.Sess(jSess).row,:),2)); diff(GM(SPM.Sess(jSess).row,:))].^2, 2));
            DVARS.Sess(jSess).DVARS_ts(1:3) = NaN; DVARS.Sess(jSess).DVARS_ts(end) = NaN;
            DVARS.Sess(jSess).FD_DVARS_corr = tmfc_corr(DVARS.Sess(jSess).DVARS_ts(4:end-1),FD(iSub).Sess(jSess).FD_ts(4:end-1));
            FD_DVARS_corr = [FD_DVARS_corr DVARS.Sess(jSess).FD_DVARS_corr];
            % -------------------------------------------------------------
            % Task-DVARS correlations
            % -------------------------------------------------------------
            task_cols = [];
            if isfield(SPM.Sess(jSess),'Fc')
                Fc = SPM.Sess(jSess).Fc;
            else
                Fc = [];
            end

            if ~isempty(Fc)
                for kCond = 1:length(Fc)
                    task_cols = [task_cols SPM.Sess(jSess).col(Fc(kCond).i)];
                end
            end

            if isempty(task_cols)
                DVARS.Sess(jSess).task_names = {};
                DVARS.Sess(jSess).taskDVARS_corr = [];
                DVARS.Sess(jSess).taskDVARS_corr_mean = NaN;
                DVARS.Sess(jSess).taskDVARS_corr_maxabs = NaN;
                DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = '';
            else
                DVARS.Sess(jSess).task_names = SPM.xX.name(task_cols);

                Xtask = SPM.xX.X(SPM.Sess(jSess).row, task_cols);
                y = DVARS.Sess(jSess).DVARS_ts(4:end-1);

                taskDV = nan(1, size(Xtask,2));
                for k = 1:size(Xtask,2)
                    x = Xtask(4:end-1, k);
                    taskDV(k) = tmfc_corr(y, x);
                end

                DVARS.Sess(jSess).taskDVARS_corr = taskDV;

                % Session summaries
                r = taskDV(:);
                ok = isfinite(r);
                if any(ok)
                    DVARS.Sess(jSess).taskDVARS_corr_mean = mean(r(ok));

                    [mx, imx] = max(abs(r(ok)));
                    idx_ok = find(ok);
                    imx = idx_ok(imx);

                    DVARS.Sess(jSess).taskDVARS_corr_maxabs = mx;
                    DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = DVARS.Sess(jSess).task_names{imx};
                else
                    DVARS.Sess(jSess).taskDVARS_corr_mean = NaN;
                    DVARS.Sess(jSess).taskDVARS_corr_maxabs = NaN;
                    DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = '';
                end

                clear Xtask taskDV r ok mx imx idx_ok x y
            end
        end

        % -----------------------------------------------------------------
        % Subject-level summaries across sessions 
        % -----------------------------------------------------------------
        DVARS.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
        [~, idx] = max(abs(FD_DVARS_corr));
        DVARS.Max_FD_DVARS_corr = FD_DVARS_corr(idx);
        clear idx
        
        % Task-DVARS correlations
        sess_mean = nan(1, length(DVARS.Sess));
        sess_max  = nan(1, length(DVARS.Sess));
        sess_name = cell(1, length(DVARS.Sess));

        for jSess = 1:length(DVARS.Sess)
            sess_mean(jSess) = DVARS.Sess(jSess).taskDVARS_corr_mean;
            sess_max(jSess)  = DVARS.Sess(jSess).taskDVARS_corr_maxabs;
            sess_name{jSess} = DVARS.Sess(jSess).taskDVARS_corr_maxabs_name;
        end

        okm = isfinite(sess_mean);
        if any(okm)
            DVARS.taskDVARS_corr_mean = mean(sess_mean(okm));
        else
            DVARS.taskDVARS_corr_mean = NaN;
        end

        okx = isfinite(sess_max);
        if any(okx)
            [DVARS.taskDVARS_corr_maxabs, imax] = max(sess_max(okx));
            idx_ok = find(okx);
            imax = idx_ok(imax);
            DVARS.taskDVARS_corr_maxabs_name = sess_name{imax};
        else
            DVARS.taskDVARS_corr_maxabs = NaN;
            DVARS.taskDVARS_corr_maxabs_name = '';
        end

        clear sess_mean sess_max sess_name okm okx imax idx_ok
        
        % Save
        save(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat'),'DVARS');
        clear SPM GM FD_DVARS_corr DVARS
    end

    % DVARS (after denoising) (use updated GLM)
    %----------------------------------------------------------------------
    if ~exist(fullfile(masks.glm_paths{iSub},DVARS2_fname),'file')
        if strcmpi(DVARS2_fname,'DVARS_[6HMP].mat')
            SPM = load(SPM_paths{iSub}).SPM;
            GM = tmfc_extract_GM_residuals(SPM,GM_path,1);
            FD_DVARS_corr = [];
            for jSess = 1:size(SPM.Sess,2)
                DVARS.Sess(jSess).DVARS_ts = sqrt(mean([zeros(1,size(GM(SPM.Sess(jSess).row,:),2)); diff(GM(SPM.Sess(jSess).row,:))].^2, 2));
                DVARS.Sess(jSess).DVARS_ts(1:3) = NaN; DVARS.Sess(jSess).DVARS_ts(end) = NaN;
                DVARS.Sess(jSess).FD_DVARS_corr = tmfc_corr(DVARS.Sess(jSess).DVARS_ts(4:end-1),FD(iSub).Sess(jSess).FD_ts(4:end-1));
                FD_DVARS_corr = [FD_DVARS_corr DVARS.Sess(jSess).FD_DVARS_corr];
                % -------------------------------------------------------------
                % Task-DVARS correlations
                % -------------------------------------------------------------
                task_cols = [];
                if isfield(SPM.Sess(jSess),'Fc')
                    Fc = SPM.Sess(jSess).Fc;
                else
                    Fc = [];
                end

                if ~isempty(Fc)
                    for kCond = 1:length(Fc)
                        task_cols = [task_cols SPM.Sess(jSess).col(Fc(kCond).i)];
                    end
                end

                if isempty(task_cols)
                    DVARS.Sess(jSess).task_names = {};
                    DVARS.Sess(jSess).taskDVARS_corr = [];
                    DVARS.Sess(jSess).taskDVARS_corr_mean = NaN;
                    DVARS.Sess(jSess).taskDVARS_corr_maxabs = NaN;
                    DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = '';
                else
                    DVARS.Sess(jSess).task_names = SPM.xX.name(task_cols);

                    Xtask = SPM.xX.X(SPM.Sess(jSess).row, task_cols);
                    y = DVARS.Sess(jSess).DVARS_ts(4:end-1);

                    taskDV = nan(1, size(Xtask,2));
                    for k = 1:size(Xtask,2)
                        x = Xtask(4:end-1, k);
                        taskDV(k) = tmfc_corr(y, x);
                    end

                    DVARS.Sess(jSess).taskDVARS_corr = taskDV;

                    % Session summaries
                    r = taskDV(:);
                    ok = isfinite(r);
                    if any(ok)
                        DVARS.Sess(jSess).taskDVARS_corr_mean = mean(r(ok));

                        [mx, imx] = max(abs(r(ok)));
                        idx_ok = find(ok);
                        imx = idx_ok(imx);

                        DVARS.Sess(jSess).taskDVARS_corr_maxabs = mx;
                        DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = DVARS.Sess(jSess).task_names{imx};
                    else
                        DVARS.Sess(jSess).taskDVARS_corr_mean = NaN;
                        DVARS.Sess(jSess).taskDVARS_corr_maxabs = NaN;
                        DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = '';
                    end

                    clear Xtask taskDV r ok mx imx idx_ok x y
                end
            end

            % -------------------------------------------------------------
            % Subject-level summaries across sessions
            % -------------------------------------------------------------
            DVARS.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
            [~, idx] = max(abs(FD_DVARS_corr));
            DVARS.Max_FD_DVARS_corr = FD_DVARS_corr(idx);
            clear idx

            % Task-DVARS correlations
            sess_mean = nan(1, length(DVARS.Sess));
            sess_max  = nan(1, length(DVARS.Sess));
            sess_name = cell(1, length(DVARS.Sess));

            for jSess = 1:length(DVARS.Sess)
                sess_mean(jSess) = DVARS.Sess(jSess).taskDVARS_corr_mean;
                sess_max(jSess)  = DVARS.Sess(jSess).taskDVARS_corr_maxabs;
                sess_name{jSess} = DVARS.Sess(jSess).taskDVARS_corr_maxabs_name;
            end

            okm = isfinite(sess_mean);
            if any(okm)
                DVARS.taskDVARS_corr_mean = mean(sess_mean(okm));
            else
                DVARS.taskDVARS_corr_mean = NaN;
            end

            okx = isfinite(sess_max);
            if any(okx)
                [DVARS.taskDVARS_corr_maxabs, imax] = max(sess_max(okx));
                idx_ok = find(okx);
                imax = idx_ok(imax);
                DVARS.taskDVARS_corr_maxabs_name = sess_name{imax};
            else
                DVARS.taskDVARS_corr_maxabs = NaN;
                DVARS.taskDVARS_corr_maxabs_name = '';
            end

            clear sess_mean sess_max sess_name okm okx imax idx_ok

            % Save
            save(fullfile(masks.glm_paths{iSub},DVARS2_fname),'DVARS');
            clear SPM GM FD_DVARS_corr DVARS

        elseif exist(fullfile(output_paths{iSub},'SPM.mat'),'file')
            SPM = load(fullfile(output_paths{iSub},'SPM.mat')).SPM;
            GM = tmfc_extract_GM_residuals(SPM,GM_path,1);
            FD_DVARS_corr = [];
            for jSess = 1:size(SPM.Sess,2)
                DVARS.Sess(jSess).DVARS_ts = sqrt(mean([zeros(1,size(GM(SPM.Sess(jSess).row,:),2)); diff(GM(SPM.Sess(jSess).row,:))].^2, 2));
                DVARS.Sess(jSess).DVARS_ts(1:3) = NaN; DVARS.Sess(jSess).DVARS_ts(end) = NaN;
                DVARS.Sess(jSess).FD_DVARS_corr = tmfc_corr(DVARS.Sess(jSess).DVARS_ts(4:end-1),FD(iSub).Sess(jSess).FD_ts(4:end-1));
                FD_DVARS_corr = [FD_DVARS_corr DVARS.Sess(jSess).FD_DVARS_corr];
                % ---------------------------------------------------------
                % Task-DVARS correlations
                % ---------------------------------------------------------
                task_cols = [];
                if isfield(SPM.Sess(jSess),'Fc')
                    Fc = SPM.Sess(jSess).Fc;
                else
                    Fc = [];
                end

                if ~isempty(Fc)
                    for kCond = 1:length(Fc)
                        task_cols = [task_cols SPM.Sess(jSess).col(Fc(kCond).i)];
                    end
                end

                if isempty(task_cols)
                    DVARS.Sess(jSess).task_names = {};
                    DVARS.Sess(jSess).taskDVARS_corr = [];
                    DVARS.Sess(jSess).taskDVARS_corr_mean = NaN;
                    DVARS.Sess(jSess).taskDVARS_corr_maxabs = NaN;
                    DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = '';
                else
                    DVARS.Sess(jSess).task_names = SPM.xX.name(task_cols);

                    Xtask = SPM.xX.X(SPM.Sess(jSess).row, task_cols);
                    y = DVARS.Sess(jSess).DVARS_ts(4:end-1);

                    taskDV = nan(1, size(Xtask,2));
                    for k = 1:size(Xtask,2)
                        x = Xtask(4:end-1, k);
                        taskDV(k) = tmfc_corr(y, x);
                    end

                    DVARS.Sess(jSess).taskDVARS_corr = taskDV;

                    % Session summaries
                    r = taskDV(:);
                    ok = isfinite(r);
                    if any(ok)
                        DVARS.Sess(jSess).taskDVARS_corr_mean = mean(r(ok));

                        [mx, imx] = max(abs(r(ok)));
                        idx_ok = find(ok);
                        imx = idx_ok(imx);

                        DVARS.Sess(jSess).taskDVARS_corr_maxabs = mx;
                        DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = DVARS.Sess(jSess).task_names{imx};
                    else
                        DVARS.Sess(jSess).taskDVARS_corr_mean = NaN;
                        DVARS.Sess(jSess).taskDVARS_corr_maxabs = NaN;
                        DVARS.Sess(jSess).taskDVARS_corr_maxabs_name = '';
                    end

                    clear Xtask taskDV r ok mx imx idx_ok x y
                end
            end
            % -------------------------------------------------------------
            % Subject-level summaries across sessions
            % -------------------------------------------------------------
            DVARS.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
            [~, idx] = max(abs(FD_DVARS_corr));
            DVARS.Max_FD_DVARS_corr = FD_DVARS_corr(idx);
            clear idx

            % Task-DVARS correlations
            sess_mean = nan(1, length(DVARS.Sess));
            sess_max  = nan(1, length(DVARS.Sess));
            sess_name = cell(1, length(DVARS.Sess));

            for jSess = 1:length(DVARS.Sess)
                sess_mean(jSess) = DVARS.Sess(jSess).taskDVARS_corr_mean;
                sess_max(jSess)  = DVARS.Sess(jSess).taskDVARS_corr_maxabs;
                sess_name{jSess} = DVARS.Sess(jSess).taskDVARS_corr_maxabs_name;
            end

            okm = isfinite(sess_mean);
            if any(okm)
                DVARS.taskDVARS_corr_mean = mean(sess_mean(okm));
            else
                DVARS.taskDVARS_corr_mean = NaN;
            end

            okx = isfinite(sess_max);
            if any(okx)
                [DVARS.taskDVARS_corr_maxabs, imax] = max(sess_max(okx));
                idx_ok = find(okx);
                imax = idx_ok(imax);
                DVARS.taskDVARS_corr_maxabs_name = sess_name{imax};
            else
                DVARS.taskDVARS_corr_maxabs = NaN;
                DVARS.taskDVARS_corr_maxabs_name = '';
            end

            clear sess_mean sess_max sess_name okm okx imax idx_ok

            % Save
            save(fullfile(masks.glm_paths{iSub},DVARS2_fname),'DVARS');
            clear SPM GM FD_DVARS_corr DVARS
        else
            disp(['To calculate DVARS after denoising, estimate the updated GLM for: ' FD(iSub).Subject]);
        end
    end
    try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
end
try, if exist('w','var') && ishghandle(w), close(w); end, end

% Load DVARS
%--------------------------------------------------------------------------
preDVARS = []; postDVARS = [];
for iSub = 1:length(SPM_paths)
    preDVARS(iSub).DVARS = load(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat')).DVARS;
    try; postDVARS(iSub).DVARS = load(fullfile(masks.glm_paths{iSub},DVARS2_fname)).DVARS; end
end
end
