function output_paths = tmfc_estimate_updated_GLMs(SPM_paths,masks,options)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% (1) Estimates updated GLMs with noise regressors. The noise regressors
%     and the updated model are saved in the TMFC_denoise subfolder.
%
% (2) Supports robust weighted least squares (rWLS) estimation.
%     It assumes that each image has its own variance parameter, i.e. some
%     scans may be disrupted by noise. With this option, SPM will 
%     estimates the noise variances in the first pass and then re-weights each
%     image by the inverse of the variance in the second pass.
%
% NOTE: The original first-level GLMs must not be estimated with rWLS.
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

% Paths to updated GLMs
%--------------------------------------------------------------------------
if strcmpi(options.motion,'6HMP') && options.spikereg == 0 && sum(options.aCompCor)==0 && strcmpi(options.WM_CSF,'none') && strcmpi(options.GSR,'none') && options.rWLS == 0 
    output_paths = [];
    shortwarn('Only 6 motion regressors are selected as noise regressors. No new models will be created.'); return;
elseif strcmpi(options.motion,'6HMP') && options.spikereg == 0 && sum(options.aCompCor)==0 && strcmpi(options.WM_CSF,'none') && strcmpi(options.GSR,'none') && options.rWLS == 1
    new_GLM_subfolder = ['GLM_[6HMP]_[rWLS]'];
elseif (~strcmpi(options.motion,'6HMP') || options.spikereg == 1 || options.rWLS == 1) && (sum(options.aCompCor)==0 && strcmpi(options.WM_CSF,'none') && strcmpi(options.GSR,'none')) 
    new_GLM_subfolder = ['GLM_[' options.motion ']'];
    if options.rWLS == 1; new_GLM_subfolder = strcat(new_GLM_subfolder,'_[rWLS]'); end
    if options.spikereg == 1; new_GLM_subfolder = strcat(new_GLM_subfolder,['_[SpikeReg_' num2str(options.spikeregFDthr,'%.2f') 'mm]']); end
elseif sum(options.aCompCor)~=0 || ~strcmpi(options.WM_CSF,'none') || ~strcmpi(options.GSR,'none')
    new_GLM_subfolder = ['[WM' num2str(round(options.WMmask.prob*100)) 'e' num2str(options.WMmask.erode) ...
        ']_[CSF' num2str(round(options.CSFmask.prob*100)) 'e' num2str(options.CSFmask.erode) ...
        ']_[GM' num2str(round(options.GMmask.prob*100)) 'd' num2str(options.GMmask.dilate) ']'];
    GLM_name = ['GLM_[' options.motion ']']; 
    if options.rWLS == 1; GLM_name = strcat(GLM_name,['_[rWLS]']); end
    if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 0
        aCompCor_fname = ['[aCompCor_' num2str(options.aCompCor(1)) 'WM_' num2str(options.aCompCor(2)) 'CSF]'];
    elseif (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 1
        aCompCor_fname = ['[aCompCor_' num2str(options.aCompCor(1)) 'WM_' num2str(options.aCompCor(2)) 'CSF_Ort]'];
    elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 0
        aCompCor_fname = '[aCompCor50]';
    elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 1
        aCompCor_fname = ['[aCompCor50_Ort]'];
    end
    if options.spikereg == 1; GLM_name = strcat(GLM_name,['_[SpikeReg_' num2str(options.spikeregFDthr,'%.2f') 'mm]']); end
    if ~strcmpi(options.GSR,'none'); GLM_name = strcat(GLM_name,['_[' options.GSR ']']); end
    if ~strcmpi(options.WM_CSF,'none'); GLM_name = strcat(GLM_name,['_[' options.WM_CSF ']']); end
    if sum(options.aCompCor)~=0; GLM_name = strcat(GLM_name,['_' aCompCor_fname]); end
    new_GLM_subfolder = fullfile(new_GLM_subfolder,GLM_name);
end

for iSub = 1:length(SPM_paths)
    old_GLM_subfolder = fileparts(SPM_paths{iSub});
    output_paths{iSub,1} = fullfile(old_GLM_subfolder,'TMFC_denoise',new_GLM_subfolder);
    clear old_GLM_subfolder
end

% Specify and estimate GLMs with noise regressors
%--------------------------------------------------------------------------
spm('defaults','fmri');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

jSub = 0;
for iSub = 1:length(SPM_paths)
    % Delete GLM folder, if model has not been estimated:
    if exist(fullfile(output_paths{iSub},'SPM.mat'),'file')
        SPMnew = load(fullfile(output_paths{iSub},'SPM.mat'));
        if ~isfield(SPMnew.SPM,'Vbeta')
            rmdir(output_paths{iSub},'s');
        end
        clear SPMnew
    end
 
    if ~exist(output_paths{iSub},'dir')
        mkdir(output_paths{iSub});
    end

    % Specify GLM with noise regressors
    %----------------------------------------------------------------------
    if ~exist(fullfile(output_paths{iSub},'SPM.mat'),'file')

        old_GLM_subfolder = fileparts(SPM_paths{iSub});

        % Load HMP
        %----------------------------------------------------------------------
        if strcmpi(options.motion,'12HMP')
            HMP = load(fullfile(old_GLM_subfolder,'TMFC_denoise','12HMP.mat')).HMP12;
        elseif strcmpi(options.motion,'24HMP')
            HMP = load(fullfile(old_GLM_subfolder,'TMFC_denoise','24HMP.mat')).HMP24;
        end
    
        % Load spike regressors
        %----------------------------------------------------------------------
        if options.spikereg == 1
            SpikeReg = load(fullfile(old_GLM_subfolder,'TMFC_denoise',sprintf('SpikeReg_[FDthr_%.2fmm].mat', options.spikeregFDthr))).SpikeReg;
        end
    
        % Load GSR
        %----------------------------------------------------------------------
        if strcmpi(options.GSR,'GSR')
            GSR = load(fullfile(masks.glm_paths{iSub},[options.GSR '.mat'])).GSR;
        elseif strcmpi(options.GSR,'2GSR')
            GSR = load(fullfile(masks.glm_paths{iSub},[options.GSR '.mat'])).GSR2;
        elseif strcmpi(options.GSR,'4GSR')
            GSR = load(fullfile(masks.glm_paths{iSub},[options.GSR '.mat'])).GSR4;    
        end
    
        % Load Phys
        %----------------------------------------------------------------------
        if strcmpi(options.WM_CSF,'2Phys')
            Phys = load(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat'])).Phys2;
        elseif strcmpi(options.WM_CSF,'4Phys')
            Phys = load(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat'])).Phys4;
        elseif strcmpi(options.WM_CSF,'8Phys')
            Phys = load(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat'])).Phys8;
        end
    
        % Load aCompCor
        %----------------------------------------------------------------------
        if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) 
            aCompCor = load(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat'])).aCompCor;
        elseif options.aCompCor(1) == 0.5
            aCompCor = load(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat'])).aCompCor50;
        end

        % Create SPM batch
        %------------------------------------------------------------------
        SPM = load(SPM_paths{iSub}).SPM;
    
        % Check if SPM.mat has concatenated sessions 
        % (if spm_fmri_concatenate.m script was used)
        if size(SPM.nscan,2) == size(SPM.Sess,2)
            SPM_concat(iSub) = -1;
        else
            SPM_concat(iSub) = 1;
        end
        concat(iSub).scans = SPM.nscan;

        matlabbatch{1}.spm.stats.fmri_spec.dir = output_paths(iSub);
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = SPM.xBF.UNITS;
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = SPM.xY.RT;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = SPM.xBF.T;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SPM.xBF.T0;

        %---------------------------------------------Loop through sessions
        for jSess = 1:length(SPM.Sess)

            % Functional images
            if SPM_concat(iSub) == -1
                for kImage = 1:SPM.nscan(jSess)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).scans{kImage,1} = [SPM.xY.VY(SPM.Sess(jSess).row(kImage)).fname ',' ...
                                                                                     num2str(SPM.xY.VY(SPM.Sess(jSess).row(kImage)).n(1))];
                end
            else
                for kImage = 1:size(SPM.xY.VY,1)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).scans{kImage,1} = [SPM.xY.VY(SPM.Sess(jSess).row(kImage)).fname ',' ...
                                                                                     num2str(SPM.xY.VY(SPM.Sess(jSess).row(kImage)).n(1))];
                end
            end
            
            % Conditions
            for kCond = 1:length(SPM.Sess(jSess).U)
                matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).name = SPM.Sess(jSess).U(kCond).name{1};
                matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).onset = SPM.Sess(jSess).U(kCond).ons;
                matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).duration = SPM.Sess(jSess).U(kCond).dur;
                matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).tmod = 0;
                if length(SPM.Sess(jSess).U(kCond).name)>1
                    for PM_number = 1:length(SPM.Sess(jSess).U(kCond).P)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).pmod(PM_number).name = SPM.Sess(jSess).U(kCond).P(PM_number).name;
                        matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).pmod(PM_number).param = SPM.Sess(jSess).U(kCond).P(PM_number).P;
                        matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).pmod(PM_number).poly = SPM.Sess(jSess).U(kCond).P(PM_number).h;
                    end
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).pmod = struct('name', {}, 'param', {}, 'poly', {});
                end
                matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).cond(kCond).orth = SPM.Sess(jSess).U(kCond).orth;
            end
            
            % Original confounds
            nOldConf = length(SPM.Sess(jSess).C.name);
            for kConf = 1:nOldConf
                matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).regress(kConf).name = SPM.Sess(jSess).C.name{1,kConf};
                matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).regress(kConf).val = SPM.Sess(jSess).C.C(:,kConf);
            end

            % Add new noise regressors-------------------------------------
            Conf = []; ConfName = {};
            % HMP
            if ~strcmpi(options.motion,'6HMP')
                Conf = [Conf HMP(jSess).Sess(:,7:end)]; % Start with 7th HMP, since 6HMP is already added
                C = cell(size(HMP(jSess).Sess(:,7:end),2),1); C(:) = {'HMP'};
                ConfName = [ConfName; C]; clear C
            end
            % SpikeReg
            if options.spikereg == 1
                Conf = [Conf SpikeReg(jSess).Sess];
                C = cell(size(SpikeReg(jSess).Sess,2),1); C(:) = {'SpikeReg'};
                ConfName = [ConfName; C]; clear C
            end
            % GSR
            if ~strcmpi(options.GSR,'none')
                Conf = [Conf GSR(jSess).Sess];
                C = cell(size(GSR(jSess).Sess,2),1); C(:) = {'GSR'};
                ConfName = [ConfName; C]; clear C
            end
            % Phys
            if ~strcmpi(options.WM_CSF,'none')
                Conf = [Conf Phys(jSess).Sess];
                C = cell(size(Phys(jSess).Sess,2),1); C(:) = {'WM_CSF'};
                ConfName = [ConfName; C]; clear C
            end
            % aCompCor
            if sum(options.aCompCor)~=0 
                Conf = [Conf aCompCor.Sess(jSess).WM_PCs aCompCor.Sess(jSess).CSF_PCs];
                C = cell(size(aCompCor.Sess(jSess).WM_PCs,2)+size(aCompCor.Sess(jSess).CSF_PCs,2),1); C(:) = {'aCompCor'};
                ConfName = [ConfName; C]; clear C
            end
            
            nNewConf = size(Conf,2);
            if nNewConf ~= 0
                for kConf = 1:nNewConf
                    matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).regress(kConf+nOldConf).name = ConfName{kConf};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).regress(kConf+nOldConf).val = Conf(:,kConf);
                end
            end

            clear Conf ConfName
            %--------------------------------------------------------------

            matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).multi_reg = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(jSess).hpf = SPM.xX.K(jSess).HParam;            
        end 
        %-----------------------------------------------End of session loop
    
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});

        % Basis functions
        if strcmpi(SPM.xBF.name,'hrf')
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        elseif strcmpi(SPM.xBF.name,'hrf (with time derivative)')
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
        elseif strcmpi(SPM.xBF.name,'hrf (with time and dispersion derivatives)')
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
        elseif strcmpi(SPM.xBF.name,'Fourier set')
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.length = SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.order = SPM.xBF.order;
        elseif strcmpi(SPM.xBF.name,'Fourier set (Hanning)')
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.length = SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.order = SPM.xBF.order;
        elseif strcmpi(SPM.xBF.name,'Gamma functions')
            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.length = SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.order = SPM.xBF.order;
        elseif strcmpi(SPM.xBF.name,'Finite Impulse Response')
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order = SPM.xBF.order;
        end
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = SPM.xGX.iGXcalc;
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = SPM.xM.gMT;
    
        try
            matlabbatch{1}.spm.stats.fmri_spec.mask = {SPM.xM.VM.fname};
        catch
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        end
    
        if strcmpi(SPM.xVi.form,'i.i.d') || strcmpi(SPM.xVi.form,'none')
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'None';
        elseif strcmpi(SPM.xVi.form,'fast') || strcmpi(SPM.xVi.form,'FAST')
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'FAST';
        else
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        end

        matlabbatch_2{1}.spm.stats.fmri_est.spmmat(1) = {fullfile(output_paths{iSub},'SPM.mat')};
        matlabbatch_2{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch_2{1}.spm.stats.fmri_est.method.Classical = 1;

        clear old_GLM_subfolder HMP SpikeReg GSR Phys aCompCor

        jSub = jSub + 1;
        save(fullfile(output_paths{iSub},'GLM_batch.mat'),'matlabbatch');
        batch{jSub} = matlabbatch;
        batch_2{jSub} = matlabbatch_2;
        clear matlabbatch matlabbatch_2 SPM 
    end  
end

% Run
%--------------------------------------------------------------------------
if exist('batch','var')
    concat(SPM_concat == 0) = [];
    SPM_concat = nonzeros(SPM_concat);
    if jSub == 1
        spm_jobman('run',batch{1});
        % Concatenated sessions
        if SPM_concat(1) == 1
            spm_fmri_concatenate(fullfile(batch{1}{1}.spm.stats.fmri_spec.dir,'SPM.mat'),concat(1).scans);
        end
        % WLS or rWLS
        if options.rWLS == 0
            spm_jobman('run',batch_2{1});
        else
            SPM = load(fullfile(output_paths{1},'SPM.mat')).SPM;
            SPM.xVi.form = 'wls';
            nScan = sum(SPM.nscan);
            for iScan = 1:nScan
                SPM.xVi.Vi{iScan} = sparse(nScan,nScan);
                SPM.xVi.Vi{iScan}(iScan,iScan) = 1;
            end
            original_dir = pwd;
            cd(output_paths{1});
            tmfc_spm_rwls_spm(SPM);
            cd(original_dir);
        end
    elseif jSub > 1
        % Waitbar
        tmfc_progress('init', jSub, 'Model estimation');

        % Detect PCT availability
        hasPCT = (exist('parfor','builtin')==5) && license('test','Distrib_Computing_Toolbox');

        % Parallel mode, PCT only 
        if options.parallel == 1 && hasPCT
            % Init parpool
            if isempty(gcp('nocreate')), parpool; end
            % DataQueue requires R2017a+ 
            D = [];
            try
                D = parallel.pool.DataQueue;
                afterEach(D, @(~) tmfc_progress('tick'));                     
            end
            % Init SPM
            spmSetup = parallel.pool.Constant(@() init_spm());
            % Run matlabbatches in parallel mode
            parfor iSub = 1:jSub 
                spmSetup.Value;
                spm_jobman('run',batch{iSub});
                % Concatenated sessions
                if SPM_concat(iSub) == 1
                    spm_fmri_concatenate(fullfile(batch{iSub}{1}.spm.stats.fmri_spec.dir,'SPM.mat'),concat(iSub).scans);
                end
                % Check for rWLS
                if options.rWLS == 0
                    spm_jobman('run',batch_2{iSub});
                else
                    tmfc_rwls(output_paths,iSub);
                end
                % Update waitbar
                try send(D,[]); end 
            end
        else 
            % Serial mode
            for iSub = 1:jSub
                spm_jobman('run', batch{iSub});
                % Concatenated sessions
                if SPM_concat(iSub) == 1
                    spm_fmri_concatenate(fullfile(batch{iSub}{1}.spm.stats.fmri_spec.dir,'SPM.mat'), concat(iSub).scans);
                end
                % Check for rWLS
                if options.rWLS == 0
                    spm_jobman('run', batch_2{iSub});
                else
                    tmfc_rwls(output_paths, iSub);
                end
                % Update waitbar
                tmfc_progress('tick');
            end
        end
        try tmfc_progress('done'); end % Close waitbar
    end
end
end

%==========================================================================

% Estimate rWLS model
%--------------------------------------------------------------------------
function tmfc_rwls(output_paths,iSub)
    SPM = load(fullfile(output_paths{iSub},'SPM.mat')).SPM;
    SPM.xVi.form = 'wls';
    nScan=sum(SPM.nscan);
    for iScan = 1:nScan
        SPM.xVi.Vi{iScan} = sparse(nScan,nScan);
        SPM.xVi.Vi{iScan}(iScan,iScan) = 1;
    end
    original_dir = pwd;
    cd(output_paths{iSub});
    tmfc_spm_rwls_spm(SPM);
    cd(original_dir);
end

% Short warning
%--------------------------------------------------------------------------
function shortwarn(msg)
    s = warning('query','backtrace');
    warning('off','backtrace');
    warning(msg);
    warning(s.state,'backtrace');
end

% SPM initialization
%--------------------------------------------------------------------------
function c = init_spm()
    spm('defaults','fmri');
    spm_jobman('initcfg');
    spm_get_defaults('cmdline', true);
    c = onCleanup(@() []); 
end