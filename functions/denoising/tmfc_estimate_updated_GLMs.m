function output_paths = tmfc_estimate_updated_GLMs(SPM_paths,masks,options)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% (1) Estimates updated GLMs with noise regressors. The noise regressors
%     and the updated model will be stored in the TMFC_denoise subfolder.
%
% (2) Can use robust weighted least squares (rWLS) for model estimation.
%     It assumes that each image has an own variance parameter, i.e. some
%     scans may be disrupted by noise. By choosing this option, SPM will 
%     estimte the noise variances in the first pass and then re-weight each
%     image by the inverse of the variance in the second pass.
%
%     NOTE: Original first-level GLMs should not use the rWLS etimation.
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

% Paths to updated GLMs
%--------------------------------------------------------------------------
if strcmp(options.motion,'6HMP') && options.spikereg == 0 && sum(options.aCompCor)==0 && strcmp(options.WM_CSF,'none') && strcmp(options.GSR,'none') && options.rWLS == 0
    output_paths = masks.glm_paths;
    warning('Only 6 motion regressors are selected as noise regressors. New models will not be created.'); return;
elseif (~strcmp(options.motion,'6HMP') || options.spikereg == 1 || options.rWLS == 1) && (sum(options.aCompCor)==0 && strcmp(options.WM_CSF,'none') && strcmp(options.GSR,'none')) 
    new_GLM_subfolder = ['GLM_[' options.motion ']'];
    if options.rWLS == 1; new_GLM_subfolder = strcat(new_GLM_subfolder,'_[rWLS]'); end
    if options.spikereg == 1; new_GLM_subfolder = strcat(new_GLM_subfolder,['_[SpikeReg_' num2str(options.spikeregFDthr) 'mm]']); end
elseif sum(options.aCompCor)~=0 || ~strcmp(options.WM_CSF,'none') || ~strcmp(options.GSR,'none')
    new_GLM_subfolder = ['[WM_' num2str(options.WMmask.prob) 'Prob_' num2str(options.WMmask.erode) ...
        'xErode]_[CSF_' num2str(options.CSFmask.prob) 'Prob_' num2str(options.CSFmask.erode) ...
        'xErode]_[GM_' num2str(options.GMmask.prob) 'Prob_' num2str(options.GMmask.dilate) 'xDilate]'];
    GLM_name = ['GLM_[' options.motion ']']; 
    if options.rWLS == 1; GLM_name = strcat(GLM_name,['_[rWLS]']); end
    if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 0
        aCompCor_fname = ['[aCompCor_' num2str(options.aCompCor(1)) 'WM_' num2str(options.aCompCor(2)) 'CSF]'];
    elseif (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 1
        aCompCor_fname = ['[aCompCor_' num2str(options.aCompCor(1)) 'WM_' num2str(options.aCompCor(2)) 'CSF_Ort_' options.motion '_HPF]'];
    elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 0
        aCompCor_fname = '[aCompCor50]';
    elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 1
        aCompCor_fname = ['[aCompCor50_Ort_' options.motion '_HPF]'];
    end
    if options.spikereg == 1; GLM_name = strcat(GLM_name,['_[SpikeReg_' num2str(options.spikeregFDthr) 'mm]']); end
    if ~strcmp(options.GSR,'none'); GLM_name = strcat(GLM_name,['_[' options.GSR ']']); end
    if ~strcmp(options.WM_CSF,'none'); GLM_name = strcat(GLM_name,['_[' options.WM_CSF ']']); end
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
jSub = 0;
for iSub = 1:length(SPM_paths)
    if exist(fullfile(output_paths{iSub},'SPM.mat'),'file')
        SPMnew = load(fullfile(output_paths{iSub},'SPM.mat'));
        if ~isfield(SPMnew.SPM,'Vbeta')
            rmdir(output_paths{iSub},'s');
        end
        clear SPMnew
    end
    if ~exist(fullfile(output_paths{iSub}),'dir')
        mkdir(fullfile(output_paths{iSub}));
    end

    % Specify GLM with noise regressors
    %----------------------------------------------------------------------
    if ~exist(fullfile(output_paths{iSub},'SPM.mat'),'file')

        old_GLM_subfolder = fileparts(SPM_paths{iSub});

        % Load HMP
        %----------------------------------------------------------------------
        if strcmp(options.motion,'12HMP')
            HMP = load(fullfile(old_GLM_subfolder,'TMFC_denoise','12HMP.mat')).HMP12;
        elseif strcmp(options.motion,'24HMP')
            HMP = load(fullfile(old_GLM_subfolder,'TMFC_denoise','24HMP.mat')).HMP24;
        end
    
        % Load spike regressors
        %----------------------------------------------------------------------
        if options.spikereg == 1
            SpikeReg = load(fullfile(old_GLM_subfolder,'TMFC_denoise',['SpikeReg_[FDthr_' num2str(options.spikeregFDthr) 'mm].mat'])).SpikeReg;
        end
    
        % Load GSR
        %----------------------------------------------------------------------
        if strcmp(options.GSR,'GSR')
            GSR = load(fullfile(masks.glm_paths{iSub},[options.GSR '.mat'])).GSR;
        elseif strcmp(options.GSR,'2GSR')
            GSR = load(fullfile(masks.glm_paths{iSub},[options.GSR '.mat'])).GSR2;
        elseif strcmp(options.GSR,'4GSR')
            GSR = load(fullfile(masks.glm_paths{iSub},[options.GSR '.mat'])).GSR4;    
        end
    
        % Load Phys
        %----------------------------------------------------------------------
        if strcmp(options.WM_CSF,'2Phys')
            Phys = load(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat'])).Phys2;
        elseif strcmp(options.WM_CSF,'4Phys')
            Phys = load(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat'])).Phys4;
        elseif strcmp(options.WM_CSF,'8Phys')
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
        % (if spm_fmri_concatenate.m sript was used)
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

        %---------------------------------------------Loop throuph sessions
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
            if ~strcmp(options.motion,'6HMP')
                Conf = [Conf HMP(jSess).Sess(:,7:end)]; % Start with 7th HMP, as 6HMP is already added
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
            if ~strcmp(options.GSR,'none')
                Conf = [Conf GSR(jSess).Sess];
                C = cell(size(GSR(jSess).Sess,2),1); C(:) = {'GSR'};
                ConfName = [ConfName; C]; clear C
            end
            % Phys
            if ~strcmp(options.WM_CSF,'none')
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
        if strcmp(SPM.xBF.name,'hrf')
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        elseif strcmp(SPM.xBF.name,'hrf (with time derivative)')
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
        elseif strcmp(SPM.xBF.name,'hrf (with time and dispersion derivatives)')
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
        elseif strcmp(SPM.xBF.name,'Fourier set')
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.length = SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier.order = SPM.xBF.order;
        elseif strcmp(SPM.xBF.name,'Fourier set (Hanning)')
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.length = SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.fourier_han.order = SPM.xBF.order;
        elseif strcmp(SPM.xBF.name,'Gamma functions')
            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.length = SPM.xBF.length;
            matlabbatch{1}.spm.stats.fmri_spec.bases.gamma.order = SPM.xBF.order;
        elseif strcmp(SPM.xBF.name,'Finite Impulse Response')
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
    
        if strcmp(SPM.xVi.form,'i.i.d') || strcmp(SPM.xVi.form,'none')
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'None';
        elseif strcmp(SPM.xVi.form,'fast') || strcmp(SPM.xVi.form,'FAST')
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
        spm('defaults','fmri');
        spm_get_defaults('cmdline',true);
        spm_jobman('run',batch{1});
        % Concatenated sessions
        if SPM_concat == 1
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
        try % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;           
            w = waitbar(0,'Please wait...','Name','Model estimation','Tag','tmfc_waitbar');
            afterEach(D, @tmfc_parfor_waitbar);    
            tmfc_parfor_waitbar(w,jSub,1);
        end
        if options.parallel == 0
            M = 0;
        else
            M = maxNumCompThreads('automatic');
        end
        parfor (iSub = 1:jSub,M) % Run matlabbatches in parallel mode
            spm('defaults','fmri');
            spm_get_defaults('cmdline',true);
            spm_jobman('run',batch{iSub});
            % Concatenated sessions
            if SPM_concat == 1
                spm_fmri_concatenate(fullfile(batch{iSub}{1}.spm.stats.fmri_spec.dir,'SPM.mat'),concat(iSub).scans);
            end
            % Check for rWLS
            if options.rWLS == 0
                spm_jobman('run',batch_2{iSub});
            else
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
            try send(D,[]); end % Update waitbar
        end
        try delete(w); end % Close waitbar
    end
end
end

%==========================================================================
% Waitbar for parallel mode
%--------------------------------------------------------------------------
function tmfc_parfor_waitbar(waitbarHandle,iterations,firstsub)
    persistent w nSub start_sub start_time count_sub 
    if nargin == 3
        w = waitbarHandle;
        nSub = iterations;
        start_sub = firstsub - 1;
        start_time = tic;
        count_sub = 1;
    else
        if isvalid(w)         
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            iSub = start_sub + count_sub;
            time_remaining = (nSub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            count_sub = count_sub + 1;
        end
    end
end

