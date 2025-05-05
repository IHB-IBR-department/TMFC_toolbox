function tmfc_physioreg(SPM_paths,subject_paths,funct_paths,masks,options)

% Extract signals from unsmoothed functional images
%--------------------------------------------------------------------------
if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 0
    aCompCor_fname = ['aCompCor_[' num2str(options.aCompCor(1)) 'WM]_[' num2str(options.aCompCor(2)) 'CSF]'];
elseif (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && options.aCompCor_ort == 1
    aCompCor_fname = ['aCompCor_[' num2str(options.aCompCor(1)) 'WM]_[' num2str(options.aCompCor(2)) 'CSF]_[ort_wrt_' options.motion '_and_HPF]'];
elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 0
    aCompCor_fname = 'aCompCor50';
elseif options.aCompCor(1) == 0.5 && options.aCompCor_ort == 1
    aCompCor_fname = ['aCompCor50_[ort_wrt_' options.motion '_and_HPF]'];
end

DVARS2_fname = ['DVARS_[' options.motion ']']; 
if options.spikereg == 1; DVARS2_fname = strcat(DVARS2_fname,['_[SpikeReg_' num2str(options.spikeregFDthr) 'mm]']); end
if ~strcmp(options.GSR,'none'); DVARS2_fname = strcat(DVARS2_fname,['_[' options.GSR ']']); end
if ~strcmp(options.WM_CSF,'none'); DVARS2_fname = strcat(DVARS2_fname,['_[' options.WM_CSF ']']); end
if sum(options.aCompCor)~=0; DVARS2_fname = strcat(DVARS2_fname,['_' aCompCor_fname]); end
DVARS2_fname = strcat(DVARS2_fname,'.mat');

w = waitbar(0,'Please wait...','Name','Calculating physiological regressors');
for iSub = 1:length(SPM_paths)
    GLM_subfolder = fileparts(SPM_paths{iSub});
    [~, sub, ~] = fileparts(subject_paths{iSub});

    % Load whole-brain mask
    %----------------------------------------------------------------------
    if ~strcmp(options.GSR,'none') && ~exist(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'file')
        WB_mask = spm_read_vols(spm_vol(masks.WB{iSub})); WB_mask = WB_mask(:); WB_mask(WB_mask == 0) = NaN;
    end

    % Load GM mask
    %----------------------------------------------------------------------
    if options.DVARS == 1 && (~exist(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat'),'file') || ~exist(fullfile(masks.glm_paths{iSub},DVARS2_fname),'file'))
        GM_mask = spm_read_vols(spm_vol(masks.GM{iSub})); GM_mask = GM_mask(:); GM_mask(GM_mask == 0) = NaN;
    end

    % Load WM and CSF masks
    %----------------------------------------------------------------------
    WM_CSF = 0; aCompCorFixed = 0; aCompCorVar50 = 0;
    if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && ~exist(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'file')
        WM_CSF = 1; aCompCorFixed = 1;
    elseif options.aCompCor(1) == 0.5 && ~exist(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'file')
        WM_CSF = 1; aCompCorVar50 = 1;
    elseif ~strcmp(options.WM_CSF,'none') && ~exist(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'file')
        WM_CSF = 1;
    end
    if WM_CSF == 1
        WM_mask = spm_read_vols(spm_vol(masks.WM{iSub})); WM_mask = WM_mask(:); WM_mask(WM_mask == 0) = NaN;
        CSF_mask = spm_read_vols(spm_vol(masks.CSF{iSub})); CSF_mask = CSF_mask(:); CSF_mask(CSF_mask == 0) = NaN; 
    end

    % Load FD
    %----------------------------------------------------------------------
    if options.DVARS == 1
        FD = load(fullfile(GLM_subfolder,'TMFC_denoise','FD.mat')).FramewiseDisplacement.Sess;
    end

    % Load HMP
    %----------------------------------------------------------------------
    if options.DVARS == 1 || (options.aCompCor_ort == 1 && sum(options.aCompCor)~=0)
        if strcmp(options.motion,'12HMP')
            load(fullfile(GLM_subfolder,'TMFC_denoise','12HMP.mat'),'HMP12');
        elseif strcmp(options.motion,'24HMP')
            load(fullfile(GLM_subfolder,'TMFC_denoise','24HMP.mat'),'HMP24');
        end
    end

    % Load spike regressors
    %----------------------------------------------------------------------
    if options.DVARS == 1 && options.spikereg == 1
        load(fullfile(GLM_subfolder,'TMFC_denoise',['SpikeReg_[FDthr_' num2str(options.spikeregFDthr) 'mm].mat']),'SpikeReg');
    end

    % Load GSR
    %----------------------------------------------------------------------
    if options.DVARS == 1 && (~strcmp(options.GSR,'none') && exist(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'file'))
        load(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']));
    end

    % Load Phys
    %----------------------------------------------------------------------
    if options.DVARS == 1 && (~strcmp(options.WM_CSF,'none') && exist(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'file'))
        load(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']));
    end

    % Load aCompCor
    %----------------------------------------------------------------------
    if options.DVARS == 1 && (sum(options.aCompCor)~=0 && exist(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'file'))
        load(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']));
    end

    % Extract signals
    %----------------------------------------------------------------------
    SPM = load(SPM_paths{iSub}).SPM;
    for jSess = 1:length(SPM.Sess)
        if size(SPM.nscan,2) == size(SPM.Sess,2)
            nScan = SPM.nscan(jSess);
        else
            nScan = length(SPM.xY.VY);
        end
        for kScan = 1:nScan
            data = spm_read_vols(spm_vol(funct_paths(iSub).fname{SPM.Sess(jSess).row(kScan)})); data = data(:);
            try; WB_data(:,kScan) = data(~isnan(data.*WB_mask)); end
            try; GM_data(:,kScan) = data(~isnan(data.*GM_mask)); end   
            try; WM_data(:,kScan) = data(~isnan(data.*WM_mask)); CSF_data(:,kScan) = data(~isnan(data.*CSF_mask)); end
        end
        clear nScan

        % Remove voxels with zero variance 
        %------------------------------------------------------------------
        if exist('WB_mask','var'); WB_data(std(WB_data,0,2) == 0,:) = []; end
        if exist('GM_mask','var'); GM_data(std(GM_data,0,2) == 0,:) = []; end
        if WM_CSF == 1
            if isempty(WM_data); error(['The eroded WM mask is empty. Create a more liberal WM mask. Check Subject No.' num2str(iSub) ': ' sub]); end
            if isempty(CSF_data); error(['The eroded CSF mask is empty. Create a more liberal CSF mask. Check Subject No.' num2str(iSub) ': ' sub]); end
            WM_data(std(WM_data,0,2) == 0,:) = []; 
            CSF_data(std(CSF_data,0,2) == 0,:) = []; 
            if isempty(WM_data); error(['The eroded WM mask is empty. Create a more liberal WM mask. Check Subject No.' num2str(iSub) ': ' sub]); end
            if isempty(CSF_data); error(['The eroded CSF mask is empty. Create a more liberal CSF mask. Check Subject No.' num2str(iSub) ': ' sub]); end
        end

        % GSR
        %------------------------------------------------------------------
        if exist('WB_data','var')
            WB_mean = spm_detrend(mean(WB_data))';
            WB_mean_diff = [0; diff(WB_mean)];
            GSR(jSess).Sess = WB_mean;
            GSR2(jSess).Sess = [WB_mean WB_mean_diff];
            GSR4(jSess).Sess = [WB_mean WB_mean_diff WB_mean.^2 WB_mean_diff.^2];
        end

        % Phys
        %------------------------------------------------------------------
        if ~strcmp(options.WM_CSF,'none') && ~exist(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'file')
            Phys = [spm_detrend(mean(WM_data))' spm_detrend(mean(CSF_data))'];
            Phys_diff = [0 0; diff(Phys)];
            Phys2(jSess).Sess = Phys;
            Phys4(jSess).Sess = [Phys Phys_diff];
            Phys8(jSess).Sess = [Phys Phys_diff Phys.^2 Phys_diff.^2];
        end

        % Prepare WM/CSF data for SVD
        %------------------------------------------------------------------
        if aCompCorFixed == 1 || aCompCorVar50 == 1
            % First dimension - time
            WM_data = WM_data'; CSF_data = CSF_data';
            % Pre-orthogonalise data w.r.t. head motion and high pass filter 
            if options.aCompCor_ort == 1
                if strcmp(options.motion,'12HMP')
                    HMP = HMP12(jSess).Sess; 
                elseif strcmp(options.motion,'24HMP')
                    HMP = HMP24(jSess).Sess;
                else
                    HMP = SPM.Sess(jSess).C.C(:,1:6);
                end
                HMP = [spm_detrend(HMP) ones(size(HMP,1),1)];
                if size(SPM.nscan,2) == size(SPM.Sess,2)
                    W = SPM.xX.W(SPM.xX.K(jSess).row,SPM.xX.K(jSess).row);
                else % concatenated GLM
                    W = SPM.xX.W;
                end
                % 'Filtered and whitened' confounds
                if size(SPM.nscan,2) == size(SPM.Sess,2)
                    xKXs = spm_sp('Set',spm_filter(SPM.xX.K(jSess),W*HMP)); 
                else % concatenated GLM
                    xKXs = spm_sp('Set',spm_filter(SPM.xX.K,W*HMP)); 
                end
                xKXs.X = full(xKXs.X);
                % 'Filtered and whitened' data
                if size(SPM.nscan,2) == size(SPM.Sess,2)
                    WM_KWY = spm_filter(SPM.xX.K(jSess),W*WM_data);
                    CSF_KWY = spm_filter(SPM.xX.K(jSess),W*CSF_data);
                else % concatenated GLM
                    WM_KWY = spm_filter(SPM.xX.K,W*WM_data);
                    CSF_KWY = spm_filter(SPM.xX.K,W*CSF_data);
                end
                % Residuals
                WM_data = spm_sp('r',xKXs,WM_KWY); 
                CSF_data = spm_sp('r',xKXs,CSF_KWY);
                clear HMP W xKXs WM_KWY CSF_KWY
            end
            % De-mean and normalize
            WM_data = zscore(WM_data);
            CSF_data = zscore(CSF_data);
        end

        % aCompCor (fixed number of PCs)
        %------------------------------------------------------------------
        if aCompCorFixed == 1
            nVox_WM = size(WM_data,2); nVox_CSF = size(CSF_data,2);
            % Check WM mask
            if options.aCompCor(1) > nVox_WM
                error('The number of desired PCs exceeds the number of voxels in the eroded WM mask.\n%s',['Create a more liberal WM mask or select fewer PCs. Check Subject No.' num2str(iSub) ': ' sub]);
            end
            % Check CSF mask
            if options.aCompCor(2) > nVox_CSF
                error('The number of desired PCs exceeds the number of voxels in the eroded CSF mask.\n%s',['Create a more liberal CSF mask or select fewer PCs. Check Subject No.' num2str(iSub) ': ' sub]);
            end
            % Calculate WM PCs
            [U,S] = svd(WM_data);
            latent = diag(S).^2/(nVox_WM-1);
            latent = 100*latent/sum(latent);
            aCompCor.Sess(jSess).WM_PCs = U(:,1:options.aCompCor(1));
            var_expl = cumsum(latent(1:options.aCompCor(1)));
            aCompCor.Sess(jSess).WM_variance_explained = var_expl(end);
            clear U S latent var_expl
            % Calculate CSF PCs
            [U,S] = svd(CSF_data);
            latent = diag(S).^2/(nVox_CSF-1);
            latent = 100*latent/sum(latent);
            aCompCor.Sess(jSess).CSF_PCs = U(:,1:options.aCompCor(2));
            var_expl = cumsum(latent(1:options.aCompCor(2)));
            aCompCor.Sess(jSess).CSF_variance_explained = var_expl(end);
            clear U S latent var_expl 
        end

        % aCompCor 50%
        %------------------------------------------------------------------
        if aCompCorVar50 == 1
            nVox_WM = size(WM_data,2); nVox_CSF = size(CSF_data,2);
            % Calculate WM PCs
            [U,S] = svd(WM_data);
            latent = diag(S).^2/(nVox_WM-1);
            latent = 100*latent/sum(latent);
            idx50 = find(cumsum(latent)>50);
            aCompCor50.Sess(jSess).WM_PCs = U(:,1:idx50(1));
            var_expl = cumsum(latent(1:idx50(1)));
            aCompCor50.Sess(jSess).WM_variance_explained = var_expl(end);
            clear U S latent var_expl idx50
            % Calculate CSF PCs
            [U,S] = svd(CSF_data);
            latent = diag(S).^2/(nVox_CSF-1);
            latent = 100*latent/sum(latent);
            idx50 = find(cumsum(latent)>50);
            aCompCor50.Sess(jSess).CSF_PCs = U(:,1:idx50(1));
            var_expl = cumsum(latent(1:idx50(1)));
            aCompCor50.Sess(jSess).CSF_variance_explained = var_expl(end);
            clear U S latent var_expl idx50
        end

        % DVARS (before noise regression)
        %------------------------------------------------------------------
        if exist('GM_data','var')
            DVARS.Sess(jSess).DVARS_ts = rms([zeros(1,size(GM_data',2)); diff(GM_data')],2); DVARS.Sess(jSess).DVARS_ts(1) = NaN;
            DVARS.Sess(jSess).FD_DVARS_corr = corr(DVARS.Sess(jSess).DVARS_ts(2:end),FD(jSess).FD_ts(2:end));
        end

        % DVARS (after noise regression)
        %------------------------------------------------------------------
        if options.DVARS == 1 && ~exist(fullfile(masks.glm_paths{iSub},DVARS2_fname),'file')
            % HMP
            if strcmp(options.motion,'12HMP')
                Conf = HMP12(jSess).Sess; 
            elseif strcmp(options.motion,'24HMP')
                Conf = HMP24(jSess).Sess;
            else
                Conf = SPM.Sess(jSess).C.C(:,1:6);
            end
            % SpikeReg
            if options.spikereg == 1
                Conf = [Conf SpikeReg(jSess).Sess];
            end
            % GSR
            if strcmp(options.GSR,'GSR')
                Conf = [Conf GSR(jSess).Sess];
            elseif strcmp(options.GSR,'2GSR')
                Conf = [Conf GSR2(jSess).Sess];
            elseif strcmp(options.GSR,'4GSR')
                Conf = [Conf GSR4(jSess).Sess];
            end
            % Phys
            if strcmp(options.WM_CSF,'2Phys')
                Conf = [Conf Phys2(jSess).Sess];
            elseif strcmp(options.WM_CSF,'4Phys')
                Conf = [Conf Phys4(jSess).Sess];
            elseif strcmp(options.WM_CSF,'8Phys')
                Conf = [Conf Phys8(jSess).Sess];
            end
            % aCompCor
            if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) 
                Conf = [Conf aCompCor.Sess(jSess).WM_PCs aCompCor.Sess(jSess).CSF_PCs];
            elseif options.aCompCor(1) == 0.5 
                Conf = [Conf aCompCor50.Sess(jSess).WM_PCs aCompCor50.Sess(jSess).CSF_PCs];
            end
            % First dimension - time
            GM_data = GM_data';
            % Noise regression
            if size(SPM.nscan,2) == size(SPM.Sess,2)
                W = SPM.xX.W(SPM.xX.K(jSess).row,SPM.xX.K(jSess).row);
            else % concatenated GLM
                W = SPM.xX.W;
            end
            % 'Filtered and whitened' confounds
            Conf = [Conf ones(size(Conf,1),1)];
            if size(SPM.nscan,2) == size(SPM.Sess,2)
                xKXs = spm_sp('Set',spm_filter(SPM.xX.K(jSess),W*Conf)); 
            else % concatenated GLM
                xKXs = spm_sp('Set',spm_filter(SPM.xX.K,W*Conf)); 
            end
            xKXs.X = full(xKXs.X);
            % 'Filtered and whitened' data
            if size(SPM.nscan,2) == size(SPM.Sess,2)
                GM_KWY = spm_filter(SPM.xX.K(jSess),W*GM_data);
            else % concatenated GLM
                GM_KWY = spm_filter(SPM.xX.K,W*GM_data);
            end
            % Residuals
            GM_data = spm_sp('r',xKXs,GM_KWY); 
            clear W xKXs GM_KWY Conf
            % DVARS
            DVARS2.Sess(jSess).DVARS_ts = rms([zeros(1,size(GM_data,2)); diff(GM_data)],2); DVARS2.Sess(jSess).DVARS_ts(1) = NaN;
            DVARS2.Sess(jSess).FD_DVARS_corr = corr(DVARS2.Sess(jSess).DVARS_ts(2:end),FD(jSess).FD_ts(2:end));
        end

        clear nVox_WM nVox_CSF
        clear WB_data GM_data WM_data CSF_data 
        clear WB_mean WB_mean_diff Phys Phys_diff
    end %-----------------------------------------------end of session loop

    % aCompCor summary info
    %----------------------------------------------------------------------
    if aCompCorFixed == 1
        WM_var_expl = []; CSF_var_expl = [];
        for jSess = 1:length(SPM.Sess)
            WM_var_expl = [WM_var_expl aCompCor.Sess(jSess).WM_variance_explained];
            CSF_var_expl = [CSF_var_expl aCompCor.Sess(jSess).CSF_variance_explained];
        end
        aCompCor.WM_mean_variance_explained = mean(WM_var_expl);
        aCompCor.CSF_mean_variance_explained = mean(CSF_var_expl);
    end
    if aCompCorVar50 == 1
        WM_nPCs = 0; CSF_nPCs = 0; WM_nPCs_per_sess = []; CSF_nPCs_per_sess = [];
        for jSess = 1:length(SPM.Sess)
            WM_nPCs = WM_nPCs + size(aCompCor50.Sess(jSess).WM_PCs,2);
            CSF_nPCs = CSF_nPCs + size(aCompCor50.Sess(jSess).CSF_PCs,2);
            WM_nPCs_per_sess = [WM_nPCs_per_sess size(aCompCor50.Sess(jSess).WM_PCs,2)];
            CSF_nPCs_per_sess = [CSF_nPCs_per_sess size(aCompCor50.Sess(jSess).CSF_PCs,2)];
        end
        aCompCor50.WM_nPCs = WM_nPCs;
        aCompCor50.CSF_nPCs = CSF_nPCs;
        aCompCor50.Mean_WM_nPCs_per_sess = mean(WM_nPCs_per_sess);
        aCompCor50.Mean_CSF_nPCs_per_sess = mean(CSF_nPCs_per_sess);
    end

    % DVARS summary info (before noise regression)
    %----------------------------------------------------------------------
    if options.DVARS == 1 && ~exist(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat'),'file')
        FD_DVARS_corr = [];
        for jSess = 1:length(SPM.Sess)
            FD_DVARS_corr = [FD_DVARS_corr DVARS.Sess(jSess).FD_DVARS_corr];
        end
        DVARS.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
        DVARS.Max_FD_DVARS_corr = mean(FD_DVARS_corr);
        clear FD_DVARS_corr
    end

    % DVARS summary info (after noise regression)
    %----------------------------------------------------------------------
    if options.DVARS == 1 && ~exist(fullfile(masks.glm_paths{iSub},DVARS2_fname),'file')
        FD_DVARS_corr = [];
        for jSess = 1:length(SPM.Sess)
            FD_DVARS_corr = [FD_DVARS_corr DVARS2.Sess(jSess).FD_DVARS_corr];
        end
        DVARS2.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
        DVARS2.Max_FD_DVARS_corr = mean(FD_DVARS_corr);
        clear FD_DVARS_corr
    end

    clear GLM_subfolder sub SPM WB_mask GM_mask WM_mask CSF_mask
    
    % Save *.mat files
    %----------------------------------------------------------------------
    % Save GSR 
    if ~strcmp(options.GSR,'none') && ~exist(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'file')
        if strcmp(options.GSR,'GSR'); save(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'GSR');
        elseif strcmp(options.GSR,'2GSR'); save(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'GSR2');
        elseif strcmp(options.GSR,'4GSR');save(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'GSR4'); end
    end
    % Save DVARS (before noise regression)
    if options.DVARS == 1 && ~exist(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat'),'file')
        save(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat'),'DVARS')
    end
    % Save DVARS (after noise regression)
    if options.DVARS == 1 && ~exist(fullfile(masks.glm_paths{iSub},DVARS2_fname),'file')
        DVARS = DVARS2;
        save(fullfile(masks.glm_paths{iSub},DVARS2_fname),'DVARS')
    end
    % Save Phys
    if ~strcmp(options.WM_CSF,'none') && ~exist(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'file')
        if strcmp(options.WM_CSF,'2Phys'); save(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'Phys2');
        elseif strcmp(options.WM_CSF,'4Phys'); save(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'Phys4'); 
        elseif strcmp(options.WM_CSF,'8Phys'); save(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'Phys8'); end    
    end
    % Save aCompCor
    if (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && ~exist(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'file')
        save(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'aCompCor');
    elseif options.aCompCor(1) == 0.5 && ~exist(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'file')
        save(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'aCompCor50');
    end
    clear FD HMP12 HMP24 SpikeReg GSR GSR2 GSR4 Phys2 Phys4 Phys8 aCompCor aCompCor50 DVARS DVARS2
    try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
end
try; close(w); end % Close waitbar

% Plot DVARS
%--------------------------------------------------------------------------
if options.DVARS == 1
    for iSub = 1:length(SPM_paths)
        GLM_subfolder = fileparts(SPM_paths{iSub});
        preDVARS(iSub) = load(fullfile(masks.glm_paths{iSub},'DVARS_before_denoising.mat'));
        postDVARS(iSub) = load(fullfile(masks.glm_paths{iSub},DVARS2_fname));
        tmp = load(fullfile(GLM_subfolder,'TMFC_denoise','FD.mat')); 
        FD(iSub) = tmp.FramewiseDisplacement;
        clear GLM_subfolder tmp
    end

    tmfc_plot_DVARS(preDVARS,postDVARS,FD);
    pause(2);
end



