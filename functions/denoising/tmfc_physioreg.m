function tmfc_physioreg(SPM_paths,subject_paths,func_paths,masks,options)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
%
% (1) Creates aCompCor regressors (Behzadi et al., 2007). Calculates either 
%     a fixed number of principal components (PCs) or a variable number of
%     PCs explaining 50% of variance, separately within eroded WM and
%     CSF masks (Muschelli et al., 2014).   
% 
% (2) Creates WM/CSF regressors (Fox et al., 2005). Calculates average
%     BOLD signals separately for eroded WM and CSF masks. Optionally
%     adds temporal derivatives and quadratic terms (Parkes et al., 2017).
%
% (3) Creates GSR regressors (Fox et al., 2005, 2009). Computes the average
%     whole-brain BOLD signal. Optionally adds temporal derivatives and
%     quadratic terms (Parkes et al., 2017).
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru


% Extract signals from unsmoothed functional images
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

w = waitbar(0,'Please wait...','Name','Calculating physiological regressors');
for iSub = 1:length(SPM_paths)
    doGSR = false; doPhys = false; doFixedaCompCor = false; do50aCompCor = false;
    doGSR = ~strcmpi(options.GSR,'none') && ~exist(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'file');
    doPhys = ~strcmpi(options.WM_CSF,'none') && ~exist(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'file');
    doFixedaCompCor = (options.aCompCor(1) > 1 || options.aCompCor(2) > 1) && ~exist(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'file');
    do50aCompCor = options.aCompCor(1) == 0.5 && ~exist(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'file');

    if ~(doGSR || doPhys || doFixedaCompCor || do50aCompCor)
        try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
        continue;
    end

    GLM_subfolder = fileparts(SPM_paths{iSub});
    [~, sub, ~] = fileparts(subject_paths{iSub});

    % Load whole-brain mask
    %----------------------------------------------------------------------
    if doGSR
        WB_mask = spm_read_vols(spm_vol(masks.WB{iSub})); WB_mask = WB_mask(:); WB_mask(WB_mask == 0) = NaN;
    end

    % Load WM and CSF masks
    %----------------------------------------------------------------------
    WM_CSF = 0; aCompCorFixed = 0; aCompCorVar50 = 0;
    if doFixedaCompCor
        WM_CSF = 1; aCompCorFixed = 1;
    elseif do50aCompCor
        WM_CSF = 1; aCompCorVar50 = 1;
    elseif doPhys
        WM_CSF = 1;
    end
    if WM_CSF == 1
        WM_mask = spm_read_vols(spm_vol(masks.WM{iSub})); WM_mask = WM_mask(:); WM_mask(WM_mask == 0) = NaN;
        CSF_mask = spm_read_vols(spm_vol(masks.CSF{iSub})); CSF_mask = CSF_mask(:); CSF_mask(CSF_mask == 0) = NaN; 
    end

    % Load HMP
    %----------------------------------------------------------------------
    if (options.aCompCor_ort == 1 && sum(options.aCompCor)~=0)
        if strcmpi(options.motion,'12HMP')
            load(fullfile(GLM_subfolder,'TMFC_denoise','12HMP.mat'),'HMP12');
        elseif strcmpi(options.motion,'24HMP')
            load(fullfile(GLM_subfolder,'TMFC_denoise','24HMP.mat'),'HMP24');
        end
    end

    % Extract signals
    %----------------------------------------------------------------------
    SPM = load(SPM_paths{iSub}).SPM;
    for jSess = 1:length(SPM.Sess)
        clear WB_data WM_data CSF_data
        nScan = numel(SPM.Sess(jSess).row);
        for kScan = 1:nScan
            data = spm_read_vols(spm_vol(func_paths(iSub).fname{SPM.Sess(jSess).row(kScan)})); data = data(:);
            try; WB_data(:,kScan) = data(~isnan(data.*WB_mask)); end  
            try; WM_data(:,kScan) = data(~isnan(data.*WM_mask)); CSF_data(:,kScan) = data(~isnan(data.*CSF_mask)); end
        end
        clear nScan

        % Remove voxels with zero variance 
        %------------------------------------------------------------------
        if exist('WB_mask','var'); WB_data(std(WB_data,0,2) == 0,:) = []; end
        if WM_CSF == 1
            if isempty(WM_data); error(['The eroded WM mask is empty. Create a more liberal WM mask. Check Subject No. ' num2str(iSub) ': ' sub]); end
            if isempty(CSF_data); error(['The eroded CSF mask is empty. Create a more liberal CSF mask. Check Subject No. ' num2str(iSub) ': ' sub]); end
            WM_data(std(WM_data,0,2) == 0,:) = []; 
            CSF_data(std(CSF_data,0,2) == 0,:) = []; 
            if isempty(WM_data); error(['The eroded WM mask is empty. Create a more liberal WM mask. Check Subject No. ' num2str(iSub) ': ' sub]); end
            if isempty(CSF_data); error(['The eroded CSF mask is empty. Create a more liberal CSF mask. Check Subject No. ' num2str(iSub) ': ' sub]); end
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
        if doPhys
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
            % Pre-orthogonalize data w.r.t. head motion and high-pass filter 
            if options.aCompCor_ort == 1
                if strcmpi(options.motion,'12HMP')
                    HMP = HMP12(jSess).Sess; 
                elseif strcmpi(options.motion,'24HMP')
                    HMP = HMP24(jSess).Sess;
                else
                    HMP = SPM.Sess(jSess).C.C(:,[options.translation_idx options.rotation_idx]);
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
            WM_data = tmfc_zscore(WM_data);
            CSF_data = tmfc_zscore(CSF_data);
        end

        % aCompCor (fixed number of PCs)
        %------------------------------------------------------------------
        if aCompCorFixed == 1
            nVox_WM = size(WM_data,2); nVox_CSF = size(CSF_data,2); nTime = size(WM_data,1);
            % Check WM mask
            if options.aCompCor(1) > nVox_WM
                error('The number of requested PCs exceeds the number of voxels in the eroded WM mask.\n%s',['Create a more liberal WM mask or select fewer PCs. Check Subject No.' num2str(iSub) ': ' sub]);
            end
            % Check CSF mask
            if options.aCompCor(2) > nVox_CSF
                error('The number of requested PCs exceeds the number of voxels in the eroded CSF mask.\n%s',['Create a more liberal CSF mask or select fewer PCs. Check Subject No.' num2str(iSub) ': ' sub]);
            end
            % Calculate WM PCs
            [U,S] = svd(WM_data,'econ');
            latent = diag(S).^2/(nTime-1);
            latent = 100*latent/sum(latent);
            aCompCor.Sess(jSess).WM_PCs = U(:,1:options.aCompCor(1));
            var_expl = cumsum(latent(1:options.aCompCor(1)));
            aCompCor.Sess(jSess).WM_variance_explained = var_expl(end);
            clear U S latent var_expl
            % Calculate CSF PCs
            [U,S] = svd(CSF_data,'econ');
            latent = diag(S).^2/(nTime-1);
            latent = 100*latent/sum(latent);
            aCompCor.Sess(jSess).CSF_PCs = U(:,1:options.aCompCor(2));
            var_expl = cumsum(latent(1:options.aCompCor(2)));
            aCompCor.Sess(jSess).CSF_variance_explained = var_expl(end);
            clear U S latent var_expl 
        end

        % aCompCor (50% variance explained)
        %------------------------------------------------------------------
        if aCompCorVar50 == 1
            nVox_WM = size(WM_data,2); nVox_CSF = size(CSF_data,2); nTime = size(WM_data,1);
            % Calculate WM PCs
            [U,S] = svd(WM_data,'econ');
            latent = diag(S).^2/(nTime-1);
            latent = 100*latent/sum(latent);
            idx50 = find(cumsum(latent)>=50, 1, 'first');
            aCompCor50.Sess(jSess).WM_PCs = U(:,1:idx50);
            var_expl = cumsum(latent(1:idx50));
            aCompCor50.Sess(jSess).WM_variance_explained = var_expl(end);
            clear U S latent var_expl idx50
            % Calculate CSF PCs
            [U,S] = svd(CSF_data,'econ');
            latent = diag(S).^2/(nTime-1);
            latent = 100*latent/sum(latent);
            idx50 = find(cumsum(latent)>=50, 1, 'first');
            aCompCor50.Sess(jSess).CSF_PCs = U(:,1:idx50);
            var_expl = cumsum(latent(1:idx50));
            aCompCor50.Sess(jSess).CSF_variance_explained = var_expl(end);
            clear U S latent var_expl idx50
        end
        clear nVox_WM nVox_CSF
        clear WB_data WM_data CSF_data 
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

    clear GLM_subfolder sub SPM WB_mask WM_mask CSF_mask
    
    % Save *.mat files
    %----------------------------------------------------------------------
    % Save GSR 
    if doGSR
        if strcmpi(options.GSR,'GSR'); save(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'GSR');
        elseif strcmpi(options.GSR,'2GSR'); save(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'GSR2');
        elseif strcmpi(options.GSR,'4GSR');save(fullfile(masks.glm_paths{iSub},[options.GSR '.mat']),'GSR4'); end
    end
    % Save Phys
    if doPhys
        if strcmpi(options.WM_CSF,'2Phys'); save(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'Phys2');
        elseif strcmpi(options.WM_CSF,'4Phys'); save(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'Phys4'); 
        elseif strcmpi(options.WM_CSF,'8Phys'); save(fullfile(masks.glm_paths{iSub},[options.WM_CSF '.mat']),'Phys8'); end    
    end
    % Save aCompCor
    if doFixedaCompCor
        save(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'aCompCor');
    elseif do50aCompCor
        save(fullfile(masks.glm_paths{iSub},[aCompCor_fname '.mat']),'aCompCor50');
    end
    clear HMP12 HMP24 GSR GSR2 GSR4 Phys2 Phys4 Phys8 aCompCor aCompCor50 
    try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
end
try, if exist('w','var') && ishghandle(w), close(w); end, end

end
