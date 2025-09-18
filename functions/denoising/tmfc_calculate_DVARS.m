function [preDVARS,postDVARS] = tmfc_calculate_DVARS(FD,SPM_paths,options,masks,output_paths)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
%
% Calculates Derivative of root mean square VARiance over voxelS (DVARS).
% DVARS is the root mean square (RMS) of the temporal derivative of
% BOLD time series within the GM mask (Muschelli et al., 2014).
%
% Computes FD-DVARS correlations. 
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
        end
        DVARS.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
        DVARS.Max_FD_DVARS_corr = max(FD_DVARS_corr);
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
            end
            DVARS.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
            DVARS.Max_FD_DVARS_corr = max(FD_DVARS_corr);
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
            end
            DVARS.Mean_FD_DVARS_corr = mean(FD_DVARS_corr);
            DVARS.Max_FD_DVARS_corr = max(FD_DVARS_corr);
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
