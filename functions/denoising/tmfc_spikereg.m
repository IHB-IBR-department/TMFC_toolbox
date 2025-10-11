function tmfc_spikereg(SPM_paths,options)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Creates spike regressors based on a user-defined FD threshold. For each
% flagged time point, a unit impulse function is included in the general
% linear model, with a value of 1 at that time point and 0 elsewhere
% (Lemieux et al., 2007; Satterthwaite et al., 2012).
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if ~isfield(options,'spikeregFDthr') || ~isscalar(options.spikeregFDthr) || ~isnumeric(options.spikeregFDthr) || options.spikeregFDthr < 0
    error('FD threshold (FDthr) for spike regression must be a non-negative scalar.');
end

for iSub = 1:length(SPM_paths)
    GLM_subfolder = fileparts(SPM_paths{iSub});
    outdir = fullfile(GLM_subfolder,'TMFC_denoise');
    outname = sprintf('SpikeReg_[FDthr_%.2fmm].mat', options.spikeregFDthr);
    if ~exist(fullfile(outdir,outname),'file')
        FD = load(fullfile(outdir,'FD.mat')).FramewiseDisplacement.Sess;
        for jSess = 1:length(FD)
            SpikeIdx = find(FD(jSess).FD_ts > options.spikeregFDthr);
            SpikeReg(jSess).Sess = zeros(length(FD(jSess).FD_ts),length(SpikeIdx));
            for kScan = 1:length(SpikeIdx)
                SpikeReg(jSess).Sess(SpikeIdx(kScan),kScan) = 1;
            end
            clear SpikeIdx
        end
        save(fullfile(outdir,outname),'SpikeReg');
        clear FD SpikeReg 
    end
    clear GLM_subfolder 
end
end
