function tmfc_spikereg(SPM_paths,options)
for iSub = 1:length(SPM_paths)
    GLM_subfolder = fileparts(SPM_paths{iSub});
    if ~exist(fullfile(GLM_subfolder,'TMFC_denoise',['SpikeReg_[FDthr_' num2str(options.spikeregFDthr) 'mm].mat']),'file')
        FD = load(fullfile(GLM_subfolder,'TMFC_denoise','FD.mat')).FramewiseDisplacement.Sess;
        for jSess = 1:length(FD)
            SpikeIdx = find(FD(jSess).FD_ts > options.spikeregFDthr);
            SpikeReg(jSess).Sess = zeros(length(FD(jSess).FD_ts),length(SpikeIdx));
            for kScan = 1:length(SpikeIdx)
                SpikeReg(jSess).Sess(SpikeIdx(kScan),kScan) = 1;
            end
            clear SpikeIdx
        end
        save(fullfile(GLM_subfolder,'TMFC_denoise',['SpikeReg_[FDthr_' num2str(options.spikeregFDthr) 'mm].mat']),'SpikeReg');
        clear FD SpikeReg 
    end
    clear GLM_subfolder 
end
end
