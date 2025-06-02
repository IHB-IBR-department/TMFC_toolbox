function tmfc_spikereg(SPM_paths,options)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Creates spike regressors based on a user-defined FD threshold. For each
% flagged time point, a unit impulse function is included in general linear
% model, which had the value of 1 at that time point and 0 elsewhere
% (Lemieux et al., 2007; Satterthwaite et al., 2013).
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
