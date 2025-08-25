function [masks] = tmfc_create_masks(SPM_paths,struct_paths,funct_paths,options)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% Creates binary masks to extract BOLD signals from GM (for DVARS calculation),
% WM and CSF (for aCompCor and Phys regressors), and the whole brain (for GSR).
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

% Number of workers in a parallel pool
%--------------------------------------------------------------------------
if options.parallel == 0
    M = 0;
else
    M = maxNumCompThreads('automatic');
end

% Create mask subfolders
%--------------------------------------------------------------------------
tmfc_dir = fileparts(which('TMFC_denoise.m'));
for iSub = 1:length(SPM_paths)
    GLM_subfolder = fileparts(SPM_paths{iSub});
    segment_paths{iSub,1} = fullfile(GLM_subfolder,'TMFC_denoise','Segment');
    glm_paths{iSub,1} = fullfile(GLM_subfolder,'TMFC_denoise', ...
        ['[WM_' num2str(options.WMmask.prob) 'Prob_' num2str(options.WMmask.erode) ...
        'xErode]_[CSF_' num2str(options.CSFmask.prob) 'Prob_' num2str(options.CSFmask.erode) ...
        'xErode]_[GM_' num2str(options.GMmask.prob) 'Prob_' num2str(options.GMmask.dilate) 'xDilate]']);
    mask_paths{iSub,1} = fullfile(glm_paths{iSub,1},'Masks');
    if ~exist(segment_paths{iSub},'dir'); mkdir(segment_paths{iSub}); end
    if ~exist(mask_paths{iSub},'dir'); mkdir(mask_paths{iSub}); end
    clear GLM_subfolder
end
masks.glm_paths = glm_paths; 
masks.mask_paths = mask_paths;

% Copy structural T1 files
%--------------------------------------------------------------------------
for iSub = 1:length(SPM_paths)
    [~,struct_name,struct_ext] = fileparts(struct_paths{iSub});
    struct_file{iSub,1} = [struct_name,struct_ext];
    struct_copy_paths{iSub,1} = fullfile(segment_paths{iSub},[struct_name,struct_ext]);
    if ~exist(struct_copy_paths{iSub},'file')
        copyfile(struct_paths{iSub},struct_copy_paths{iSub});
    end
    clear struct_name struct_ext
end

% Segment T1 image
%--------------------------------------------------------------------------
spm_jobman('initcfg');

jSub = 0;
for iSub = 1:length(SPM_paths)
    forward_field{iSub,1} = fullfile(segment_paths{iSub},['y_',struct_file{iSub}]);
    bias_corrected{iSub,1} = fullfile(segment_paths{iSub},['m',struct_file{iSub}]);
    GM{iSub,1} = fullfile(segment_paths{iSub},['c1',struct_file{iSub}]);
    WM{iSub,1} = fullfile(segment_paths{iSub},['c2',struct_file{iSub}]);
    CSF{iSub,1} = fullfile(segment_paths{iSub},['c3',struct_file{iSub}]);
    if ~exist(bias_corrected{iSub},'file') || ~exist(GM{iSub},'file') || ~exist(WM{iSub},'file') || ~exist(CSF{iSub},'file')
        matlabbatch{1}.spm.spatial.preproc.channel.vols = struct_copy_paths(iSub);
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(fileparts(which('spm.m')),'tpm','TPM.nii,1')};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(fileparts(which('spm.m')),'tpm','TPM.nii,2')};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(fileparts(which('spm.m')),'tpm','TPM.nii,3')};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(fileparts(which('spm.m')),'tpm','TPM.nii,4')};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(fileparts(which('spm.m')),'tpm','TPM.nii,5')};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(fileparts(which('spm.m')),'tpm','TPM.nii,6')};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
        matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
        matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                                      NaN NaN NaN];
        jSub = jSub + 1;
        batch{jSub} = matlabbatch;
        clear matlabbatch
    end
end

if jSub == 1
    spm('defaults','fmri');
    spm_get_defaults('cmdline',true);
    spm_jobman('run',batch{1});
elseif jSub > 1
    try % Waitbar for MATLAB R2017a and higher
        D = parallel.pool.DataQueue;           
        w = waitbar(0,'Please wait...','Name','Segmentation','Tag','tmfc_waitbar');
        afterEach(D, @tmfc_parfor_waitbar);    
        tmfc_parfor_waitbar(w,jSub,1);
    end
    parfor (iSub = 1:jSub,M) % Run matlabbatches in parallel mode
        spm('defaults','fmri');
        spm_get_defaults('cmdline',true);
        spm_jobman('run',batch{iSub});
        try send(D,[]); end % Update waitbar
    end
    try delete(w); end % Close waitbar
end
clear batch

% Apply thresholds to tissue probability maps
%--------------------------------------------------------------------------
w = waitbar(0,'Please wait...','Name','Applying tissue probability thresholds');
for iSub = 1:length(SPM_paths)
    skull_stripped{iSub,1} = fullfile(mask_paths{iSub},'Skull_stripped_T1.nii');
    WB_mask{iSub,1} = fullfile(mask_paths{iSub},'Whole_brain_mask.nii');
    GM_mask{iSub,1} = fullfile(mask_paths{iSub},'GM_mask.nii');
    WM_mask{iSub,1} = fullfile(mask_paths{iSub},'WM_mask.nii');
    CSF_mask{iSub,1} = fullfile(mask_paths{iSub},'CSF_mask.nii');
    % Skull-stripped T1
    if ~exist(skull_stripped{iSub},'file')        
        input_images{1,1} = bias_corrected{iSub};
        input_images{2,1} = GM{iSub};
        input_images{3,1} = WM{iSub};
        input_images{4,1} = CSF{iSub};
        spm_imcalc(input_images,skull_stripped{iSub}, ...
            ['i1.*(((i2>0)+(i3>' num2str(options.WMmask.prob) ')+(i4>' num2str(options.CSFmask.prob) '))>0)'],{0,0,1,4});
        clear input_images
    end
    % Whole-brain binary mask
    if ~exist(WB_mask{iSub},'file')           
        input_images{1,1} = GM{iSub};
        input_images{2,1} = WM{iSub};
        input_images{3,1} = CSF{iSub};
        spm_imcalc(input_images,WB_mask{iSub}, ...
            ['((i1>0)+(i2>' num2str(options.WMmask.prob) ')+(i3>' num2str(options.CSFmask.prob) '))>0'],{0,0,0,2});
        clear input_images
    end
    % GM binary mask
    if ~exist(GM_mask{iSub},'file')           
        spm_imcalc(GM{iSub},GM_mask{iSub},['i1>' num2str(options.GMmask.prob)],{0,0,0,2});
    end
    % WM binary mask
    if ~exist(WM_mask{iSub},'file')           
        spm_imcalc(WM{iSub},WM_mask{iSub},['i1>' num2str(options.WMmask.prob)],{0,0,0,2});
    end
    % CSF binary mask 
    if ~exist(CSF_mask{iSub},'file')             
        spm_imcalc(CSF{iSub},CSF_mask{iSub},['i1>' num2str(options.CSFmask.prob)],{0,0,0,2});
    end
    try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
end
try; close(w); end % Close waitbar

%==========================================================================
% Creating WM and CSF masks
%--------------------------------------------------------------------------
if sum(options.aCompCor)~=0 || ~strcmp(options.WM_CSF,'none')

    % Erode WM masks
    %----------------------------------------------------------------------
    if options.WMmask.erode == 0
        WM_eroded{iSub,1} = WM_mask{iSub};
    else
        w = waitbar(0,'Please wait...','Name','Erosion of WM masks');
        for iSub = 1:length(SPM_paths)
            WM_eroded{iSub,1} = fullfile(mask_paths{iSub},'WM_mask_eroded.nii');
            if options.WMmask.erode > 0
                if ~exist(WM_eroded{iSub},'file')
                    V = spm_vol(WM_mask{iSub});
                    ima = spm_read_vols(V);
                    for iCycle = 1:options.WMmask.erode
                        ima = spm_erode(ima);
                    end
                    V.fname = WM_eroded{iSub,1};
                    spm_write_vol(V,ima);
                    clear V ima
                end
            end
            try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
        end
        try; close(w); end % Close waitbar
    end
    
    % Dilate liberal GM mask
    %----------------------------------------------------------------------
    if options.GMmask.dilate == 0
        GM_dilated{iSub,1} = GM_mask{iSub};
    else
        w = waitbar(0,'Please wait...','Name','Dilation of GM masks');
        for iSub = 1:length(SPM_paths)
            GM_dilated{iSub,1} = fullfile(mask_paths{iSub},'GM_mask_dilated.nii');
            if options.GMmask.dilate > 0
                if ~exist(GM_dilated{iSub},'file')
                    V = spm_vol(GM_mask{iSub});
                    ima = spm_read_vols(V);
                    for iCycle = 1:options.GMmask.dilate
                        ima = spm_dilate(ima);
                    end
                    V.fname = GM_dilated{iSub,1};
                    spm_write_vol(V,ima);
                    clear V ima
                end
            end
            try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
        end
        try; close(w); end % Close waitbar
    end
    
    % Subtract dilated GM mask from CSF mask
    %----------------------------------------------------------------------
    w = waitbar(0,'Please wait...','Name','Subtract dilated GM masks from CSF masks');
    for iSub = 1:length(SPM_paths)
        CSF_mask_GM_removed{iSub,1} = fullfile(mask_paths{iSub},'CSF_mask_GM_removed.nii');
        if ~exist(CSF_mask_GM_removed{iSub},'file')        
            input_images{1,1} = CSF_mask{iSub};
            input_images{2,1} = GM_dilated{iSub};
            spm_imcalc(input_images,CSF_mask_GM_removed{iSub},'(i1-i2)>0',{0,0,0,2});
            clear input_images
        end
        try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
    end
    try; close(w); end % Close waitbar
    
    % Erode CSF masks
    %----------------------------------------------------------------------
    if options.CSFmask.erode == 0
        CSF_eroded{iSub,1} = CSF_mask{iSub};
    else
        w = waitbar(0,'Please wait...','Name','Erosion of CSF masks');
        for iSub = 1:length(SPM_paths)
            CSF_eroded{iSub,1} = fullfile(mask_paths{iSub},'CSF_mask_eroded.nii');
            if options.CSFmask.erode > 0
                if ~exist(CSF_eroded{iSub},'file')
                    V = spm_vol(CSF_mask{iSub});
                    ima = spm_read_vols(V);
                    for iCycle = 1:options.CSFmask.erode
                        ima = spm_erode(ima);
                    end
                    V.fname = CSF_eroded{iSub,1};
                    spm_write_vol(V,ima);
                    clear V ima
                end
            end
            try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
        end
        try; close(w); end % Close waitbar
    end
    
    % Normalization and resampling to fMRI resolution
    %----------------------------------------------------------------------
    jSub = 0;
    for iSub = 1:length(SPM_paths)
        V = spm_vol(funct_paths(iSub).fname{1});
        [BB,vx] = spm_get_bbox(V);
        vx = abs(vx);
        coreg = 0;
        wWM_eroded{iSub,1} = fullfile(mask_paths{iSub},   'w_WM_mask_eroded.nii');
        wCSF_eroded{iSub,1} = fullfile(mask_paths{iSub},  'w_CSF_mask_eroded.nii');
        wskull_stripped{iSub} = fullfile(mask_paths{iSub},'w_Skull_stripped_T1.nii');
        if ~exist(wWM_eroded{iSub},'file') || ~exist(wCSF_eroded{iSub},'file') || ~exist(wskull_stripped{iSub},'file')
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = forward_field(iSub);
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1,1} = WM_eroded{iSub};
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{2,1} = CSF_eroded{iSub};
            if ~exist(wskull_stripped{iSub},'file')
                matlabbatch{1}.spm.spatial.normalise.write.subj.resample{3,1} = skull_stripped{iSub};
                coreg = 1;
            end
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BB;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_';
            matlabbatch{2}.spm.spatial.coreg.write.ref{1,1} = funct_paths(iSub).fname{1};
            matlabbatch{2}.spm.spatial.coreg.write.source{1,1} = wWM_eroded{iSub,1};
            matlabbatch{2}.spm.spatial.coreg.write.source{2,1} = wCSF_eroded{iSub,1};
            if coreg == 1
                matlabbatch{2}.spm.spatial.coreg.write.source{3,1} = wskull_stripped{iSub};
            end
            matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
            matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
            jSub = jSub + 1;
            batch{jSub} = matlabbatch;
            clear matlabbatch
        end
        wWM_eroded{iSub,1} = fullfile(mask_paths{iSub},   'rw_WM_mask_eroded.nii');
        wCSF_eroded{iSub,1} = fullfile(mask_paths{iSub},  'rw_CSF_mask_eroded.nii');
        wskull_stripped{iSub} = fullfile(mask_paths{iSub},'rw_Skull_stripped_T1.nii');
        clear V BB vx coreg
    end
    
    if jSub == 1
        spm('defaults','fmri');
        spm_get_defaults('cmdline',true);
        spm_jobman('run',batch{1});
    elseif jSub > 1
        try % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;           
            w = waitbar(0,'Please wait...','Name','Normalization: WM and CSF masks','Tag','tmfc_waitbar');
            afterEach(D, @tmfc_parfor_waitbar);    
            tmfc_parfor_waitbar(w,jSub,1);
        end
        parfor (iSub = 1:jSub,M) % Run matlabbatches in parallel mode
            spm('defaults','fmri');
            spm_get_defaults('cmdline',true);
            spm_jobman('run',batch{iSub});
            try; send(D,[]); end % Update waitbar
        end
        try; delete(w); end % Close waitbar
    end
    clear batch
    
    % Remove brainstem voxels from WM mask and applying implicit SPM mask
    %----------------------------------------------------------------------
    w = waitbar(0,'Please wait...','Name','Removing brainstem from WM masks');
    for iSub = 1:length(SPM_paths)
        wWM_eroded_final{iSub,1} = fullfile(mask_paths{iSub},'rw_WM_mask_eroded_no_brainstem.nii');
        if ~exist(wWM_eroded_final{iSub},'file')  
            SPM = load(SPM_paths{iSub}).SPM;
            input_images{1,1} = wWM_eroded{iSub};
            input_images{2,1} = fullfile(tmfc_dir,'functions','masks','Harvard_Oxford_Brainstem_mask.nii');
            input_images{3,1} = fullfile(SPM.swd,SPM.VM.fname); 
            spm_imcalc(input_images,wWM_eroded_final{iSub},'(i1-i2).*i3>0.5',{0,0,0,2});
            clear input_images SPM
        end
        try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
    end
    try; close(w); end % Close waitbar
    
    % Remove nonventricle voxels from CSF mask and applying implicit SPM mask
    %----------------------------------------------------------------------
    w = waitbar(0,'Please wait...','Name','Removing nonventricle voxels from CSF masks');
    for iSub = 1:length(SPM_paths)
        wCSF_eroded_final{iSub,1} = fullfile(mask_paths{iSub},'rw_CSF_mask_eroded_only_ventricles.nii');
        if ~exist(wCSF_eroded_final{iSub},'file')
            SPM = load(SPM_paths{iSub}).SPM;
            input_images{1,1} = wCSF_eroded{iSub};
            input_images{2,1} = fullfile(tmfc_dir,'functions','masks','ALVIN_mask.nii');
            input_images{3,1} = fullfile(SPM.swd,SPM.VM.fname); 
            spm_imcalc(input_images,wCSF_eroded_final{iSub},'(i1.*(i2>50).*i3)>0.5',{0,0,0,2});
            clear input_images SPM
        end
        try; waitbar(iSub/length(SPM_paths),w,['Subject No. ' num2str(iSub)]); end % Update waitbar
    end

    % Save paths
    %----------------------------------------------------------------------
    masks.WM = wWM_eroded_final;
    masks.CSF = wCSF_eroded_final;
    try; close(w); end % Close waitbar
end

%==========================================================================
% Creating GM masks
%--------------------------------------------------------------------------
if options.DVARS == 1

    % Normalization and resampling to fMRI resolution
    %----------------------------------------------------------------------
    jSub = 0;
    for iSub = 1:length(SPM_paths)
        V = spm_vol(funct_paths(iSub).fname{1});
        [BB,vx] = spm_get_bbox(V);
        vx = abs(vx);
        coreg = 0;
        wGM_mask{iSub,1} = fullfile(mask_paths{iSub},     'w_GM_mask.nii');
        wskull_stripped{iSub} = fullfile(mask_paths{iSub},'w_Skull_stripped_T1.nii');
        if ~exist(wGM_mask{iSub},'file') || ~exist(wskull_stripped{iSub},'file')
            SPM = load(SPM_paths{iSub}).SPM;
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = forward_field(iSub);
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1,1} = GM_mask{iSub};
            if ~exist(wskull_stripped{iSub},'file')
                matlabbatch{1}.spm.spatial.normalise.write.subj.resample{2,1} = skull_stripped{iSub};
                coreg = 1;
            end
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BB;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_';
            matlabbatch{2}.spm.spatial.coreg.write.ref{1,1} = funct_paths(iSub).fname{1};
            matlabbatch{2}.spm.spatial.coreg.write.source{1,1} = wGM_mask{iSub,1};
            if coreg == 1
                matlabbatch{2}.spm.spatial.coreg.write.source{2,1} = wskull_stripped{iSub};
            end
            matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
            matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
            wGM_mask{iSub,1} = fullfile(mask_paths{iSub},     'rw_GM_mask.nii');
            matlabbatch{3}.spm.util.imcalc.input{1,1} = wGM_mask{iSub};
            matlabbatch{3}.spm.util.imcalc.input{2,1} = fullfile(SPM.swd,SPM.VM.fname); 
            matlabbatch{3}.spm.util.imcalc.output = 'rw_GM_mask';
            matlabbatch{3}.spm.util.imcalc.outdir = mask_paths(iSub);
            matlabbatch{3}.spm.util.imcalc.expression = '(i1.*i2)>0.5';
            matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{3}.spm.util.imcalc.options.mask = 0;
            matlabbatch{3}.spm.util.imcalc.options.interp = 0;
            matlabbatch{3}.spm.util.imcalc.options.dtype = 2;
            jSub = jSub + 1;
            batch{jSub} = matlabbatch;
            clear matlabbatch SPM
        end
        wGM_mask{iSub,1} = fullfile(mask_paths{iSub},     'rw_GM_mask.nii');
        wskull_stripped{iSub} = fullfile(mask_paths{iSub},'rw_Skull_stripped_T1.nii');
        clear V BB vx coreg
    end

    if jSub == 1
        spm('defaults','fmri');
        spm_get_defaults('cmdline',true);
        spm_jobman('run',batch{1});
    elseif jSub > 1
        try % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;           
            w = waitbar(0,'Please wait...','Name','Creating GM masks','Tag','tmfc_waitbar');
            afterEach(D, @tmfc_parfor_waitbar);    
            tmfc_parfor_waitbar(w,jSub,1);
        end
        parfor (iSub = 1:jSub,M) % Run matlabbatches in parallel mode
            spm('defaults','fmri');
            spm_get_defaults('cmdline',true);
            spm_jobman('run',batch{iSub});
            try; send(D,[]); end % Update waitbar
        end
        try; delete(w); end % Close waitbar
    end
    clear batch

    % Save paths
    %----------------------------------------------------------------------
    masks.GM = wGM_mask;
end

%==========================================================================
% Creating whole-brain masks
%--------------------------------------------------------------------------
if ~strcmp(options.GSR,'none')

    % Normalization and resampling to fMRI resolution
    %----------------------------------------------------------------------
    jSub = 0;
    for iSub = 1:length(SPM_paths)
        V = spm_vol(funct_paths(iSub).fname{1});
        [BB,vx] = spm_get_bbox(V);
        vx = abs(vx);
        coreg = 0;
        wWB_mask{iSub,1} = fullfile(mask_paths{iSub},     'w_Whole_brain_mask.nii');
        wskull_stripped{iSub} = fullfile(mask_paths{iSub},'w_Skull_stripped_T1.nii');
        if ~exist(wWB_mask{iSub},'file') || ~exist(wskull_stripped{iSub},'file')
            SPM = load(SPM_paths{iSub}).SPM;
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = forward_field(iSub);
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample{1,1} = WB_mask{iSub};
            if ~exist(wskull_stripped{iSub},'file')
                matlabbatch{1}.spm.spatial.normalise.write.subj.resample{2,1} = skull_stripped{iSub};
                coreg = 1;
            end
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BB;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vx;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_';
            matlabbatch{2}.spm.spatial.coreg.write.ref{1,1} = funct_paths(iSub).fname{1};
            matlabbatch{2}.spm.spatial.coreg.write.source{1,1} = wWB_mask{iSub,1};
            if coreg == 1
                matlabbatch{2}.spm.spatial.coreg.write.source{2,1} = wskull_stripped{iSub};
            end
            matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
            matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
            wWB_mask{iSub,1} = fullfile(mask_paths{iSub},     'rw_Whole_brain_mask.nii');
            matlabbatch{3}.spm.util.imcalc.input{1,1} = wWB_mask{iSub};
            matlabbatch{3}.spm.util.imcalc.input{2,1} = fullfile(SPM.swd,SPM.VM.fname); 
            matlabbatch{3}.spm.util.imcalc.output = 'rw_Whole_brain_mask';
            matlabbatch{3}.spm.util.imcalc.outdir = mask_paths(iSub);
            matlabbatch{3}.spm.util.imcalc.expression = '(i1.*i2)>0.5';
            matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{3}.spm.util.imcalc.options.mask = 0;
            matlabbatch{3}.spm.util.imcalc.options.interp = 0;
            matlabbatch{3}.spm.util.imcalc.options.dtype = 2;
            jSub = jSub + 1;
            batch{jSub} = matlabbatch;
            clear matlabbatch SPM
        end
        wWB_mask{iSub,1} = fullfile(mask_paths{iSub},     'rw_Whole_brain_mask.nii');
        wskull_stripped{iSub} = fullfile(mask_paths{iSub},'rw_Skull_stripped_T1.nii');
        clear V BB vx
    end

    if jSub == 1
        spm('defaults','fmri');
        spm_get_defaults('cmdline',true);
        spm_jobman('run',batch{1});
    elseif jSub > 1
        try % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;           
            w = waitbar(0,'Please wait...','Name','Creating whole-brain masks','Tag','tmfc_waitbar');
            afterEach(D, @tmfc_parfor_waitbar);    
            tmfc_parfor_waitbar(w,jSub,1);
        end
        parfor (iSub = 1:jSub,M) % Run matlabbatches in parallel mode
            spm('defaults','fmri');
            spm_get_defaults('cmdline',true);
            spm_jobman('run',batch{iSub});
            try; send(D,[]); end % Update waitbar
        end
        try; delete(w); end % Close waitbar
    end
    clear batch

    % Save paths
    %----------------------------------------------------------------------
    masks.WB = wWB_mask;
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

