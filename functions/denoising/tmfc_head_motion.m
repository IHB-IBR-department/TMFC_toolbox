function FD = tmfc_head_motion(SPM_paths,subject_paths,options)

% =======[ Task-Modulated Functional Connectivity Denoise Toolbox ]========
% 
% (1) Calculates head motion parameters (temporal derivatives and quadratic
%     terms). Temporal derivatives are calculated as backward differences
%     (Van Dijk et al., 2012). Quadratic terms represent 6 squared motion
%     parameters and 6 squared temporal derivatives (Satterthwaite et al., 2013).
%
% (2) Calculates framewise displacement (FD) as the sum of absolute
%     derivatives of translational and rotational parameters
%     (Power et al., 2012).
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

% Load head motion parameters (HMP)
%--------------------------------------------------------------------------
disp('Loading head motion parameters...');
for iSub = 1:length(SPM_paths)
    SPM = load(SPM_paths{iSub}).SPM;
    for jSess = 1:length(SPM.Sess)
        if size(SPM.Sess(jSess).C.C,2) < 6
            error('The original model contains fewer than six confound regressors. It must include six head motion regressors. Please check:\n%s',SPM_paths{iSub});
        else
            group_HMP(iSub).Sess{jSess} = SPM.Sess(jSess).C.C(:,[options.translation_idx options.rotation_idx]);
        end
    end
    clear SPM
end

% Calculate 12HMP, 24HMP, and framewise displacement (FD)
%--------------------------------------------------------------------------
for iSub = 1:length(group_HMP)
    GLM_subfolder = fileparts(SPM_paths{iSub});
    if ~exist(fullfile(GLM_subfolder,'TMFC_denoise'),'dir')
        mkdir(fullfile(GLM_subfolder,'TMFC_denoise'));
    end
    sub_FD_max = []; sub_FD_mean = [];
    for jSess = 1:length(group_HMP(iSub).Sess)
        HMP = group_HMP(iSub).Sess{jSess};
        HMP_diff = [zeros(1,6); diff(HMP)];
        HMP12(jSess).Sess = [HMP HMP_diff];
        HMP24(jSess).Sess = [HMP HMP_diff HMP.^2 HMP_diff.^2];

        % Calculate FD
        HMP_diff_xyz = HMP_diff(:,1:3);
        HMP_diff_rot = HMP_diff(:,4:6);

        if strcmpi(options.rotation_unit,'rad')        % Convert radians to mm
            HMP_diff_rot = options.head_radius*HMP_diff_rot;
        elseif strcmpi(options.rotation_unit,'deg')    % Convert degrees to mm
            HMP_diff_rot = options.head_radius*pi/180*HMP_diff_rot;
        end
        
        ts_FD = sum(abs(HMP_diff_xyz),2)+sum(abs(HMP_diff_rot),2);
        sub_FD_mean =   [sub_FD_mean   mean(ts_FD)];
        sub_FD_max =    [sub_FD_max    max(ts_FD)];
        
        FD(iSub).SPM_path = SPM_paths{iSub};
        [~, sub, ~] = fileparts(subject_paths{iSub});
        FD(iSub).Subject = sub;
        FD(iSub).Sess(jSess).FD_ts = ts_FD;
        FD(iSub).Sess(jSess).mean = mean(ts_FD);
        FD(iSub).Sess(jSess).max = max(ts_FD);

        clear HMP HMP_diff* ts_FD 
    end

    FD(iSub).FD_mean = mean(sub_FD_mean);
    FD(iSub).FD_max = max(sub_FD_max);
    clear sub_FD_max sub_FD_mean sub

    % Save 12HMP.mat and 24HMP.mat files
    if strcmpi(options.motion,'12HMP')
        if ~exist(fullfile(GLM_subfolder,'TMFC_denoise','12HMP.mat'),'file')
            save(fullfile(GLM_subfolder,'TMFC_denoise','12HMP.mat'),'HMP12');
        end
    elseif strcmpi(options.motion,'24HMP')
        if ~exist(fullfile(GLM_subfolder,'TMFC_denoise','24HMP.mat'),'file')
            save(fullfile(GLM_subfolder,'TMFC_denoise','24HMP.mat'),'HMP24');
        end
    end
    clear HMP12 HMP24
    
    % Save FD.mat files
    if ~exist(fullfile(GLM_subfolder,'TMFC_denoise','FD.mat'),'file')
        FramewiseDisplacement = FD(iSub);
        save(fullfile(GLM_subfolder,'TMFC_denoise','FD.mat'),'FramewiseDisplacement');
        clear FramewiseDisplacement
    end
    clear GLM_subfolder
end
end
