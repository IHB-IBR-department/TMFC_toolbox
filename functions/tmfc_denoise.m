function output_paths = tmfc_denoise(SPM_paths,options,struct_paths,funct_paths,display_FD,estimate_GLMs,clear_all)

% =[ Task-Modulated Functional Connectivity (TMFC) Denoise Toolbox v1.3 ]=
% 
% The TMFC denoise toolbox updates a selected general linear model with
% the addition of noise regressors. It can be used prior to TMFC analysis 
% (gPPI or BSC) or standard task activation analysis. General linear model
% should be specified and estimated in the SPM8/12/25 software (user needs
% to select the corresponding SPM.mat files).
% 
% The TMFC denoise toolbox requires that the first six unconvolved 
% regressors in each session represent head motion parameters. 
% 
% Extraction of BOLD signals from whole-brain, GM, WM, and CSF masks 
% requires structural T1 images in native space and unsmoothed, realigned,
% functional images in MNI space. If the source SPM.mat files specify paths
% to smoothed functional images, then unsmoothed functional images should
% be stored in the same folders and have a shorter prefix. 
%
% NOTE: All regressors specified in the original general linear model will
% be included in the updated model along with noise regressors. That is, 
% if the original model already contains expansion of head motion parameters
% or physiological regressors, they can be duplicated in the updated model.
% Thus, it is necessary to select models that include six standard motion
% regressors and other confound regressors that will not be calculated by
% the TMFC denoise toolbox.
%
% Functionality of the TMFC denoise toolbox:
%
% (1) Calculates head motion parameters (temporal derivatives and quadratic
%     terms). Temporal derivatives are calculated as backwards differences
%     (Van Dijk et al., 2012). Quadratic terms represent 6 squared motion
%     parameters and 6 squared temporal derivatives (Satterthwaite et al., 2013).
%
% (2) Calculates framewise displacement (FD) as the sum of the absolute values
%     of the derivatives of translational and rotational motion parameters
%     (Power et al., 2012).
%
% (3) Creates spike regressors based on a user-defined FD threshold. For each
%     flagged time point, a unit impulse function is included in general linear
%     model, which had the value of 1 at that time point and 0 elsewhere
%     (Lemieux et al., 2007; Satterthwaite et al., 2013).
%
% (4) Creates aCompCor regressors (Behzadi et al., 2007). Calculates fixed
%     number of principal components (PCs) or variable number of PCs
%     explaining 50% of the signal variability separately for eroded WM
%     and CSF masks (Muschelli et al., 2014).   
% 
% (5) Creates WM/CSF regressors (Fox and Raichle, 2007). Calculates average
%     BOLD signals separately for eroded WM and CSF masks. Optionally
%     calculates derivatives and quadratic terms (Parkes et al., 2017).
%
% (6) Creates GSR regressor (Fox et al, 2009). Calculates the average
%     BOLD signal for a whole-brain mask. Optionally calculates
%     derivatives and quadratic terms (Parkes et al., 2017).
%
% (7) Calculates Derivative of root mean square VARiance over voxelS (DVARS).
%     DVARS is computed as the root mean square (RMS) of the differentiated
%     BOLD time series within the GM mask (Muschelli et al., 2014).
%     Calculates FD/DVARS correlations. 
%     DVARS is computed before and after noise regression 
%     (for the original and updated GLM, respectively).
%
% (8) Adds noise regressiors to the original model and estimates it. The noise
%     regressors and the updated model will be stored in the TMFC_denoise subfolder.
%
% (9) Can use robust weighted least squares (rWLS) for model estimation.
%     It assumes that each image has an own variance parameter, i.e. some
%     scans may be disrupted by noise. By choosing this option, SPM will 
%     estimte the noise variances in the first pass and then re-weight each
%     image by the inverse of the variance in the second pass.
%
% -------------------------------------------------------------------------
% FORMAT: output_paths = tmfc_denoise
% Will call GUIs to select SPM.mat files, define denoising options, select
% structural and functional files, define FD threshold for spike regression.
%
% FORMAT: output_paths = tmfc_denoise(SPM_paths,options,struct_paths,funct_paths,display_FD,estimate_GLMs,clear_all)
% Performes noise regression without calling the GUI.
% 
% INPUTS: 
% SPM_paths             - Cell array containing paths to SPM.mat files that
%                         need to be re-estimated with noise regressors
%         
% options.motion        - '6HMP' : do not add additional motion regressors
%                       - '12HMP': add 6 temporal derivatives
%                       - '24HMP': add squared regressors
%
% options.rotation_indx - The order of motion regressors:
%                       - [4 5 6] (last three regressors - rotation, e.g., SPM12, HCP, fMRIPrep)
%                       - [1 2 3] (first three regressors - rotation, e.g., FSL, AFNI)
%
% options.rotation_unit - Rotation units:
%                       - 'rad' (radians, e.g., SPM12, FSL, fMRIPrep)
%                       - 'deg' (degrees, e.g., HCP, AFNI)
%
% options.head_radius   - Approximate head radius in mm (Default: 50)
%
% options.DVARS         - 0 (none) or 1 (calculate)
%
% options.aCompCor      - Number of the aCompCor regressors for WM mask
%                         and CSF mask (Default: [5 5]; None: [0 0];
%                         aCompCor50%: [0.5 0.5])
% options.aCompCor_ort  - Pre-orthogonalize WM and CSF signals w.r.t. high
%                         pass filter (HPF) regressors and head motion
%                         regressors prior to PC calculation (Default: 1)
%
% options.rWLS          - 0 (none) or 1 (apply rWLS)
%
% options.spikereg      - 0 (none) or 1 (add spike regressors)
% options.spikeregFDthr - FD threshold for creating spike regressors in mm (Default: 0.5)
%
% options.WM_CSF        - 'none' : do not add WM and CSF regressors
%                       - '2Phys': add WM and CSF signals 
%                       - '4Phys': add WM and CSF signals, and 2 temporal derivatives
%                       - '8Phys': add WM and CSF signals, 2 temporal derivatives and 4 quadratic terms
%
% options.GSR           - 'none': do not add whole-brain signal
%                       - 'GSR' : add whole-brain signal
%                       - '2GSR': add whole-brain signal and its temporal derivative
%                       - '4GSR': add whole-brain signal, its temporal derivative and 2 quadratic terms
%
% options.parallel      - 0 or 1: Sequential or parallel computations
%
% options.GMmask.prob   - Probability threshold for the liberal GM mask (Default: 0.95)
% options.WMmask.prob   - Probability threshold for the WM mask (Default: 0.99)
% options.CSFmask.prob  - Probability threshold for the CSF mask (Default: 0.99)
% options.GMmask.dilate - Dilation cycles for the GM mask (Default: 2 voxels)
% options.WMmask.erode  - Erosion cycles for the WM mask (Default: 2 voxels)
% options.CSFmask.erode - Erosion cycles for the CSF mask (Default: 2 voxels)
%
% struct_paths{iSub}.fname - Cell array containing paths to structural T1 images
% funct_paths{iSub}.fname  - Cell array containing paths to realigned, normalized and unsmoothed functional images
% (NOTE: Structural images should be in the native space. Functional images should be in the MNI space)
%
% display_FD    - 1 or 0 : Display individual FD plots (Default: 1)
% estimate_GLMs - 1 or 0 : Estimate GLMs with noise regeressors (Default: 1)
% clear_all     - 1 or 0 : Delete all already created files in 'TMFC_denoise'
%                          subfolders before creating new files (Default: 0)  
%
% OUTPUT: 
% output_paths  - Cell array containing paths to estimated GLMs
%                 with added noise regeressors. 
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

disp('=============[TMFC denoise: v1.3]=============');

output_paths = [];

%-Prepare inputs
%--------------------------------------------------------------------------

% Select SPM.mat files
if ~nargin
    [SPM_paths,subject_paths] = tmfc_select_subjects_GUI(0);
end
if isempty(SPM_paths); error('Select SPM.mat files.'); end

% Check SPM.mat files
for iSub = 1:length(SPM_paths)
    if ~exist(SPM_paths{iSub},'file')
        error(['SPM file not found: ' SPM_paths{iSub}]); 
    end
end

% Subject paths
if ~exist('subject_paths','var')
    subject_paths = tmfc_SPM_subject_paths(SPM_paths);
end

% Define denoising options
if nargin<2
    options = tmfc_denoise_options_GUI;
end
if isempty(options); error('Denoising options not selected.'); end

% Select structural T1 images
if nargin<3
    if (sum(options.aCompCor)~=0 || ~strcmp(options.WM_CSF,'none') || ~strcmp(options.GSR,'none') || options.DVARS == 1)
        struct_paths = tmfc_select_struct_GUI(subject_paths);
        if isempty(struct_paths); error('Select structural T1 files.'); end
    else
        struct_paths = [];
    end
end

% Select realigned and unsmoothed functional images
if nargin<4
    if (sum(options.aCompCor)~=0 || ~strcmp(options.WM_CSF,'none') || ~strcmp(options.GSR,'none') || options.DVARS == 1)
        funct_paths = tmfc_select_funct_GUI(SPM_paths,subject_paths);
        if isempty(funct_paths); error('Select unsmoothed functional files.'); end
    else
        funct_paths = [];
    end
end

% Display individual FD plots
if nargin<5
    display_FD = 1;
end

% Estimate GLMs with noise regressors 
if nargin<6
    estimate_GLMs = 1;
end

% Clear "TMFC_denoise" subfolders
if nargin<7
    answer = questdlg('Delete previosly created TMFC denoise files?', ...
    'TMFC denoise', ...
    'Do not delete','Delete','Do not delete');
    switch answer
        case 'Do not delete'
            clear_all = 0;
        case 'Delete'
            clear_all = 1;
    end
end

%-Create TMFC denoise subfolders
%--------------------------------------------------------------------------
for iSub = 1:length(SPM_paths)
    GLM_subfolder = fileparts(SPM_paths{iSub});
    if clear_all == 1 && exist(fullfile(GLM_subfolder,'TMFC_denoise'),'dir')
        rmdir(fullfile(GLM_subfolder,'TMFC_denoise'),'s');
    end
    if ~exist(fullfile(GLM_subfolder,'TMFC_denoise'),'dir')
        mkdir(fullfile(GLM_subfolder,'TMFC_denoise'));
    end
    clear GLM_subfolder
end

%-Calculate head motion parameters (HMP) and framewise displacement (FD)
%--------------------------------------------------------------------------
if ~strcmp(options.motion,'6HMP') || options.spikereg == 1 || display_FD == 1 || options.DVARS == 1
    disp('Head motion assessment...'); tic;
    FD = tmfc_head_motion(SPM_paths,subject_paths,options);
    hms = fix(mod((toc), [0, 3600, 60]) ./ [3600, 60, 1]);
    disp(['Done in ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec].']);
end

%-Plot FD time series and select FDthr for spike regression
%--------------------------------------------------------------------------
if display_FD == 1
    FDthr = tmfc_plot_FD(FD,options,SPM_paths,subject_paths,struct_paths,funct_paths);
    options.spikeregFDthr = FDthr;
end

%-Create spike regressors
%--------------------------------------------------------------------------
if options.spikereg == 1
    disp('----------------------------------------');
    disp('Creating spike regressors...'); tic;
    tmfc_spikereg(SPM_paths,options);
    hms = fix(mod((toc), [0, 3600, 60]) ./ [3600, 60, 1]);
    disp(['Done in ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec].']);
end

%-Create GM/WM/CSF and whole-brain masks
%--------------------------------------------------------------------------
if sum(options.aCompCor)~=0 || ~strcmp(options.WM_CSF,'none') || ~strcmp(options.GSR,'none')  || options.DVARS == 1
    if isempty(struct_paths); error('Select structural T1 files.'); end
    if isempty(funct_paths); error('Select unsmoothed functional files.'); end
    disp('----------------------------------------');
    if ~isfield(options,'GMmask') || ~isfield(options,'WMmask') || ~isfield(options,'CSFmask') 
        [options.GMmask.prob, options.WMmask.prob, options.CSFmask.prob, ...
         options.GMmask.dilate, options.WMmask.erode, options.CSFmask.erode] = tmfc_masks_GUI();
    end
    disp('Creating binary masks...'); tic;
    masks = tmfc_create_masks(SPM_paths,struct_paths,funct_paths,options);
    hms = fix(mod((toc), [0, 3600, 60]) ./ [3600, 60, 1]);
    disp(['Done in ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec].']);
end

%-Calculate physiological regressors
%--------------------------------------------------------------------------
if sum(options.aCompCor)~=0 || ~strcmp(options.WM_CSF,'none') || ~strcmp(options.GSR,'none')
    disp('----------------------------------------');
    disp('Calculating physiological regressors...'); tic;
    tmfc_physioreg(SPM_paths,subject_paths,funct_paths,masks,options);
    hms = fix(mod((toc), [0, 3600, 60]) ./ [3600, 60, 1]);
    disp(['Done in ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec].']);
end

%-Estimate updated GLMs with noise regressors
%--------------------------------------------------------------------------
output_paths = [];
if estimate_GLMs == 1
    disp('----------------------------------------');
    disp('Estimating GLMs with noise regressors...'); tic;
    if ~exist('masks','var'); masks = []; end
    output_paths = tmfc_estimate_updated_GLMs(SPM_paths,masks,options);
    hms = fix(mod((toc), [0, 3600, 60]) ./ [3600, 60, 1]);
    disp(['Done in ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec].']);
end

%-Calculate and plot DVARS
%--------------------------------------------------------------------------
if options.DVARS == 1
    disp('----------------------------------------');
    disp('Calculating DVARS...'); tic;
    [preDVARS,postDVARS] = tmfc_calculate_DVARS(FD,SPM_paths,options,masks,output_paths);
    tmfc_plot_DVARS(preDVARS,postDVARS,FD,options,SPM_paths,subject_paths,struct_paths,funct_paths,masks);
    hms = fix(mod((toc), [0, 3600, 60]) ./ [3600, 60, 1]);
    disp(['Done in ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec].']);
end

end

%%============================[Subfunctions]===============================
function [subject_paths] = tmfc_SPM_subject_paths(SPM_paths)
subfolders = regexp(SPM_paths, filesep, 'split');
if length(SPM_paths) > 1
    nLevel1 = length(subfolders{1});
    nLevel2 = length(subfolders{2});
    for iLevel = 1:nLevel1
        if ~strcmp(subfolders{1}{nLevel1-iLevel},subfolders{2}{nLevel2-iLevel})
            break;
        end
    end
    for iSub = 1:length(SPM_paths)
        nLevel = length(subfolders{iSub});
        subject_paths{iSub,1} = fullfile(subfolders{iSub}{1:nLevel-iLevel});
    end
else
    subject_paths{1,1} = spm_select(1,'dir','Select "Subject" folder');
end
end
