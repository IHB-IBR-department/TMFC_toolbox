clc
clear
close all

% BEFORE RUNNING THIS SCRIPT:
% 1) Set path to SPM12
% 2) Set path to TMFC_toolbox (Add with subfolders)
% 3) Change current working directory to: '...\TMFC_toolbox\examples'

cd(fileparts(matlab.desktop.editor.getActiveFilename)); % Set path to '...\TMFC_toolbox\examples'

%% Prepare example data and calculate basic first-level GLMs

data.SF  = 1;         % Scaling Factor (SF) for co-activations: SF = SD_oscill/SD_coact
data.SNR = 1;         % Signal-to-noise ratio (SNR): SNR = SD_signal/SD_noise
data.STP_delay = 0.2; % Short-term synaptic plasticity (STP) delay, [s]
data.N = 20;          % Sample size (Select 20 subjects out of 100 to reduce computations)
data.N_ROIs = 100;    % Number of ROIs
data.dummy = 3;       % Remove first M dummy scans
data.TR = 2;          % Repetition time (TR), [s]
data.model = 'AR(1)'; % Autocorrelation modeling

% Set path for stat folder 
spm_jobman('initcfg');
data.stat_path = spm_select(1,'dir','Select a folder for data extraction and statistical analysis');

% Set path for simulated BOLD time series *.mat file
data.sim_path = fullfile(pwd,'data','SIMULATED_BOLD_EVENT_RELATED_[2s_TR]_[1s_DUR]_[6s_ISI]_[40_TRIALS].mat');

% Set path for task design *.mat file (stimulus onset times, SOTs)
data.sots_path = fullfile(pwd,'data','TASK_DESIGN_EVENT_RELATED_[2s_TR]_[1s_DUR]_[6s_ISI]_[40_TRIALS].mat');

% Generate *.nii images and calculate GLMs
prepare_example_data(data)

% Change current directory to new TMFC project folder
cd(data.stat_path)


%% Setting up computation parameters

% Sequential or parallel computing (0 or 1)
tmfc.defaults.parallel = 1;         % Parallel
% Store temporaty files during GLM estimation in RAM or on disk
tmfc.defaults.resmem = true;        % RAM
% How much RAM can be used at the same time during GLM estimation
tmfc.defaults.maxmem = 2^33;        % 8 GB
% Seed-to-voxel and ROI-to-ROI analyses
tmfc.defaults.analysis = 1;


%% Setting up paths

% The path where all results will be saved
tmfc.project_path = data.stat_path;

% Define paths to individual subject SPM.mat files
% tmfc.subjects(1).path = '...\Your_study\Subjects\sub_001\stat\Standard_GLM\SPM.mat';
% tmfc.subjects(2).path = '...\Your_study\Subjects\sub_002\stat\Standard_GLM\SPM.mat';
% tmfc.subjects(3).path = '...\Your_study\Subjects\sub_003\stat\Standard_GLM\SPM.mat';
% etc

% Alternativelly, use tmfc_select_subjects_GUI to select subjects
% Go to GLMs subfolder and select 20 subjects 
SPM_check = 1;                      % Check SPM.mat files
[paths] = tmfc_select_subjects_GUI(SPM_check);

for i = 1:length(paths)
    tmfc.subjects(i).path = paths{i};
end

clear SPM_check paths

%% Select ROIs

% Use tmfc_select_ROIs_GUI to select ROIs
%
% The tmfc_select_ROIs_GUI function creates group binary mask based on
% 1st-level masks (SPM.VM) and applies it to all selected ROIs. Empty ROIs
% will be removed. Masked ROIs will be limited to only voxels which have 
% data for all subjects. The dimensions, orientation, and voxel sizes of 
% the masked ROI images will be adjusted according to the group binary mask
%
% Go to ROI_masks subfolder and select 100 ROIs

[ROI_set] = tmfc_select_ROIs_GUI(tmfc);
tmfc.ROI_set(1) = ROI_set;

clear ROI_set


%% LSS regression

% Define conditions of interest
% tmfc.LSS.conditions(1).sess   = 1;   
% tmfc.LSS.conditions(1).number = 1;
% tmfc.LSS.conditions(2).sess   = 1;
% tmfc.LSS.conditions(2).number = 2;

% Alternatively, use tmfc_LSS_GUI to select conditions of interest
[conditions] = tmfc_LSS_GUI(tmfc.subjects(1).path);
tmfc.LSS.conditions = conditions;

% Run LSS regression
start_sub = 1;                      % Start from the 1st subject
[sub_check] = tmfc_LSS(tmfc,start_sub);

clear conditions


%% BSC-LSS

% Extract and correlate mean beta series for conditions of interest
ROI_set_number = 1;                 % Select ROI set
[sub_check,contrasts] = tmfc_BSC(tmfc,ROI_set_number);

% Update contrasts info
% The tmfc_BSC function creates default contrasts for each
% condition of interest (i.e., Condition > Baseline)
tmfc.ROI_set(ROI_set_number).contrasts.BSC = contrasts;

% Define new contrasts:
tmfc.ROI_set(ROI_set_number).contrasts.BSC(3).title = 'TaskA_vs_TaskB';
tmfc.ROI_set(ROI_set_number).contrasts.BSC(4).title = 'TaskB_vs_TaskA';
tmfc.ROI_set(ROI_set_number).contrasts.BSC(3).weights = [1 -1];
tmfc.ROI_set(ROI_set_number).contrasts.BSC(4).weights = [-1 1];

% Calculate new contrasts
type = 3;                           % BSC-LSS
contrast_number = [3,4];            % Calculate contrasts #3 and #4
[sub_check] = tmfc_ROI_to_ROI_contrast(tmfc,type,contrast_number,ROI_set_number);
[sub_check] = tmfc_seed_to_voxel_contrast(tmfc,type,contrast_number,ROI_set_number);

% Load BSC-LSS matrices for the 'TaskA_vs_TaskB' contrast (contrast # 3)
for i = 1:data.N 
    M(i).paths = struct2array(load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS','ROI_to_ROI',...
        ['Subject_' num2str(i,'%04.f') '_Contrast_0003_[TaskA_vs_TaskB].mat'])));
end
matrices = cat(3,M(:).paths);

% Perform one-sample t-test (two-sided, FDR-correction) 
contrast = 1;                       % TaskA_vs_TaskB effect
alpha = 0.001/2;                    % alpha = 0.001 thredhold corrected for two-sided comparison
correction = 'FDR';                 % False Discovery Rate (FDR) correction (Benjamini–Hochberg procedure)
[thresholded,pval,tval,conval] = tmfc_ttest(matrices,contrast,alpha,correction); 

% Plot BSC-LSS results
figure(1);
try
    sgtitle('BSC-LSS results');
catch
    suptitle('BSC-LSS results');
end
subplot(1,2,1); imagesc(conval);        title('Group mean'); axis square; colorbar; caxis(tmfc_axis(conval,1));
subplot(1,2,2); imagesc(thresholded);   title('pFDR<0.001'); axis square; colorbar;
colormap(subplot(1,2,2),'parula')
colormap(subplot(1,2,1),'redblue')
set(findall(gcf,'-property','FontSize'),'FontSize',16)

clear type contrasts contrast_number


%% FIR task regression (regress out co-activations and save residual time series)

% FIR window length in [s]
tmfc.FIR.window = 24;
% Number of FIR time bins
tmfc.FIR.bins = 24;

% Run FIR task regression
[sub_check] = tmfc_FIR(tmfc,start_sub);


%% LSS regression after FIR task regression (use residual time series)

% Define conditions of interest
tmfc.LSS_after_FIR.conditions = tmfc.LSS.conditions;

% Run LSS regression
[sub_check] = tmfc_LSS_after_FIR(tmfc,start_sub);


%% BSC-LSS after FIR task regression (use residual time series)

% Extract and correlate mean beta series for conditions of interest
ROI_set_number = 1;                 % Select ROI set
[sub_check,contrasts] = tmfc_BSC_after_FIR(tmfc,ROI_set_number);

% Update contrasts info
% The tmfc_BSC_after_FIR function creates default contrasts for each
% condition of interest (i.e., Condition > Baseline)
tmfc.ROI_set(ROI_set_number).contrasts.BSC_after_FIR = contrasts;

% Define new contrast:
tmfc.ROI_set(ROI_set_number).contrasts.BSC_after_FIR(3).title = 'TaskA_vs_TaskB';
tmfc.ROI_set(ROI_set_number).contrasts.BSC_after_FIR(4).title = 'TaskB_vs_TaskA';
tmfc.ROI_set(ROI_set_number).contrasts.BSC_after_FIR(3).weights = [1 -1];
tmfc.ROI_set(ROI_set_number).contrasts.BSC_after_FIR(4).weights = [-1 1];

% Calculate new contrast
type = 4;                           % BSC-LSS after FIR
contrast_number = [3,4];            % Calculate contrast #3 and #4
[sub_check] = tmfc_ROI_to_ROI_contrast(tmfc,type,contrast_number,ROI_set_number);
[sub_check] = tmfc_seed_to_voxel_contrast(tmfc,type,contrast_number,ROI_set_number);

% Load BSC-LSS (after FIR) matrices for the 'TaskA_vs_TaskB' contrast (contrast # 3)
clear M marices conval thresholded
for i = 1:data.N 
    M(i).paths = struct2array(load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','ROI_to_ROI',...
        ['Subject_' num2str(i,'%04.f') '_Contrast_0003_[TaskA_vs_TaskB].mat'])));
end
matrices = cat(3,M(:).paths);

% Perform one-sample t-test (two-sided, FDR-correction) 
[thresholded,pval,tval,conval] = tmfc_ttest(matrices,contrast,alpha,correction); 

% Plot BSC-LSS (after FIR) results
figure(2);
try
    sgtitle('BSC-LSS (after FIR task regression) results');
catch
    suptitle('BSC-LSS (after FIR task regression) results');
end
subplot(1,2,1); imagesc(conval);        title('Group mean'); axis square; colorbar; caxis(tmfc_axis(conval,1));
subplot(1,2,2); imagesc(thresholded);   title('pFDR<0.001'); axis square; colorbar;
colormap(subplot(1,2,2),'parula') 
colormap(subplot(1,2,1),'redblue')
set(findall(gcf,'-property','FontSize'),'FontSize',16)

clear type contrasts contrast_number

%% BGFC

% Calculate background functional connectivity (BGFC)
[sub_check] = tmfc_BGFC(tmfc,ROI_set_number,start_sub);


%% gPPI

% Define conditions of interest
[conditions] = tmfc_gPPI_GUI(tmfc.subjects(1).path);
tmfc.ROI_set(ROI_set_number).gPPI.conditions = conditions;
clear conditions

% VOI extraction
[sub_check] = tmfc_VOI(tmfc,ROI_set_number,start_sub);

% PPI calculation
[sub_check] = tmfc_PPI(tmfc,ROI_set_number,start_sub);

% gPPI calculation
[sub_check,contrasts] = tmfc_gPPI(tmfc,ROI_set_number,start_sub);

% Update contrasts info
% The tmfc_gPPI function creates default contrasts for each
% condition of interest (i.e., Condition > Baseline)
tmfc.ROI_set(ROI_set_number).contrasts.gPPI = contrasts;

% Define new contrasts:
tmfc.ROI_set(ROI_set_number).contrasts.gPPI(3).title = 'TaskA_vs_TaskB';
tmfc.ROI_set(ROI_set_number).contrasts.gPPI(4).title = 'TaskB_vs_TaskA';
tmfc.ROI_set(ROI_set_number).contrasts.gPPI(3).weights = [1 -1];
tmfc.ROI_set(ROI_set_number).contrasts.gPPI(4).weights = [-1 1];

% Calculate new contrasts
type = 1;                           % gPPI
contrast_number = [3,4];            % Calculate contrasts #3 and #4
[sub_check] = tmfc_ROI_to_ROI_contrast(tmfc,type,contrast_number,ROI_set_number);
[sub_check] = tmfc_seed_to_voxel_contrast(tmfc,type,contrast_number,ROI_set_number);

% Load gPPI matrices for the 'TaskA_vs_TaskB' contrast (contrast # 3)
clear M marices conval thresholded
for i = 1:data.N 
    M(i).paths = struct2array(load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI','ROI_to_ROI','symmetrical',...
        ['Subject_' num2str(i,'%04.f') '_Contrast_0003_[TaskA_vs_TaskB].mat'])));
end
matrices = cat(3,M(:).paths);

% Perform one-sample t-test (two-sided, FDR-correction) 
[thresholded,pval,tval,conval] = tmfc_ttest(matrices,contrast,alpha,correction); 

% Plot gPPI results
figure(3);
try
    sgtitle('gPPI results');
catch
    suptitle('gPPI results');
end
subplot(1,2,1); imagesc(conval);        title('Group mean'); axis square; colorbar; caxis(tmfc_axis(conval,1));
subplot(1,2,2); imagesc(thresholded);   title('pFDR<0.001'); axis square; colorbar;
colormap(subplot(1,2,2),'parula')
colormap(subplot(1,2,1),'redblue')
set(findall(gcf,'-property','FontSize'),'FontSize',16)

clear type contrasts contrast_number


%% gPPI-FIR (gPPI model with psychological regressors defined by FIR functions)

% Define FIR parameters for gPPI-FIR
tmfc.ROI_set(ROI_set_number).gPPI_FIR.window = 24;   % FIR window length in [s]
tmfc.ROI_set(ROI_set_number).gPPI_FIR.bins = 24;     % Number of FIR time bins

% gPPI-FIR calculation
[sub_check,contrasts] = tmfc_gPPI_FIR(tmfc,ROI_set_number,start_sub);

% Update contrasts info
% The tmfc_gPPI_FIR function creates default contrasts for each
% condition of interest (i.e., Condition > Baseline)
tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR = contrasts;

% Define new contrasts:
tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(3).title = 'TaskA_vs_TaskB';
tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(4).title = 'TaskB_vs_TaskA';
tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(3).weights = [1 -1];
tmfc.ROI_set(ROI_set_number).contrasts.gPPI_FIR(4).weights = [-1 1];

% Calculate new contrasts
type = 2;                           % gPPI-FIR
contrast_number = [3,4];            % Calculate contrasts #3 and #4
[sub_check] = tmfc_ROI_to_ROI_contrast(tmfc,type,contrast_number,ROI_set_number);
[sub_check] = tmfc_seed_to_voxel_contrast(tmfc,type,contrast_number,ROI_set_number);

% Load gPPI-FIR matrices for the 'TaskA_vs_TaskB' contrast (contrast # 3)
clear M marices conval thresholded
for i = 1:data.N 
    M(i).paths = struct2array(load(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'gPPI_FIR','ROI_to_ROI','symmetrical',...
        ['Subject_' num2str(i,'%04.f') '_Contrast_0003_[TaskA_vs_TaskB].mat'])));
end
matrices = cat(3,M(:).paths);

% Perform one-sample t-test (two-sided, FDR-correction) 
[thresholded,pval,tval,conval] = tmfc_ttest(matrices,contrast,alpha,correction); 

% Plot gPPI-FIR results
figure(4);
try
    sgtitle('gPPI-FIR results');
catch
    suptitle('gPPI-FIR results');
end
subplot(1,2,1); imagesc(conval);        title('Group mean'); axis square; colorbar; caxis(tmfc_axis(conval,1));
subplot(1,2,2); imagesc(thresholded);   title('pFDR<0.001'); axis square; colorbar;
colormap(subplot(1,2,2),'parula')
colormap(subplot(1,2,1),'redblue')
set(findall(gcf,'-property','FontSize'),'FontSize',16)

clear type contrasts contrast_number

%% Save TMFC project *.mat file
save(fullfile(data.stat_path,'TMFC_project.mat'),'tmfc');