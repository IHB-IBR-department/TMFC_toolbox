function tmfc_prepare_example_data(data,parallel)

% ========================================================================
% Ruslan Masharipov, 2025
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

% Experiment folder
exp_folder = ['Example_data_SF_[' num2str(data.SF,'%.2f') ']_SNR_[' num2str(data.SNR,'%.2f') ']_STP_[' num2str(data.STP_delay,'%.2f') ']_'  data.model];

% Apply STP delay
load(data.sots_path);
onsets{1,1} = onsets{1,1} - data.STP_delay;
onsets{1,2} = onsets{1,2} - data.STP_delay;

sots_path = data.sots_path;
sots_path = sots_path(1:end-4);
sots_path = [sots_path '_[' num2str(data.STP_delay,'%.2f') 's_STP].mat'];
[~,sots] = fileparts(sots_path);
sots_path = fullfile(data.stat_path,exp_folder,[sots '.mat']);
mkdir(fullfile(data.stat_path,exp_folder));
save(sots_path,'activations','onsets','durations','names','rest_matrix','task_matrices');

% Generate .nii functional images
fprintf('Generating *.nii images... \n');
tmfc_generate_funct_images(data.stat_path,data.sim_path,exp_folder,data.SF,data.SNR,data.N,data.N_ROIs,data.dummy)

% Generate .nii ROI binary masks 
fprintf('Generating ROI masks... \n');
tmfc_generate_ROI_masks(data.stat_path,exp_folder,data.N_ROIs)

% Estimate GLM
tic
fprintf('Estimating basic GLMs... \n');
tmfc_estimate_GLM(data.stat_path,sots_path,exp_folder,data.N,data.TR,data.model,parallel)
fprintf(['Estimating basic GLMs :: Done in: ' num2str(toc) 's \n']);