function [sub_check,contrasts] = tmfc_BSC_after_FIR(tmfc,ROI_set_number,clear_BSC)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Extracts average (mean or first eigenvariate) beta series from selected
% ROIs. Correlates beta series for conditions of interest. Saves individual
% correlational matrices (ROI-to-ROI analysis) and correlational images
% (seed-to-voxel analysis)for each condition of interest. These refer to
% default contrasts, which can then be multiplied by linear contrast weights.
%
% FORMAT [sub_check,contrasts] = tmfc_BSC_after_FIR(tmfc)
%
%   tmfc.subjects.path     - Paths to individual SPM.mat files
%   tmfc.subjects.name     - Subject names within the TMFC project
%                           ('Subject_XXXX' naming will be used by default)
%   tmfc.project_path      - Path where all results will be saved
%   tmfc.defaults.analysis - 1 (Seed-to-voxel and ROI-to-ROI analyses)
%                          - 2 (ROI-to-ROI analysis only)
%                          - 3 (Seed-to-voxel analysis only)
%
%   tmfc.ROI_set                  - List of selected ROIs
%   tmfc.ROI_set.BSC_after_FIR    - 'mean' or 'first_eigenvariate'(default)
%   tmfc.ROI_set.type             - Type of the ROI set
%   tmfc.ROI_set.set_name         - Name of the ROI set
%   tmfc.ROI_set.ROIs.name        - Name of the selected ROI
%   tmfc.ROI_set.ROIs.path_masked - Paths to the ROI images masked by group
%                                   mean binary mask 
%
%   tmfc.LSS_after_FIR.conditions        - List of conditions of interest
%   tmfc.LSS_after_FIR.conditions.sess   - Session number
%                                          (as specified in SPM.Sess)
%   tmfc.LSS_after_FIR.conditions.number - Condition number 
%                                          (as specified in SPM.Sess.U)
%   tmfc.LSS_after_FIR.conditions.name   - Condition name
%                                          (as specified in SPM.Sess.U.name)
%   tmfc.LSS_after_FIR.conditions.file_name - Condition-specific file names:
%   (['[Sess_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
%    regexprep(char(SPM.Sess(iSess).U(jCond).name(1)),' ','_') ']'];)
%
% Session number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two sessions. Both sessions 
% contains three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified (see tmfc_conditions_GUI, nested function:
% [cond_list] = generate_conditions(SPM_path)):
%
%   tmfc.LSS_after_FIR.conditions(1).sess   = 1;   
%   tmfc.LSS_after_FIR.conditions(1).number = 1; 
%   tmfc.LSS_after_FIR.conditions(1).name = 'Cond_A'; 
%   tmfc.LSS_after_FIR.conditions(1).file_name = '[Sess_1]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS_after_FIR.conditions(2).sess   = 1;
%   tmfc.LSS_after_FIR.conditions(2).number = 2;
%   tmfc.LSS_after_FIR.conditions(2).name = 'Cond_B';
%   tmfc.LSS_after_FIR.conditions(2).file_name = '[Sess_1]_[Cond_2]_[Cond_B]';  
%
%   tmfc.LSS_after_FIR.conditions(3).sess   = 2;
%   tmfc.LSS_after_FIR.conditions(3).number = 1;
%   tmfc.LSS_after_FIR.conditions(3).name = 'Cond_A';
%   tmfc.LSS_after_FIR.conditions(3).file_name = '[Sess_2]_[Cond_1]_[Cond_A]';  
%
%   tmfc.LSS_after_FIR.conditions(4).sess   = 2;
%   tmfc.LSS_after_FIR.conditions(4).number = 2;
%   tmfc.LSS_after_FIR.conditions(4).name = 'Cond_B';
%   tmfc.LSS_after_FIR.conditions(4).file_name = '[Sess_2]_[Cond_2]_[Cond_B]';  
%
% Example of the ROI set (see tmfc_select_ROIs_GUI):
%
%   tmfc.ROI_set(1).set_name = 'two_ROIs';
%   tmfc.ROI_set(1).type = 'binary_images';
%   tmfc.ROI_set(1).ROIs(1).name = 'ROI_1';
%   tmfc.ROI_set(1).ROIs(2).name = 'ROI_2';
%   tmfc.ROI_set(1).ROIs(1).path_masked = 'C:\ROI_set\two_ROIs\ROI_1.nii';
%   tmfc.ROI_set(1).ROIs(2).path_masked = 'C:\ROI_set\two_ROIs\ROI_2.nii';
%
% FORMAT [sub_check,contrasts] = tmfc_BSC_after_FIR(tmfc,ROI_set_number)
% Run the function for the selected ROI set.
%
%   tmfc                   - As above
%   ROI_set_number         - Number of the ROI set in the tmfc structure
%                            (by default, ROI_set_number = 1)
%
% FORMAT [sub_check,contrasts] = tmfc_BSC_after_FIR(tmfc,ROI_set_number,clear_BSC)
% Run the function for the selected ROI set.
%
%   clear_BSC              - Clear previosly created BSC (after FIR) folders
%                            (0 - do not clear, 1 - clear)
%                            (by default, clear_BSC = 1)
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
    
if nargin == 1
    ROI_set_number = 1;
    clear_BSC = 1;
elseif nargin == 2
    clear_BSC = 1;
end

if ~isfield(tmfc.ROI_set(ROI_set_number),'BSC_after_FIR')
    tmfc.ROI_set(ROI_set_number).BSC_after_FIR = 'first_eigenvariate';
elseif isempty(tmfc.ROI_set(ROI_set_number).BSC_after_FIR)
    tmfc.ROI_set(ROI_set_number).BSC_after_FIR = 'first_eigenvariate';
end

% Check subject names
if ~isfield(tmfc.subjects,'name')
    for iSub = 1:length(tmfc.subjects)
        tmfc.subjects(iSub).name = ['Subject_' num2str(iSub,'%04.f')];
    end
end

nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
nSub = length(tmfc.subjects);
cond_list = tmfc.LSS_after_FIR.conditions;
nCond = length(cond_list);

% Clear previosly created BSC (after FIR) folders
if clear_BSC == 1
    if isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR'))
        rmdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR'),'s');
    end
end

% Create BSC (after FIR) folders
if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','Beta_series'))
    mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','Beta_series'));
end

if tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 2
    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','ROI_to_ROI'))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','ROI_to_ROI'));
    end
end

if tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 3
    for iROI = 1:nROI
        if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','Seed_to_voxel',tmfc.ROI_set(ROI_set_number).ROIs(iROI).name))
            mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','Seed_to_voxel',tmfc.ROI_set(ROI_set_number).ROIs(iROI).name));
        end
    end
end

SPM = load(tmfc.subjects(1).path).SPM;
XYZ  = SPM.xVol.XYZ;
iXYZ = cumprod([1,SPM.xVol.DIM(1:2)'])*XYZ - sum(cumprod(SPM.xVol.DIM(1:2)'));
hdr.dim = SPM.Vbeta(1).dim;
hdr.dt = SPM.Vbeta(1).dt;
hdr.pinfo = SPM.Vbeta(1).pinfo;
hdr.mat = SPM.Vbeta(1).mat;

% Loading ROIs
switch tmfc.ROI_set(ROI_set_number).type
    case {'binary_images','fixed_spheres'}
        w = waitbar(0,'Please wait...','Name','Loading ROIs');
        for iROI = 1:nROI
            ROIs(iROI).mask = spm_data_read(spm_data_hdr_read(tmfc.ROI_set(ROI_set_number).ROIs(iROI).path_masked),'xyz',XYZ);
            ROIs(iROI).mask(ROIs(iROI).mask == 0) = NaN;
            try
                waitbar(iROI/nROI,w,['ROI No ' num2str(iROI,'%.f')]);
            end
        end
        try
            delete(w)
        end
    otherwise
        ROIs = [];
end

% Sequential or parallel computing
switch tmfc.defaults.parallel
    % ----------------------- Sequential Computing ------------------------
    case 0

        % Create waitbar
        w = waitbar(0,'Please wait...','Name','Extract and correlate beta series','Tag','tmfc_waitbar');
        start_time = tic;
        count_sub = 1;
        cleanupObj = onCleanup(@unfreeze_after_ctrl_c);  

        for iSub = 1:nSub
            tmfc_extract_betas_after_FIR(tmfc,ROI_set_number,ROIs,nROI,nCond,cond_list,XYZ,iXYZ,hdr,iSub);
            sub_check(iSub) = 1;

            % Update waitbar
            elapsed_time = toc(start_time);
            time_per_sub = elapsed_time/count_sub;
            count_sub = count_sub + 1;
            time_remaining = (nSub-iSub)*time_per_sub;
            hms = fix(mod((time_remaining), [0, 3600, 60]) ./ [3600, 60, 1]);
            try
                waitbar(iSub/nSub, w, [num2str(iSub/nSub*100,'%.f') '%, ' num2str(hms(1),'%02.f') ':' num2str(hms(2),'%02.f') ':' num2str(hms(3),'%02.f') ' [hr:min:sec] remaining']);
            end
        end

    % ------------------------ Parallel Computing -------------------------
    case 1
        
        % Create waitbar
        try   % Waitbar for MATLAB R2017a and higher
            D = parallel.pool.DataQueue;            
            w = waitbar(0,'Please wait...','Name','Extract and correlate beta series','Tag','tmfc_waitbar');
            afterEach(D, @tmfc_parfor_waitbar);    
            tmfc_parfor_waitbar(w,nSub,1);     
        catch % No waitbar for MATLAB R2016b and earlier
            D = [];
            opts = struct('WindowStyle','non-modal','Interpreter','tex');
            w = warndlg({'\fontsize{12}Sorry, waitbar progress update is not available for parallel computations in MATLAB R2016b and earlier.',[],...
                'Please wait until all computations are completed.',[],...
                'If you want to interrupt computations:',...
                '   1) Do not close this window;',...
                '   2) Select MATLAB main window;',...
                '   3) Press Ctrl+C.'},'Please wait...',opts);
        end
        
        cleanupObj = onCleanup(@unfreeze_after_ctrl_c);
        
        try
            parpool;
            figure(findobj('Tag','TMFC_GUI'));
        end

        parfor iSub = 1:nSub
            tmfc_extract_betas_after_FIR(tmfc,ROI_set_number,ROIs,nROI,nCond,cond_list,XYZ,iXYZ,hdr,iSub);
            sub_check(iSub) = 1;

            % Update waitbar 
            try
                send(D,[]); 
            end
        end    
end

% Default contrasts info
for iCond = 1:nCond
    contrasts(iCond).title = cond_list(iCond).file_name;
    contrasts(iCond).weights = zeros(1,nCond);
    contrasts(iCond).weights(1,iCond) = 1;
end

% Close waitbar
try
    delete(w)
end

function unfreeze_after_ctrl_c()    
    try
        delete(findall(0,'type','figure','Tag', 'tmfc_waitbar'));
        GUI = guidata(findobj('Tag','TMFC_GUI')); 
        set([GUI.TMFC_GUI_B1, GUI.TMFC_GUI_B2, GUI.TMFC_GUI_B3, GUI.TMFC_GUI_B4,...
           GUI.TMFC_GUI_B5a, GUI.TMFC_GUI_B5b, GUI.TMFC_GUI_B6, GUI.TMFC_GUI_B7,...
           GUI.TMFC_GUI_B8, GUI.TMFC_GUI_B9, GUI.TMFC_GUI_B10, GUI.TMFC_GUI_B11,...
           GUI.TMFC_GUI_B12a,GUI.TMFC_GUI_B12b,GUI.TMFC_GUI_B13a,GUI.TMFC_GUI_B13b,...
           GUI.TMFC_GUI_B14a, GUI.TMFC_GUI_B14b], 'Enable', 'on');
    end
end
end

%% ========================================================================

% Extract and correlate betas
function tmfc_extract_betas_after_FIR(tmfc,ROI_set_number,ROIs,nROI,nCond,cond_list,XYZ,iXYZ,hdr,iSub)
    SPM = load(tmfc.subjects(iSub).path).SPM;

    % Load individual ROIs
    if isempty(ROIs)
        for iROI = 1:nROI
            ROIs(iROI).mask = spm_data_read(spm_data_hdr_read(tmfc.ROI_set(ROI_set_number).ROIs(iROI).path_masked(iSub).subjects),'xyz',XYZ);
            ROIs(iROI).mask(ROIs(iROI).mask == 0) = NaN;
        end
    end

    % Number of trials per condition
    nTrialCond = [];
    for jCond = 1:nCond
        nTrialCond(jCond) = length(SPM.Sess(cond_list(jCond).sess).U(cond_list(jCond).number).ons);
    end
    
    % Conditions of interest
    for jCond = 1:nCond

        % Extract average beta series from ROIs
        disp(['Extracting average beta series: Subject: ' num2str(iSub) ' || Condition: ' num2str(jCond)]);
        
        for kTrial = 1:nTrialCond(jCond)
            betas(kTrial,:) = spm_data_read(spm_data_hdr_read(fullfile(tmfc.project_path,'LSS_regression_after_FIR',tmfc.subjects(iSub).name,'Betas', ...
                ['Beta_' cond_list(jCond).file_name '_[Trial_' num2str(kTrial) '].nii'])),'xyz',XYZ);
        end

        for kROI = 1:nROI
            betas_masked = betas;
            betas_masked(:,isnan(ROIs(kROI).mask)) = []; 
            if strcmp(tmfc.ROI_set(ROI_set_number).BSC_after_FIR,'mean')
                beta_series(jCond).ROI_average(:,kROI) = mean(betas_masked,2);
            elseif strcmp(tmfc.ROI_set(ROI_set_number).BSC_after_FIR,'first_eigenvariate')
                betas_masked = betas_masked - mean(betas_masked,1);
                [m,n]   = size(betas_masked);
                if m > n
                    [v,s,v] = svd(betas_masked'*betas_masked);
                    s       = diag(s);
                    v       = v(:,1);
                    u       = betas_masked*v/sqrt(s(1));
                else
                    [u,s,u] = svd(betas_masked*betas_masked');
                    s       = diag(s);
                    u       = u(:,1);
                    v       = betas_masked'*u/sqrt(s(1));
                end
                d       = sign(sum(v));
                u       = u*d;
                beta_series(jCond).ROI_average(:,kROI) = (u*sqrt(s(1)/n))';    
                clear betas_masked v s u d
            end
        end

        % ROI-to-ROI correlation
        if tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 2
            z_matrix = atanh(corr(beta_series(jCond).ROI_average));
            z_matrix(1:size(z_matrix,1)+1:end) = nan;     

            % Save BSC matrices
            save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','ROI_to_ROI', ...
                [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.mat']),'z_matrix');

            clear z_matrix
        end

        % Seed-to-voxel correlation
        if tmfc.defaults.analysis == 1 || tmfc.defaults.analysis == 3
            for kROI = 1:nROI
                BSC_image(kROI).z_value = atanh(corr(beta_series(jCond).ROI_average(:,kROI),betas));
            end

            % Save BSC images
            for kROI = 1:nROI
                hdr.fname = fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR', ...
                    'Seed_to_voxel',tmfc.ROI_set(ROI_set_number).ROIs(kROI).name, ...
                    [tmfc.subjects(iSub).name '_Contrast_' num2str(jCond,'%04.f') '_' cond_list(jCond).file_name '.nii']);
                hdr.descrip = ['z-value map: ' cond_list(jCond).file_name];    
                image = NaN(SPM.xVol.DIM');
                image(iXYZ) = BSC_image(kROI).z_value;
                spm_write_vol(hdr,image);
            end

            clear BSC_image
        end

        clear betas  
    end

    % Save mean beta-series
    save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'BSC_LSS_after_FIR','Beta_series', ...
        [tmfc.subjects(iSub).name '_beta_series.mat']),'beta_series');
end

% Waitbar for parallel mode
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




           