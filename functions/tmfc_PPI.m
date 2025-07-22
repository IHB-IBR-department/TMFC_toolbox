function [sub_check] = tmfc_PPI(tmfc,ROI_set_number,start_sub)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Calculates psycho-physiological interactions (PPIs).
%
% FORMAT [sub_check] = tmfc_PPI_after_FIR(tmfc)
% Run a function starting from the first subject in the list.
%
%   tmfc.subjects.path            - Paths to individual SPM.mat files
%   tmfc.subjects.name            - Subject names within the TMFC project
%                           ('Subject_XXXX' naming will be used by default)
%   tmfc.project_path             - Path where all results will be saved
%   tmfc.defaults.parallel        - 0 or 1 (sequential/parallel computing)
%
%   tmfc.ROI_set                  - List of selected ROIs
%   tmfc.ROI_set.PPI_centering    - Apply mean centering of psychological
%                                   regressor (PSY) prior to deconvolution:
%                                   'with_mean_centering' (default)
%                                   or 'no_mean_centering'
%                                   (Di, Reynolds & Biswal, 2017; Masharipov et al., 2024)
%   tmfc.ROI_set.PPI_whitening    - Apply whitening inversion of the seed 
%                                   time series prior to the deconvolution
%                                   and PPI term calculation to avoid double
%                                   prewhitening (He et al., 2025):
%                                   'inverse' (default) or 'none'
%
%   tmfc.ROI_set.type             - Type of the ROI set
%   tmfc.ROI_set.set_name         - Name of the ROI set
%   tmfc.ROI_set.ROIs.name        - Name of the selected ROI
%   tmfc.ROI_set.ROIs.path_masked - Paths to the ROI images masked by group
%                                   mean binary mask 
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions        - List of conditions of interest
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.sess   - Session number (as specified in SPM.Sess)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.number - Condition number (as specified in SPM.Sess.U)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.pmod   - Parametric/Time modulator number (see SPM.Sess.U.P)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.name   - Condition name (as specified in SPM.Sess.U.name(kPmod))
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.file_name - Condition-specific file names:
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
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).sess   = 1;   
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).number = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).name = 'Cond_A'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).file_name = '[Sess_1]_[Cond_1]_[Cond_A]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).sess   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).name = 'Cond_B'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).file_name = '[Sess_1]_[Cond_1]_[Cond_B]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).number = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).name = 'Cond_A'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).file_name = '[Sess_2]_[Cond_1]_[Cond_A]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).name = 'Cond_B'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).file_name = '[Sess_2]_[Cond_2]_[Cond_B]';
%
% If GLMs contain parametric or time modulators, add the following fields:
% e.g. first modulator for fourth condition:
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).sess   = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).pmod   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).name = 'Cond_BxModulator1^1';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).file_name = '[Sess_2]_[Cond_2]_[Cond_BxModulator1^1]'; 
% e.g. second modulator for second condition:
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).sess   = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).number = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).pmod = 3; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).name = 'Cond_BxModulator2^1'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).file_name = '[Sess_2]_[Cond_2]_[Cond_BxModulator2^1]'; 
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
% FORMAT [sub_check] = tmfc_PPI(tmfc,ROI_set_number,start_sub)
% Run the function starting from a specific subject in the path list for
% the selected ROI set.
%
%   tmfc                   - As above
%   ROI_set_number         - Number of the ROI set in the tmfc structure
%   start_sub              - Subject number on the path list to start with
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
   start_sub = 1;
elseif nargin == 2
   start_sub = 1;
end

if ~isfield(tmfc.ROI_set(ROI_set_number),'PPI_centering')
    tmfc.ROI_set(ROI_set_number).PPI_centering = 'with_mean_centering';
elseif isempty(tmfc.ROI_set(ROI_set_number).PPI_centering)
    tmfc.ROI_set(ROI_set_number).PPI_centering = 'with_mean_centering';
end

if ~isfield(tmfc.ROI_set(ROI_set_number),'PPI_whitening')
    tmfc.ROI_set(ROI_set_number).PPI_whitening = 'inverse';
elseif isempty(tmfc.ROI_set(ROI_set_number).PPI_whitening)
    tmfc.ROI_set(ROI_set_number).PPI_whitening = 'inverse';
end

% Check subject names
if ~isfield(tmfc.subjects,'name')
    for iSub = 1:length(tmfc.subjects)
        tmfc.subjects(iSub).name = ['Subject_' num2str(iSub,'%04.f')];
    end
end

try
    main_GUI = guidata(findobj('Tag','TMFC_GUI'));                           
    set(main_GUI.TMFC_GUI_S4,'String', 'Updating...','ForegroundColor',[0.772, 0.353, 0.067]);       
end

nSub = length(tmfc.subjects);
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);
cond_list = tmfc.ROI_set(ROI_set_number).gPPI.conditions;
sub_check = zeros(1,nSub);
if start_sub > 1
    sub_check(1:start_sub) = 1;
end

% Initialize waitbar
w = waitbar(0,'Please wait...','Name','PPI regressors calculation','Tag','tmfc_waitbar');
start_time = tic;
count_sub = 1;
cleanupObj = onCleanup(@unfreeze_after_ctrl_c);

spm('defaults','fmri');
spm_jobman('initcfg');

for iSub = start_sub:nSub

    if isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'PPIs',tmfc.subjects(iSub).name))
        rmdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'PPIs',tmfc.subjects(iSub).name),'s');
    end

    if ~isdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'PPIs',tmfc.subjects(iSub).name))
        mkdir(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'PPIs',tmfc.subjects(iSub).name));
    end

    % Conditions of interest
    for jCond = 1:length(cond_list)
        switch tmfc.defaults.parallel
            case 0  % Sequential
                for kROI = 1:nROI
                    tmfc_PEB_PPI(tmfc,ROI_set_number,cond_list,iSub,jCond,kROI);
                end
            case 1  % Parallel
                try
                    parpool;
                    figure(findobj('Tag','TMFC_GUI'));
                end
                parfor kROI = 1:nROI
                    tmfc_PEB_PPI(tmfc,ROI_set_number,cond_list,iSub,jCond,kROI);
                end
        end
    end
    
    sub_check(iSub) = 1;

    % Update main TMFC GUI
    try  
        main_GUI = guidata(findobj('Tag','TMFC_GUI'));                        
        set(main_GUI.TMFC_GUI_S4,'String', strcat(num2str(iSub), '/', num2str(nSub), ' done'),'ForegroundColor',[0.219, 0.341, 0.137]);    
    end
    
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

% Close waitbar
try
    delete(w)
end

% -------------------------------------------------------------------------
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

%==========================================================================
function tmfc_PEB_PPI(tmfc,ROI_set_number,cond_list,iSub,jCond,kROI)
% Adapted from spm_peb_ppi.m
    
% Load SPM 
SPM = load(tmfc.subjects(iSub).path).SPM;

% Load VOI
VOI = fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'VOIs', ... 
      tmfc.subjects(iSub).name,['VOI_' tmfc.ROI_set(ROI_set_number).ROIs(kROI).name '_' num2str(cond_list(jCond).sess) '.mat']);
p   = load(deblank(VOI(1,:)),'xY');

xY(1) = p.xY;
Sess  = SPM.Sess(xY(1).Sess);

PPI.name = ['[' regexprep(tmfc.ROI_set(ROI_set_number).ROIs(kROI).name,' ','_') ']_' cond_list(jCond).file_name];

% Define Uu
% Matrix of input variables and contrast weights. This is an
% [n x 3] matrix. The first column indexes SPM.Sess.U(i). The
% second column indexes the name of the input or cause, see
% SPM.Sess.U(i).name{j}. The third column is the contrast
% weight. Unless there are parametric effects the second
% column will generally be a 1.
U.name = {};
U.u    = [];
U.w    = [];
try
    Uu = [cond_list(jCond).number cond_list(jCond).pmod 1];
catch
    Uu = [cond_list(jCond).number 1 1];
end
for i = 1:size(Uu,1)
    U.u           = [U.u Sess.U(Uu(i,1)).u(33:end,Uu(i,2))];
    U.name{end+1} = Sess.U(Uu(i,1)).name{Uu(i,2)};
    U.w           = [U.w Uu(i,3)];
end

% See spm_peb_ppi:

% Setup variables
%--------------------------------------------------------------------------
RT      = SPM.xY.RT;
dt      = SPM.xBF.dt;
NT      = round(RT/dt);
fMRI_T0 = SPM.xBF.T0;
N       = length(xY(1).u);
k       = 1:NT:N*NT;                       % microtime to scan time indices

% Setup other output variables
%--------------------------------------------------------------------------
PPI.xY = xY;
PPI.RT = RT;
PPI.dt = dt;  

% Create basis functions and hrf in scan time and microtime
%--------------------------------------------------------------------------
hrf = spm_hrf(dt);

% Create convolved explanatory {Hxb} variables in scan time
%--------------------------------------------------------------------------
xb  = spm_dctmtx(N*NT + 128,N);
Hxb = zeros(N,N);
for i = 1:N
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb = xb(129:end,:);

% Get confounds (in scan time) and constant term
%--------------------------------------------------------------------------
X0 = xY(1).X0;
M  = size(X0,2);

% Get response variable
%--------------------------------------------------------------------------
if strcmp(tmfc.ROI_set(ROI_set_number).PPI_whitening,'inverse') % Whitening inversion to avoid double whitening (see He et al., 2025)
    % Check if SPM.mat has concatenated sessions 
    % (if spm_fmri_concatenate.m sript was used)
    if size(SPM.nscan,2) == size(SPM.Sess,2) 
        W = SPM.xX.W(SPM.xX.K(xY.Sess).row,SPM.xX.K(xY.Sess).row);
        Y = inv(W)*xY.u;
    else
        Y = inv(SPM.xX.W)*xY.u;
    end
else
    Y = xY.u;
end

% Remove confounds and save Y in ouput structure
%--------------------------------------------------------------------------
Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
PPI.Y = Yc(:,1);
if size(Y,2) == 2
    PPI.P = Yc(:,2);
end

% Specify covariance components; assume neuronal response is white
% treating confounds as fixed effects
%--------------------------------------------------------------------------
Q = speye(N,N)*N/trace(Hxb'*Hxb);
Q = blkdiag(Q, speye(M,M)*1e6  );

% Get whitening matrix (NB: confounds have already been whitened)
%--------------------------------------------------------------------------
W = SPM.xX.W(Sess.row,Sess.row);

% Create structure for spm_PEB
%--------------------------------------------------------------------------
clear P
P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
P{1}.C = speye(N,N)/4;      % i.i.d assumptions
P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
P{2}.C = Q;

C  = spm_PEB(Y,P);
xn = xb*C{2}.E(1:N);
xn = spm_detrend(xn);

% Setup psychological variable from inputs and contrast weights
%----------------------------------------------------------------------
PSY = zeros(N*NT,1);
for i = 1:size(U.u,2)
    PSY = PSY + full(U.u(:,i) * U.w(i));
end

% Mean centering
if strcmp(tmfc.ROI_set(ROI_set_number).PPI_centering,'with_mean_centering') || strcmp(tmfc.ROI_set(ROI_set_number).PPI_centering,'mean_centering') 
    PSY = spm_detrend(PSY);
end

% Multiply psychological variable by neural signal
%----------------------------------------------------------------------
PSYxn = PSY.*xn;

% Convolve, convert to scan time, and account for slice timing shift
%----------------------------------------------------------------------
ppi = conv(PSYxn,hrf);
ppi = ppi((k-1) + fMRI_T0);

% Convolve psych effect, convert to scan time, and account for slice
% timing shift
%----------------------------------------------------------------------
PSYHRF = conv(PSY,hrf);
PSYHRF = PSYHRF((k-1) + fMRI_T0);

% Save psychological variables
%----------------------------------------------------------------------
PPI.psy = U;
PPI.P   = PSYHRF;
PPI.xn  = xn;
PPI.ppi = spm_detrend(ppi);

% Save PPI *.mat file
%----------------------------------------------------------------------
save(fullfile(tmfc.project_path,'ROI_sets',tmfc.ROI_set(ROI_set_number).set_name,'PPIs',tmfc.subjects(iSub).name, ...
    ['PPI_[' regexprep(tmfc.ROI_set(ROI_set_number).ROIs(kROI).name,' ','_') ']_' cond_list(jCond).file_name '.mat']),'PPI');
end
