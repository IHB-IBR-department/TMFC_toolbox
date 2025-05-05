function [funct_paths] = tmfc_select_funct_GUI(SPM_paths,subject_paths)

warning('backtrace','off')
funct_paths = [];
nPrefix = 1;
unsmoothed_path = [];
no_files = [];

SF_MW = figure('Name','Select unsmoothed functional images','NumberTitle','off','Units','normalized','Position',[0.325 0.202 0.35 0.575],'MenuBar','none','ToolBar','none','color','w','CloseRequestFcn',@SF_MW_exit);

SF_S1_str = {'Functional images must be realigned, normalized and unsmoothed.','',...
    'NOTE: Unsmoothed functional images are expected to be in the same folder as the functional images defined in the selected SPM.mat files and have a shorter prefix (e.g. ''war'' instead of ''swar'').'};

SF_txt = uicontrol(SF_MW,'Style','text','String',SF_S1_str,'Units','normalized','Position',[0.025 0.82 0.95 0.14],'fontunits','normalized','FontSize',0.16,'HorizontalAlignment','left','backgroundcolor','w');
SF_MW_B1 = uicontrol(SF_MW,'Style','pushbutton','String','Number of letters to remove from prefix:','Units','normalized','Position',[0.025 0.74 0.555 0.080],'FontUnits','normalized','FontSize',0.295,'callback',@prefix_filter);
SF_MW_B1_E = uicontrol(SF_MW,'Style','Edit','String',nPrefix,'Units','normalized','Position',[0.64 0.74 0.335 0.080],'FontUnits','normalized','FontSize',0.32,'backgroundcolor','w');

SF_MT_LB1_txt = uicontrol(SF_MW,'Style','text','String','The first selected functional image for each subject:','Units','normalized','Position',[0.025 0.67 0.95 0.040],'fontunits','normalized','FontSize',0.62,'HorizontalAlignment','center','backgroundcolor','w');
SF_MW_LB1 = uicontrol(SF_MW,'Style','listbox','String','','Max',100000,'Units','normalized','Position',[0.025 0.12 0.95 0.54],'FontUnits','points','FontSize',10,'Value',[]);
SF_MW_OK = uicontrol(SF_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.355 0.03 0.3 0.069],'FontUnits','normalized','FontSize',0.295,'callback',@export_paths);
movegui(SF_MW,'center');

subject_paths = strtrim(subject_paths(:,:));

% Apply default prefix
for iSub = 1:length(SPM_paths)
     SPM = load(SPM_paths{iSub}).SPM;
     [orig_path,orig_file,orig_ext] =  fileparts(SPM.xY.VY(1).fname);
     unsmoothed_file = orig_file(nPrefix+1:end);
     unsmoothed_path{iSub,1} = fullfile(orig_path,[unsmoothed_file orig_ext ',' num2str(SPM.xY.VY(1).n(1))]);
end
clear unsmoothed_file orig_path orig_file orig_ext iSub;
if ~isempty(unsmoothed_path)
    set(SF_MW_LB1,'String',unsmoothed_path);
end

% Сlose GUI
function SF_MW_exit(~,~)
    funct_paths = [];
    warning('Functional images are not selected.');
    delete(SF_MW);
end

% Apply prefix filter 
function prefix_filter(~,~)     
    temp_prefix = str2double(get(SF_MW_B1_E,'String'));
    if isnan(temp_prefix)
        warning('Please enter a natural number for prefix.');
    elseif ~(temp_prefix >= 0 && floor(temp_prefix) == temp_prefix)
        warning('Please enter a natural number for prefix.');
    else
        nPrefix = temp_prefix;
    end
    for jSub = 1:length(SPM_paths)
         SPM = load(SPM_paths{jSub}).SPM;
         [orig_path,orig_file,orig_ext] =  fileparts(SPM.xY.VY(1).fname);
         unsmoothed_file = orig_file(nPrefix+1:end);
         unsmoothed_path{jSub,1} = fullfile(orig_path,[unsmoothed_file orig_ext ',' num2str(SPM.xY.VY(1).n(1))]);
    end
    clear unsmoothed_file orig_path orig_file orig_ext SPM    
    set(SF_MW_LB1,'String',unsmoothed_path);
end

% Check and export paths
function export_paths(~,~)
    no_files = [];
    if ~isempty(unsmoothed_path)
        for jSub = 1:length(SPM_paths)
             SPM = load(SPM_paths{jSub}).SPM;
             for kScan = 1:length(SPM.xY.VY)
                  [orig_path,orig_file,orig_ext] =  fileparts(SPM.xY.VY(kScan).fname);
                  unsmoothed_file = orig_file(nPrefix+1:end);
                  funct_paths(jSub).fname{kScan,1} = fullfile(orig_path,[unsmoothed_file orig_ext ',' num2str(SPM.xY.VY(kScan).n(1))]);
                  funct_paths2(jSub).fname{kScan,1} = fullfile(orig_path,[unsmoothed_file orig_ext]);
             end
        end
    else
        warning('Functional images are not selected, please try again.');
    end
    
    disp('Checking realigned and unsmoothed functional images...')
    sub_check = zeros(1,length(funct_paths));
    for jSub = 1:length(funct_paths)
        image_check = zeros(1,length(funct_paths(jSub).fname));
        for kImage = 1:length(funct_paths(jSub).fname)
            if exist(funct_paths2(jSub).fname{kImage,1},'file')
               image_check(kImage) = 1; 
            end
        end
        if ~any(image_check==0)
            sub_check(jSub) = 1;
        end
    end
    
    if any(sub_check==0)
        missing_idx = find(sub_check==0);
        no_files = vertcat(no_files,subject_paths(missing_idx,:));
        missing_images_GUI(no_files);
    else
        delete(SF_MW);
        disp('Functional images selected.')
    end
end

% Warning window: missing images
function missing_images_GUI(no_files)
    SF_WW = figure('Name','Select subjects','NumberTitle','off','Units','normalized','Position',[0.32 0.30 0.35 0.28],'color','w','MenuBar','none','ToolBar','none','WindowStyle','Modal');
    SF_WW_LB = uicontrol(SF_WW,'Style','listbox','String',no_files,'Max',inf,'Units','normalized','Position',[0.032 0.250 0.940 0.520],'FontUnits','points','FontSize',10,'Value',[]);
    SF_WW_S1 = uicontrol(SF_WW,'Style','text','String','Warning, functional images are missing for the following subjects:','Units','normalized','Position',[0.15 0.820 0.720 0.095],'FontUnits','normalized','FontSize',0.5,'backgroundcolor','w');
    SF_WW_close = uicontrol(SF_WW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.415 0.06 0.180 0.120],'FontUnits','normalized','FontSize',0.30,'callback',@close_SF_WW);
    movegui(SF_WW,'center');
    uiwait();        
    function close_SF_WW(~,~)
        close(SF_WW);
    end
end

uiwait();

end