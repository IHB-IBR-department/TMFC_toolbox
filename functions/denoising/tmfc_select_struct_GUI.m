function [struct_paths] =  tmfc_select_struct_GUI(subject_paths)

warning('backtrace','off')
struct_paths = [];
txt_filter = '';
struct_subfolder = '';
no_file = '';

ST_MW = figure('Name','Select structual images','NumberTitle','off','Units','normalized','Position',[0.36 0.25 0.35 0.575],'MenuBar','none','ToolBar','none','color','w','CloseRequestFcn',@ST_MW_exit);

ST_S1_str = {'Structural images must be in native space (not normalized).'};
ST_txt_1 = uicontrol(ST_MW,'Style','text','String',ST_S1_str,'Units','normalized','Position',[0.03 0.90 0.95 0.05],'fontunits','normalized','FontSize',0.48,'HorizontalAlignment','left','backgroundcolor','w');

ST_MW_S1 = uicontrol(ST_MW,'Style','pushbutton','String','Select ''struct'' subfolder for the Subject No 1','Units','normalized','Position',[0.025 0.80 0.555 0.080],'FontUnits','normalized','FontSize',0.295,'callback',@select_subfolder);
ST_MW_S1_panel = uipanel(ST_MW,'Units','normalized','Position',[0.64 0.800 0.335 0.078],'HighLightColor',[0.78 0.78 0.78],'BorderType','line','backgroundcolor','w');
ST_MW_S1_txt = uicontrol(ST_MW,'Style','text','String','Not Selected','ForegroundColor','red','Units','normalized','Position',[0.660 0.813 0.300 0.040],'FontUnits','normalized','FontSize',0.60,'backgroundcolor','w');

ST_MW_S2 = uicontrol(ST_MW,'Style','pushbutton','String','Enter a text filter unique for structural images:','Units','normalized','Position',[0.025 0.689 0.555 0.080],'FontUnits','normalized','FontSize',0.295,'callback',@apply_filter);
ST_MW_S2_E = uicontrol(ST_MW,'Style','Edit','String','T1*.nii','Units','normalized','Position',[0.64 0.689 0.335 0.080],'FontUnits','normalized','FontSize',0.32,'backgroundcolor','w');

ST_MW_LB1 = uicontrol(ST_MW,'Style','listbox','String','','Max',100000,'Units','normalized','Position',[0.025 0.2400 0.95 0.420],'FontUnits','points','FontSize',10,'Value',[]);

ST_MW_S3 = uicontrol(ST_MW,'Style','pushbutton','String','Select all structural images manually','Units','normalized','Position',[0.0245 0.130 0.45 0.080],'FontUnits','normalized','FontSize',0.30,'callback',@add_images);
ST_MW_S4 = uicontrol(ST_MW,'Style','pushbutton','String','Clear all','Units','normalized','Position',[0.524 0.130 0.451 0.080],'FontUnits','normalized','FontSize',0.30,'callback',@clear_paths);

ST_MW_OK = uicontrol(ST_MW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.345 0.025 0.3 0.080],'FontUnits','normalized','FontSize',0.30,'callback',@export_paths);
movegui(ST_MW,'center');


% Close GUI 
function ST_MW_exit(~,~)
    struct_paths = [];
    warning('Structural images are not selected.');
    delete(ST_MW);
end

% Select 'struct' subfolder for the first subject
function select_subfolder(~,~)
    set(ST_MW_S1_txt,'String','Not selected','ForegroundColor','red');
    struct_subfolder = spm_select(1,'dir','Select ''struct'' subfolder',{},subject_paths{1},'..');
    if ~isempty(struct_subfolder)
       set(ST_MW_LB1,'String','');
       set(ST_MW_S1_txt,'String','Selected','ForegroundColor',[0.219 0.341 0.137])
    else
       warning('Please select ''structural'' subfolder for 1st subject.'); 
    end
end

% Apply text filter and generate paths 
function apply_filter(~,~)
    if ~isempty(struct_subfolder)
        txt_filter = get(ST_MW_S2_E,'String');
        txt_filter = strrep(txt_filter,' ','');
        if ~isempty(txt_filter) 
            struct_paths = [];
            no_file = [];     

            % Prepare subject paths
            subject_paths = strtrim(subject_paths(:,:));

            % Stuctural subfolder for 1st subject
            struct_subfolder = strrep(struct_subfolder,strtrim(subject_paths(1,:)),'');

            % Prepare struct paths for all subjects
            for iSub = 1:size(subject_paths,1)
                struct_file = dir(fullfile(subject_paths{iSub},struct_subfolder{1},txt_filter));
                if ~isempty(struct_file)
                    struct_paths = vertcat(struct_paths,fullfile(subject_paths{iSub},struct_subfolder,struct_file(1).name));
                else
                    no_file = vertcat(no_file,subject_paths{iSub});
                end
            end
            set(ST_MW_LB1,'String',struct_paths);
            clear struct_file
        else
            warning('Filter is empty or invalid, please re-enter.');
        end
    else
        warning('Please select ''structural'' subfolder for 1st subject before applying text filter.');
    end
end

% Clear struct paths
function clear_paths(~,~)
    if isempty(struct_paths)
        disp('No images to remove.');
    else
        struct_path = {};
        struct_paths = {};
        no_file = {};
        set(ST_MW_LB1,'String',struct_path,'Value',[]);
        disp('All images have been removed.');
    end
end

% Add all structural images manually
function add_images(~,~)
    struct_paths = '';
    no_file = [];
    struct_paths = spm_select(inf,'any','Select structural images for all subjects',{},subject_paths{1},'T1.*\.nii$');
    struct_paths = cellstr(struct_paths);
    if isempty(struct_paths)
        warning('Structural images are not selected, please try again.');
    end
    set(ST_MW_LB1,'String',struct_paths);
end

% Export struct paths
function export_paths(~,~)
    if ~isempty(struct_paths)
        if size(struct_paths,1) == size(subject_paths,1) 
            disp('Structural images have been selected.');
            delete(ST_MW);
        elseif size(struct_paths,1) > size(subject_paths,1) 
            warning('The number of structural images should be equal to the number of selected SPM.mat files. Please try again.');
        else
            if isempty(no_file)
                no_file = [];
                for iSub =  size(struct_paths,1)+1:size(subject_paths,1)
                    no_file = vertcat(no_file,subject_paths{iSub});
                end
            end
            missing_file_GUI(no_file);
         end
    else 
        warning('No images selected, please try again.');
    end
end

% Warning window: missing images
function missing_file_GUI(file_missing)
    ST_WW = figure('Name','Select subjects','NumberTitle','off','Units','normalized','Position',[0.32 0.30 0.35 0.28],'color','w','MenuBar','none','ToolBar','none','WindowStyle','Modal');
    ST_WW_LB = uicontrol(ST_WW,'Style','listbox','String',file_missing,'Max',inf,'Units','normalized','Position',[0.032 0.250 0.940 0.520],'FontUnits','points','FontSize',10,'Value',[]);
    ST_WW_S1 = uicontrol(ST_WW,'Style','text','String','Warning, structural images are missing for the following subjects:','Units','normalized','Position',[0.15 0.820 0.720 0.095],'FontUnits','normalized','FontSize',0.5,'backgroundcolor','w');
    ST_WW_close = uicontrol(ST_WW,'Style','pushbutton','String','OK','Units','normalized','Position',[0.415 0.06 0.180 0.120],'FontUnits','normalized','FontSize',0.30,'callback',@close_SS_WW);
    movegui(ST_WW,'center');
    uiwait();        
    function close_SS_WW(~,~)
        close(ST_WW);
    end
end

uiwait();
end