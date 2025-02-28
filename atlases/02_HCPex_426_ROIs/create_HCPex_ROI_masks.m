clear

path = fileparts(which('create_HCPex_ROI_masks.m'));
load(fullfile(path,'HCPex_sorted.mat'));
mkdir(fullfile(path,'ROI_masks_2mm'));

HCPex{1,1} = fullfile(path,'HCPex_2mm.nii');

for iROI = 1:426
    ROI_mask = fullfile(path,'ROI_masks_2mm',[num2str(iROI,'%03.f') '_' ROIs(iROI).Short_label '.nii']);
    spm_imcalc(HCPex,ROI_mask,['i1==' num2str(ROIs(iROI).HCPex_ID)],{0,0,1,2});
    clear ROI_mask
end