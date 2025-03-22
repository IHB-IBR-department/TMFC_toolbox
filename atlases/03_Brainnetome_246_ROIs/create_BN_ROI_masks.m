clear

path = fileparts(which('create_BN_ROI_masks.m'));
load(fullfile(path,'BN_atlas.mat'));
mkdir(fullfile(path,'ROI_masks_2mm'));

BN{1,1} = fullfile(path,'BN_Atlas_246_2mm.nii');

for iROI = 1:246
    ROI_mask = fullfile(path,'ROI_masks_2mm',[num2str(iROI,'%03.f') '_' ROIs(iROI).Label '.nii']);
    spm_imcalc(BN,ROI_mask,['i1==' num2str(ROIs(iROI).BN_ID)],{0,0,1,2});
    clear ROI_mask
end