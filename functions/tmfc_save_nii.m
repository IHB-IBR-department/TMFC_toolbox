function paths = tmfc_save_nii(timeseries,outputdir,filename,options)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
% 
% Creates NIfTI images for selected time series.
%
% Each voxel represents one ROI, and each volume represents one time point:
%   - In 3D mode: each time point is saved as a separate NIfTI image 
%     (ROI values stored as voxels in a 1D image).
%   - In 4D mode: all time points are stored as volumes in a single NIfTI file.
%
% FORMAT paths = tmfc_save_nii(timeseries,outputdir,filename)
% FORMAT paths = tmfc_save_nii(timeseries,outputdir,filename,options)
%
% INPUTS:
% timeseries    - Time series matrix (time x ROIs) for single subject
% outputdir     - Path to save *.nii files
% filename      - Name for *.nii files
%
% options.type   - '3D' (default): Save *.nii images as separate 3D files
%                - '4D': Save all *.nii images into single 4D file  
% options.format - 'float32' (default)
%                - 'float64'
%                - 'uint8'
%                - 'int16'
%                - 'int32'
%
% OUTPUTS:
% paths - Full paths to *.nii files.
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin < 4, options = struct; end
if ~isfield(options,'type'),   options.type = '3D'; end
if ~isfield(options,'format'), options.format = 'float32'; end

if strcmp(options.format,'float32')
    datatype = 16;
elseif strcmp(options.format,'float64')
    datatype = 64;
elseif strcmp(options.format,'uint8')
    datatype = 2;
elseif strcmp(options.format,'int16')
    datatype = 4;
elseif strcmp(options.format,'int32')
    datatype = 8;
else
    error('Unsupported format: %s', options.format);    
end

if ~isdir(outputdir)
    mkdir(outputdir);
end

switch options.type 
    case '3D'
        for i = 1:size(timeseries,1)
            nii_image = [];
            nii_image = make_nii(timeseries(i,:)',[],[],datatype);
            paths(i).image = fullfile(outputdir,[filename '_' num2str(i,'%05.f') '.nii']);
            save_nii(nii_image,paths(i).image);
        end
    case '4D'
        for i = 1:size(timeseries,1)
            timeseries_4D = zeros(size(timeseries,2),1,1,size(timeseries,1));
            for j = 1:size(timeseries,2)
                timeseries_4D(j,1,1,i) = timeseries(i,j);
            end
            paths(i).image = fullfile(outputdir,[filename '.nii,' num2str(i)]);
        end
        nii_image = make_nii(timeseries_4D,[],[],datatype);        
        save_nii(nii_image,fullfile(outputdir,[filename '.nii']));
end
end




