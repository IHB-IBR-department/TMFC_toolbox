function paths = tmfc_save_nii(timeseries,outputdir,filename,options)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
% 
% Creates NIfTI images for selected time-series.
%
% FORMAT paths = tmfc_save_nii(time_series,output_dir,filename)
% FORMAT paths = tmfc_save_nii(time_series,output_dir,filename,options)
%
% INPUTS:
% time_series   - Time series matrix (time x ROIs) for single subject
% output_dir    - Path to save *.nii files
% filename      - Name for *.nii files
%
% options.type  - '3D' (default): Save *.nii images as separate 3D files
%               - '4D': Save all *.nii images into single 4D file  
% option.format - 'float32' (default)
%               - 'float64'
%               - 'uint8'
%               - 'int16'
%               - 'int32'
%
% OUTPUTS:
% paths - Full paths to *.nii files.
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

if nargin < 4
    options.type = '3D';
    options.format = 'float32';
end

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
            for j = 1:size(timeseries,2)
                timeseries_4D(j,1,1,i) = timeseries(i,j);
            end
            paths(i).image = fullfile(outputdir,[filename '.nii,' num2str(i)]);
        end
        nii_image = make_nii(timeseries_4D,[],[],datatype);        
        save_nii(nii_image,fullfile(outputdir,[filename '.nii']));
end
end




