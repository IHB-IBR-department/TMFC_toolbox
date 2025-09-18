function min_max_ax = tmfc_axis(matrix,type)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Defines colorbar limits. The color scale will be adjusted based on the
% maximum absolute value and is guaranteed to be symmetric about zero.
%
% min_max_ax = tmfc_axis(matrix,type)
%
% matrix      - connectivity matrix
% type        - defines the output:
%               1: [-max_ax, max_ax]
%               0: [-max_ax, 0, max_ax]
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

matrix(1:1+size(matrix,1):end) = NaN;
max_ax = max(abs(matrix),[],'all','omitnan');

if ~isfinite(max_ax) || max_ax <= 0
    max_ax = 1e-6;     
end

if max_ax >= 1e-4
    max_ax = round(max_ax*1e4)/1e4;
end

if type == 1
    min_max_ax = [-max_ax, max_ax];
else   
    min_max_ax = [-max_ax, 0, max_ax]; 
end 
end
