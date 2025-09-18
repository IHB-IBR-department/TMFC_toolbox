function z = tmfc_zscore(x,flag,dim)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Standardize data to zero mean and unit variance (no NaN handling).
%
% Z = ZSCORE(X) returns a centered, scaled copy of X, where each column
% has a mean of 0 and a standard deviation of 1.
%
% Z = ZSCORE(X,FLAG,DIM) standardizes along dimension DIM using
% normalization FLAG for the standard deviation:
%     FLAG = 0 (default) → normalize by N-1 (sample std)
%     FLAG = 1           → normalize by N   (population std)
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin < 2 || isempty(flag), flag = 0; end
if nargin < 3 || isempty(dim)
    dim = find(size(x) ~= 1, 1); 
    if isempty(dim), dim = 1; end
end

try 
    %cTry built-in zscore first (requires  Statistics and Machine Learning Toolbox)
    z = zscore(x, flag, dim);
catch
    mu    = mean(x, dim);
    sigma = std(x, flag, dim);
    sigma(sigma == 0) = 1;  % avoid division by zero
    
    % Standardize
    z = bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma);
end
end
