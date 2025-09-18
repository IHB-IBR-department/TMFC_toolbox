function tmfc_progress(cmd, varargin)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Minimal, consistent waitbar for serial + parfor.
% Usage:
%   tmfc_progress('init', totalN, 'Model estimation');   % once
%   tmfc_progress('tick');                               % each iteration
%   tmfc_progress('done');                               % at the end
%
% For parfor: pair with a DataQueue and call tmfc_progress('tick') in afterEach.
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru


persistent w nTot kk t0
switch lower(cmd)
    case 'init'
        nTot  = varargin{1};               % total iterations
        ttl = varargin{2};                 % window title
        kk  = 0;                           % completed
        t0 = tic;
        if isempty(w) || ~ishandle(w)
            w = waitbar(0, 'Please wait...', 'Name', ttl, 'Tag', 'tmfc_waitbar');
        else
            waitbar(0, w, 'Please wait...');
            set(w, 'Name', ttl);
        end

    case 'tick'
        if isempty(nTot) || nTot<=0 || isempty(w) || ~ishandle(w), return; end
        kk = min(kk+1, nTot);
        frac    = kk / nTot;
        elapsed = toc(t0);
        % Stable ETA using simple average time per iteration:
        eta     = elapsed * (nTot - kk) / max(kk,1);

        msg = sprintf('%.0f%% — %s [hr:min:sec] remaining', 100*frac, tmfc_sec2hms(eta));
        try
            waitbar(frac, w, msg);
        catch
            w = [];
        end

    case 'done'
        if ~isempty(w) && ishandle(w), delete(w); end
        w = []; nTot = []; kk = []; t0 = [];

    otherwise
        % no-op
end
end

function s = tmfc_sec2hms(sec)
if ~isfinite(sec) || sec < 0, sec = 0; end
h = floor(sec/3600); sec = sec - 3600*h;
m = floor(sec/60);   sec = sec - 60*m;
s = sprintf('%02d:%02d:%02.0f', h, m, sec);
end
