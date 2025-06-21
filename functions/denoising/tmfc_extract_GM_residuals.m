function GM_res = tmfc_extract_GM_residuals(SPM,GM_path,c)

% This is a modification of the original 'spm_write_residuals' function.
%
% Extract residuals from GM mask.
% FORMAT y = spm_write_residuals(SPM,GM_path,Ic)
% SPM    - structure containing generic analysis details
% c      - contrast weights to adjust data (0:   no adjustment)
%                                          (1:   adjust for covariates)
%                                          (NaN: adjust for everything) 
%
% y      - matrix with GM residuals (time x voxel)  
%__________________________________________________________________________
% Copyright (C) 2012-2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_write_residuals.m 6656 2015-12-24 16:49:52Z guillaume $

cwd = pwd; 
cd(SPM.swd);

%-Compute and write residuals
%--------------------------------------------------------------------------
DIM = SPM.xY.VY(1).dim(1:min(numel(SPM.xY.VY(1).dim),3));
[nScan, nBeta] = size(SPM.xX.X);

%-Loop over chunks
%--------------------------------------------------------------------------
chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

GM_res = [];

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
    
    %-Get GM mask
    %----------------------------------------------------------------------
    GM_mask = spm_data_hdr_read(GM_path);
    m = spm_data_read(GM_mask,chunk) > 0;
    m = m(:)';
    
    %-Get raw data, whiten and filter
    %----------------------------------------------------------------------
    y = zeros(nScan,numel(chunk));
    for j=1:nScan
        y(j,:) = spm_data_read(SPM.xY.VY(j),chunk);
    end
    y(:,~m) = [];
    
    y = spm_filter(SPM.xX.K,SPM.xX.W*y);
    
    if c ~= 0
        
        %-Parameter estimates: beta = xX.pKX*xX.K*y
        %------------------------------------------------------------------
        beta = zeros(nBeta,numel(chunk));
        for j=1:nBeta
            beta(j,:) = spm_data_read(SPM.Vbeta(j),chunk);
        end
        beta(:,~m) = [];
        
        %-Subtract Y0 = XO*beta,  Y = Yc + Y0 + e
        %------------------------------------------------------------------
        if ~isnan(c)
            % Contrast weights for Omnibus F-contrast
            weights = tmfc_omnibus_F_contrast(SPM);
            sX = SPM.xX.xKXs;
            ukX1o = spm_SpUtil('+c->Tsp',sX,weights);
            hsqr = spm_sp('ox',spm_sp('set',ukX1o))' * spm_sp('cukx',sX);
            H = hsqr' * hsqr;
            Y0 = sX.X*(eye(spm_sp('size',sX,2)) - spm_sp('xpx-',sX)*H)*beta; % Y0 = spm_FcUtil('Y0',SPM.xCon(Ic),SPM.xX.xKXs,beta);
            y = y - Y0; 
        else
            y = y - SPM.xX.xKXs.X * beta;
        end 
    end

    GM_res = [GM_res y];
end

cd(cwd);

end

%==========================================================================
function weights = tmfc_omnibus_F_contrast(SPM)
kCond = 1;
for iSess = 1:length(SPM.Sess)
    for jCond = 1:length(SPM.Sess(iSess).U)
        for kPmod = 1:length(SPM.Sess(iSess).U(jCond).name)
            cond_list(kCond).sess = iSess;
            cond_list(kCond).number = jCond;
            cond_list(kCond).pmod = kPmod;
            kCond = kCond + 1;
        end
    end 
end

cond_col = [];
for iCond = 1:length(cond_list)
    FCi = [];
    FCi = SPM.Sess(cond_list(iCond).sess).Fc(cond_list(iCond).number).i; 
    try
        FCp = []; 
        FCp = SPM.Sess(cond_list(iCond).sess).Fc(cond_list(iCond).number).p; 
        FCi = FCi(FCp==cond_list(iCond).pmod);
    end
    cond_col = [cond_col SPM.Sess(cond_list(iCond).sess).col(FCi)];
end 
weights = zeros(length(cond_col),size(SPM.xX.X,2));
for iCond = 1:length(cond_col)
    weights(iCond,cond_col(iCond)) = 1;
end
weights = weights';
if ~spm_sp('isinspp',SPM.xX.xKXs,weights)
    weights = spm_sp('oPp:',SPM.xX.xKXs,weights);
end
end
