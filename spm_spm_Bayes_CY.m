function CY = spm_spm_Bayes_CY(SPM)
% Estimation of the average sample covariance of whole-brain data.
%
% SPM - Structure (see spm_spm.m)
% CY  - Matrix of dimension [P x P] where P is the number of volumes
%__________________________________________________________________________
%
% Normalisation (and where appropriate, temporal pre-whitening) are 
% performed prior to estimation. Voxels are included that exceed a liberal
% statistical threshold, by default p < 0.001 uncorrected for effects of
% interest. The resulting matrix is used in spm_spm_Bayes.m.
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2002-2024 Wellcome Centre for Human Neuroimaging

xX  = SPM.xX;
DIM = SPM.xVol.DIM;

CY    = 0;                                       % <(Y - <Y>) * (Y - <Y>)'>
EY    = 0;                                       % <Y>    for ReML
nScan = size(xX.X,1);
xVi   = SPM.xVi;

%-Compute Hsqr and F-threshold under i.i.d.
%--------------------------------------------------------------------------
xX.xKXs      = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));
xX.xKXs.X    = full(xX.xKXs.X);
xX.pKX       = spm_sp('x-',xX.xKXs);

if isfield(xVi,'Fcontrast')
    Fcname   = 'User-specified contrast';
    xCon     = spm_FcUtil('Set',Fcname,'F','c',xVi.Fcontrast,xX.xKXs);
else
    Fcname   = 'effects of interest';
    iX0      = [xX.iB xX.iG];
    xCon     = spm_FcUtil('Set',Fcname,'F','iX0',iX0,xX.xKXs);
end

if ~isempty(xCon(1).c)
    X1o      = spm_FcUtil('X1o', xCon(1),xX.xKXs);
    Hsqr     = spm_FcUtil('Hsqr',xCon(1),xX.xKXs);
    trMV     = spm_SpUtil('trMV',X1o);
else
    % Force all voxels to enter non-sphericity
    trMV     = 1;
    Hsqr     = Inf;
end
trRV         = spm_SpUtil('trRV',xX.xKXs);

%-Threshold for voxels entering non-sphericity estimates
%--------------------------------------------------------------------------
try
    modality = lower(spm_get_defaults('modality'));
    UFp      = spm_get_defaults(['stats.' modality '.ufp']);
catch
    UFp      = 0.001;
end
xVi.UFp      = UFp;
UF           = spm_invFcdf(1 - UFp,[trMV,trRV]);

%-Split data into chunks
%--------------------------------------------------------------------------
VY        = SPM.xY.VY;
mask      = logical(spm_read_vols(SPM.VM));

chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
                       
    %-Get data & construct analysis mask
    %----------------------------------------------------------------------
    Y       = zeros(nScan,numel(chunk));
    cmask   = mask(chunk);
    for j=1:nScan
        if ~any(cmask), break, end                    %-Break if empty mask
        Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));%-Read chunk of data
    end
    mask(chunk)  = cmask;
    if ~any(cmask), continue, end
    Y       = Y(:,cmask);                             %-Data within mask

    %-Remove filter confounds
    %----------------------------------------------------------------------
    KWY     = spm_filter(xX.K,xX.W*Y);
    
    %-Ordinary Least Squares estimation
    %----------------------------------------------------------------------
    beta    = xX.pKX*KWY;                             %-Parameter estimates
    if any(cmask)
        res = spm_sp('r',xX.xKXs,KWY);                %-Residuals
    else
        res = zeros(nScan,0);
    end
    ResSS   = sum(res.^2);                            %-Residual SSQ
    clear res
    
    %-F-threshold & accumulate spatially whitened Y*Y'
    %----------------------------------------------------------------------
    j       = sum((Hsqr*beta).^2,1)/trMV > UF*ResSS/trRV;
    if nnz(j)
        Y   = Y(:,j);
        CY  = CY + Y*Y';
        EY  = EY + sum(Y,2);
    end
    
end

%-average sample covariance and mean of Y (over voxels)
%--------------------------------------------------------------------------
S  = nnz(mask);
CY = CY/S;
EY = EY/S;
CY = CY - EY*EY';
