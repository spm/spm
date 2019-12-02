function [SPM] = spm_spm_Bayes(SPM)
% Conditional parameter estimation of a General Linear Model
% FORMAT [SPM] = spm_spm_Bayes(SPM)
%__________________________________________________________________________
%
% spm_spm_Bayes returns to voxels identified by spm_spm (ML parameter
% estimation) to get conditional parameter estimates and ReML hyper-
% parameter estimates.  These estimates use prior covariances, on the
% parameters, from empirical Bayes.  These PEB prior variances come from
% the hierarchical model that obtains by considering voxels as providing a
% second level.  Put simply, the variance in parameters, over voxels, is
% used as a prior variance from the point of view of any one voxel. The
% error covariance hyperparameters are re-estimated in the light of these
% priors.  The approach adopted is essentially a fully Bayesian analysis at
% each voxel, using empirical Bayesian prior variance estimators over
% voxels.
%
% Each separable partition (i.e. session) is assigned its own
% hyperparameter but within session covariance components are lumped
% together, using their relative expectations over voxels.  This makes
% things much more computationally efficient and avoids inefficient
% voxel-specific multiple hyperparameter estimates.
%
% spm_spm_Bayes adds the following fields to SPM:
%
%                           ----------------
%
%
%   SPM.PPM.l      = session-specific hyperparameter means
%   SPM.PPM.Cb     = empirical prior parameter covariances
%   SPM.PPM.C      = conditional covariances of parameters
%   SPM.PPM.dC{i}  = dC/dl;
%   SPM.PPM.ddC{i} = ddC/dldl
%
% The derivatives are used to compute the conditional variance of various
% contrasts in spm_getSPM, using a Taylor expansion about the hyperparameter
% means.
%
%
%                           ----------------
%
%   SPM.VCbeta     - Handles of conditional parameter estimates
%   SPM.VHp        - Handles of hyperparameter estimates
%
%                           ----------------
%
% Cbeta_????.<ext>                     - conditional parameter images
% These are 32-bit (float) images of the conditional estimates. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.<ext>) are given
% value NaN.
%
%                           ----------------
%
% CHp_????.<ext>              - error covariance hyperparameter images
% This is a 64-bit (double) image of the ReML error variance estimate.
% for each separable partition (Session).  Voxels outside the analysis 
% mask are given value NaN.
%__________________________________________________________________________
% 
% For single subject fMRI analysis there is an alternative function
% using voxel-wise GLM-AR models that are spatially regularised
% using the VB framework. This is implemented using spm_spm_vb.m.
%__________________________________________________________________________
% Copyright (C) 2002-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spm_Bayes.m 7739 2019-12-02 14:00:18Z guillaume $


%-Say hello
%--------------------------------------------------------------------------
Finter = spm('FigName','Stats: Bayesian estimation...');

%-Select SPM.mat & change directory
%--------------------------------------------------------------------------
if ~nargin
    [Pf, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, return; end
    swd = spm_file(Pf,'fpath');
    load(fullfile(swd,'SPM.mat'))
    cd(swd)
end

try
    M    = SPM.xVol.M;
    DIM  = SPM.xVol.DIM;
    xdim = DIM(1); ydim = DIM(2); zdim = DIM(3);
    XYZ  = SPM.xVol.XYZ;
catch
    helpdlg({   'Please do a ML estimation first.',...
            'This identifies the voxels to analyse.'});
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    return
end


%==========================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%==========================================================================

%-Initialise output images
%==========================================================================
fprintf('%-40s: %30s','Output images','...initialising')                %-#

%-Initialise conditional estimate image files
%--------------------------------------------------------------------------
xX             = SPM.xX;
[nScan,nBeta]  = size(xX.X);
Vbeta(1:nBeta) = deal(struct(...
            'fname',   [],...
            'dim',     DIM',...
            'dt',      [spm_type('float32'), spm_platform('bigend')],...
            'mat',     M,...
            'pinfo',   [1 0 0]',...
            'descrip', ''));
for i = 1:nBeta
    Vbeta(i).fname   = [sprintf('Cbeta_%04d',i) spm_file_ext];
    Vbeta(i).descrip = sprintf('Cond. beta (%04d) - %s',i,xX.name{i});
    spm_unlink(Vbeta(i).fname)
end
Vbeta = spm_create_vol(Vbeta);

%-Initialise ReML hyperparameter image files
%--------------------------------------------------------------------------
try
    nHp       = length(SPM.nscan);
catch
    nHp       = nScan;
    SPM.nscan = nScan;
end

VHp(1:nHp)    = deal(struct(...
            'fname',   [],...
            'dim',     DIM',...
            'dt',      [spm_type('float64'), spm_platform('bigend')],...
            'mat',     M,...
            'pinfo',   [1 0 0]',...
            'descrip', ''));
for i = 1:nHp
    VHp(i).fname   = [sprintf('Hp_%04d',i) spm_file_ext];
    VHp(i).descrip = sprintf('Hyperparameter (%04d)',i);
    spm_unlink(VHp(i).fname)
end
VHp   = spm_create_vol(VHp);

fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...initialised')           %-#


%==========================================================================
% - A V E R A G E   S A M P L E   C O V A R I A N C E   M A T R I X
%==========================================================================
fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...estimating CY')         %-#
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

SPM.xVi.CY = CY;
clear CY EY


%==========================================================================
% - E M P I R I C A L  B A Y E S  F O R  P R I O R  V A R I A N C E
%==========================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...estimating priors')   %-#

% get row u{i} and column v{i}/v0{i} indices for separable designs
%--------------------------------------------------------------------------
s = nHp;
if isfield(SPM,'Sess')
    for i = 1:s
         u{i} = SPM.Sess(i).row;
         v{i} = SPM.Sess(i).col;
        v0{i} = xX.iB(i);
    end
else
     u{1} = [1:nScan];
     v{1} = [xX.iH xX.iC];
    v0{1} = [xX.iB xX.iG];
end

% cycle over separarable partitions
%--------------------------------------------------------------------------
for i = 1:s

    % Get design X and confounds X0
    %----------------------------------------------------------------------
    fprintf('%-30s\n',sprintf('  ReML Session %i',i));                  %-#
    X     = xX.X(u{i}, v{i});
    X0    = xX.X(u{i},v0{i});
    [m,n] = size(X);

    % add confound in 'filter'
    %----------------------------------------------------------------------
    if isstruct(xX.K)
        X0 = full([X0 xX.K(i).X0]);
    end

    % orthogonalize X w.r.t. X0
    %----------------------------------------------------------------------
    X     = X - X0*(pinv(X0)*X);

    % covariance components induced by parameter variations {Q}
    %----------------------------------------------------------------------
    for j = 1:n
        Q{j} = X*sparse(j,j,1,n,n)*X';
    end

    % covariance components induced by error non-sphericity {V}
    %----------------------------------------------------------------------
    Q{n + 1} = SPM.xVi.V(u{i},u{i});

    % ReML covariance component estimation
    %----------------------------------------------------------------------
    [C,h]   = spm_reml(SPM.xVi.CY,X0,Q);

    % check for negative variance components
    %----------------------------------------------------------------------
    h       = abs(h);

    % 2-level model for this partition using prior variances sP(i)
    % treat confounds as fixed (i.e. infinite prior variance)
    %----------------------------------------------------------------------
    n0      = size(X0,2);
    Cb      = blkdiag(diag(h(1:n)),speye(n0,n0)*1e8);
    P{1}.X  = [X X0];
    P{1}.C  = {SPM.xVi.V};
    P{2}.X  = sparse(size(P{1}.X,2),1);
    P{2}.C  = Cb;

    sP(i).P = P;
    sP(i).u = u{:};
    sP(i).v = [v{:} v0{:}];
end


%==========================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%==========================================================================

%-Cycle to avoid memory problems (plane by plane)
%==========================================================================
spm_progress_bar('Init',100,'Bayesian estimation','');
spm('Pointer','Watch')

%-maxMem is the maximum amount of data processed at a time (bytes)
%--------------------------------------------------------------------------
MAXMEM = spm_get_defaults('stats.maxmem');
blksz  = ceil(MAXMEM/8/nScan);
SHp    = 0;             % sum of hyperparameters
for  z = 1:zdim

    % current plane-specific parameters
    %----------------------------------------------------------------------
    U       = find(XYZ(3,:) == z);
    nbch    = ceil(length(U)/blksz);
    CrBl    = zeros(nBeta,length(U)); %-conditional parameter estimates
    CrHp    = zeros(nHp,  length(U)); %-ReML hyperparameter estimates
    for bch = 1:nbch                  %-loop over bunches of lines (planks)
        
        %-construct list of voxels in this block
        %------------------------------------------------------------------
        I     = [1:blksz] + (bch - 1)*blksz;
        I     = I(I <= length(U));
        xyz   = XYZ(:,U(I));
        nVox  = size(xyz,2);
        
        %-Get response variable
        %------------------------------------------------------------------
        Y     = spm_get_data(SPM.xY.VY,xyz);
        
        %-Conditional estimates (per partition, per voxel)
        %------------------------------------------------------------------
        beta  = zeros(nBeta,nVox);
        Hp    = zeros(nHp,  nVox);
        for j = 1:s
            P     = sP(j).P;
            u     = sP(j).u;
            v     = sP(j).v;
            for i = 1:nVox
                C         = spm_PEB(Y(u,i),P);
                beta(v,i) = C{2}.E(1:length(v));
                Hp(j,i)   = C{1}.h;
            end
        end
        
        %-Save for current plane in memory as we go along
        %------------------------------------------------------------------
        CrBl(:,I) = beta;
        CrHp(:,I) = Hp;
        SHp       = SHp + sum(Hp,2);
        
    end % (bch)


    %-write out plane data to image files
    %======================================================================

    %-Write conditional beta images
    %----------------------------------------------------------------------
    for i = 1:nBeta
    tmp       = sparse(XYZ(1,U),XYZ(2,U),CrBl(i,:),xdim,ydim);
    tmp(~tmp) = NaN;
    Vbeta(i)  = spm_write_plane(Vbeta(i),tmp,z);
    end

    %-Write hyperparameter images
    %----------------------------------------------------------------------
    for i = 1:nHp
    tmp       = sparse(XYZ(1,U),XYZ(2,U),CrHp(i,:),xdim,ydim);
    tmp(~tmp) = NaN;
    VHp(i)    = spm_write_plane(VHp(i),tmp,z);
    end


    %-Report progress
    %----------------------------------------------------------------------
    spm_progress_bar('Set',100*(z - 1)/zdim);


end % (for z = 1:zdim)
fprintf('\n')                                                           %-#
spm_progress_bar('Clear')

%==========================================================================
% - P O S T   E S T I M A T I O N
%==========================================================================

% Taylor expansion for conditional covariance
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Non-sphericity','...REML estimation')          %-#

% expansion point (mean hyperparameters)
%--------------------------------------------------------------------------
l     = SHp/SPM.xVol.S;

% change in conditional coavriance w.r.t. hyperparameters
%--------------------------------------------------------------------------
n     = size(xX.X,2);
PPM.l = l;
for i = 1:s
    PPM.dC{i}  = sparse(n,n);
    PPM.ddC{i} = sparse(n,n);
end
for i = 1:s

    P     = sP(i).P;
    u     = sP(i).u;
    v     = sP(i).v;

    % derivatives of conditional covariance w.r.t. hyperparameters
    %----------------------------------------------------------------------
    d     = P{1}.X'*inv(P{1}.C{1})*P{1}.X;
    Cby   = inv(d/l(i) + inv(P{2}.C));
    d     = d*Cby;
    dC    = Cby*d/(l(i)^2);
    ddC   = 2*(dC/(l(i)^2) - Cby/(l(i)^3))*d;

    % place in output structure
    %----------------------------------------------------------------------
    j               = 1:length(v);
    PPM.Cb(v,v)     = P{2}.C(j,j);
    PPM.Cby(v,v)    = Cby(j,j);
    PPM.dC{i}(v,v)  = dC(j,j);
    PPM.ddC{i}(v,v) = ddC(j,j);        
    
end

%-Save remaining results files and analysis parameters
%==========================================================================
fprintf('%-40s: %30s','Saving results','...writing')                    %-#

%-Save analysis parameters in SPM.mat file
%--------------------------------------------------------------------------
SPM.VCbeta = Vbeta;         % Filenames - parameters
SPM.VHp    = VHp;           % Filenames - hyperparameters
SPM.PPM    = PPM;           % PPM structure

fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end
save('SPM.mat', 'SPM', fmt);

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#


%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
