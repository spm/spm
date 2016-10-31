function [xVi, mask] = spm_est_non_sphericity(SPM)
% Non-sphericity estimation using ReML
% FORMAT [xVi, mask] = spm_est_non_sphericity(SPM)
% 
% Required fields of SPM structure (see spm_spm):
% SPM.xY.VY  - nScan x 1 struct array of file handles
% SPM.xX     - structure containing design matrix information
% SPM.xX.W   - optional whitening/weighting matrix
% SPM.xVi    - structure describing intrinsic non-sphericity
% SPM.xM     - structure containing masking information
% 
% Return xVi from SPM.xVi with extra fields:
% xVi.V      - estimated non-sphericity, trace(V) = rank(V)
% xVi.h      - hyperparameters  xVi.V = xVi.h(1)*xVi.Vi{1} + ...
% xVi.Cy     - spatially whitened <Y*Y'> (used by ReML to estimate h)
% 
% mask       - logical array of voxels within analysis mask
%__________________________________________________________________________
%
% In a first pass, voxels over which non-sphericity will be estimated are
% selected using an 'effects of interest' F-contrast (can be specified in 
% SPM.xVi.Fcontrast) and critical threshold taken from SPM defaults
% stats.<modality>.UFp.
% The sample covariance matrix (xVi.Cy) is then estimated by pooling over
% these voxels, assuming V is constant over them.
% Finally, SPM will invoke ReML to estimate hyperparameters (xVi.h) of an
% array of non-sphericity components (xVi.Vi), providing a high precise
% estimate of the non-sphericity matrix (xVi.V).
%__________________________________________________________________________
% Copyright (C) 1994-2016 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Guillaume Flandin
% $Id: spm_est_non_sphericity.m 6913 2016-10-31 10:12:27Z guillaume $


SVNid = '$Rev: 6913 $';

%-Say hello
%--------------------------------------------------------------------------
spm('FnBanner',mfilename,SVNid);


%-Get data, design, mask and variance components details
%--------------------------------------------------------------------------
VY           = SPM.xY.VY;
M            = VY(1).mat;
DIM          = VY(1).dim;
YNaNrep      = spm_type(VY(1).dt(1),'nanrep');

xX           = SPM.xX;
nScan        = size(xX.X,1);
if ~isfield(xX,'W')
    xX.W     = speye(nScan,nScan);
end

xM           = SPM.xM;

xVi          = SPM.xVi;

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


%==========================================================================
%-         P O O L E D   V A R I A N C E   E S T I M A T I O N
%==========================================================================

%-Get explicit mask(s)
%==========================================================================
mask = true(DIM);
for i = 1:numel(xM.VM)
    if ~(isfield(SPM,'xVol') && isfield(SPM.xVol,'G'))
        %-Assume it fits entirely in memory
        C = spm_bsplinc(xM.VM(i), [0 0 0 0 0 0]');
        v = true(DIM);
        [x1,x2] = ndgrid(1:DIM(1),1:DIM(2));
        for x3 = 1:DIM(3)
            M2  = inv(M\xM.VM(i).mat);
            y1 = M2(1,1)*x1+M2(1,2)*x2+(M2(1,3)*x3+M2(1,4));
            y2 = M2(2,1)*x1+M2(2,2)*x2+(M2(2,3)*x3+M2(2,4));
            y3 = M2(3,1)*x1+M2(3,2)*x2+(M2(3,3)*x3+M2(3,4));
            v(:,:,x3) = spm_bsplins(C, y1,y2,y3, [0 0 0 0 0 0]') > 0;
        end
        mask = mask & v;
        clear C v x1 x2 x3 M2 y1 y2 y3
    else
        if spm_mesh_detect(xM.VM(i))
            v = xM.VM(i).private.cdata() > 0;
        else
            v = spm_mesh_project(gifti(SPM.xVol.G), xM.VM(i)) > 0;
        end
        mask = mask & v(:);
        clear v
    end
end

Cy        = 0;                                  %-<Y*Y'> spatially whitened
Cm        = mask;                               %-mask of selected voxels

%-Split data into chunks
%==========================================================================
chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

spm_progress_bar('Init',nbchunks,'Hyperparameter estimation','Chunks');

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
    
    %-Report progress
    %======================================================================
    if i > 1, fprintf(repmat(sprintf('\b'),1,72)); end                  %-# 
    fprintf('%-40s: %30s', sprintf('Chunk %3d/%-3d',i,nbchunks),...
                           '...processing');                            %-#
                       
    %-Get data & construct analysis mask
    %----------------------------------------------------------------------
    Y       = zeros(nScan,numel(chunk));
    cmask   = mask(chunk);
    for j=1:nScan
        if ~any(cmask), break, end                 %-Break if empty mask
        
        Y(j,cmask) = spm_data_read(VY(j),chunk(cmask));%-Read chunk of data
        
        cmask(cmask) = Y(j,cmask) > xM.TH(j);      %-Threshold (& NaN) mask
        if xM.I && ~YNaNrep && xM.TH(j) < 0        %-Use implicit mask
            cmask(cmask) = abs(Y(j,cmask)) > eps;
        end
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));        %-Mask constant data
    mask(chunk)  = cmask;
    Cm(chunk)    = cmask;
    if ~any(cmask), continue, end
    Y       = Y(:,cmask);                          %-Data within mask

    %-Remove filter confounds
    %----------------------------------------------------------------------
    KWY      = spm_filter(xX.K,xX.W*Y);
    
    %-Ordinary Least Squares estimation
    %----------------------------------------------------------------------
    beta    = xX.pKX*KWY;                          %-Parameter estimates
    if any(cmask)
        res = spm_sp('r',xX.xKXs,KWY);             %-Residuals
    else
        res = zeros(nScan,0);
    end
    ResSS   = sum(res.^2);                         %-Residual SSQ
    clear res
    
    %-F-threshold & accumulate spatially whitened Y*Y'
    %----------------------------------------------------------------------
    j       = sum((Hsqr*beta).^2,1)/trMV > UF*ResSS/trRV;
    Cm(chunk(cmask)) = j;
    q       = nnz(j);
    if q
        q   = spdiags(sqrt(trRV./ResSS(j)'),0,q,q);
        Y   = Y(:,j)*q;
        Cy  = Cy + Y*Y';
    end
    
    %-Report progress
    %======================================================================
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');             %-#
    spm_progress_bar('Set',i);
end

fprintf('\n');                                                          %-#
spm_progress_bar('Clear');

s = nnz(Cm);                                    %-Number of selected voxels
if ~s
    error('Please check your data: There are no significant voxels.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Estimate sample covariance matrix Y*Y' from voxels in Cm
% Cy        = 0;
% nbchunks  = ceil(s / chunksize);
% chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),s+1);
% vx        = find(Cm);
% for i=1:nbchunks
%     chunk = chunks(i):chunks(i+1)-1;
%     n     = numel(chunk);
%     Y     = zeros(nScan,n);
%     for j=1:nScan
%         Y(j,:) = spm_data_read(VY(j),vx(chunk));
%     end
%     % store ResSS above to prevent recomputing it here:
%     %--------------------------------------------------
%     KY    = spm_filter(xX.K,Y);
%     ResSS = spm_sp('r',xX.xKXs,KY);
%     ResSS = sum(ResSS.^2); 
%     %--------------------------------------------------
%     q     = spdiags(sqrt(trRV./ResSS'),0,n,n);
%     Y     = Y*q;
%     Cy    = Cy + Y*Y';
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cy = Cy / s;                                    %-Sample covariance matrix


%==========================================================================
%-                     R e M L   E S T I M A T I O N
%==========================================================================

%-ReML estimate of residual correlations through hyperparameters (h)
%--------------------------------------------------------------------------
str    = sprintf('Temporal non-sphericity (%d voxels)',s);
fprintf('%-40s: %30s\n',str,'...ReML estimation');                      %-#

%-ReML for separable designs and covariance components
%--------------------------------------------------------------------------
if isstruct(xX.K)
    m     = length(xVi.Vi);
    h     = zeros(m,1);
    V     = sparse(nScan,nScan);
    for i = 1:length(xX.K)
        
        % extract blocks from bases
        %------------------------------------------------------------------
        q     = xX.K(i).row;
        p     = [];
        Qp    = {};
        for j = 1:m
            if nnz(xVi.Vi{j}(q,q))
                Qp{end + 1} = xVi.Vi{j}(q,q);
                p           = [p j];
            end
        end
        
        % design space for ReML (with confounds in filter)
        %------------------------------------------------------------------
        Xp     = xX.X(q,:);
        try
            Xp = [Xp xX.K(i).X0];
        end
        
        % ReML
        %------------------------------------------------------------------
        fprintf('%-30s\n',sprintf('  ReML Block %i',i));
        [Vp,hp] = spm_reml(Cy(q,q),Xp,Qp);
        V(q,q)  = V(q,q) + Vp;
        h(p)    = hp;
    end
else
    [V,h] = spm_reml(Cy,xX.X,xVi.Vi);
end

%-Normalize non-sphericity and save hyperparameters
%--------------------------------------------------------------------------
V      = V*nScan/trace(V);
xVi.h  = h;
xVi.V  = V;                                     % Save non-sphericity xVi.V
xVi.Cy = Cy;                                    % spatially whitened <Y*Y'>
