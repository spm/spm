function [M,h] = spm_maff(varargin)
% Affine registration to MNI space using mutual information
% FORMAT M = spm_maff
% FORMAT M = spm_maff(V)
% FORMAT M = spm_maff(V,opts)
% V    - image filename/handle
% opts -  a structure containing optional fields
%         M       - starting estimate
%         tpm     - filenames of belonging probability images
%         regtype - regularisation type, a string of either
%                   'mni'   - registration with MNI space
%                   'rigid' - rigid(ish)-body registration
%                   'subj'  - inter-subject registration
%                   'none'  - no regularisation
%         fudge   - regularisation fudge factor
%         samp    - approximate sampling distance (mm)
%
% M    - resulting affine transform
%
% FORMAT M = spm_maff(buf,x,b0,MG,MF,M,regtyp,ff)
% buf     - a structure array (one per plane) containing
%           msk - a mask (logical) showing included voxels
%           g   - data (unsigned byte) at each position
%           nm  - number of data points
% x       - cell array of {x1,x2,x3}, where x1 and x2 are
%           co-ordinates (from ndgrid), and x3 is a list of
%           slice numbers to use
% b0      - a cell array of belonging probability images
%           (see spm_load_priors.m).
% MG      - voxel-to-world transform of image
% MF      - voxel-to-world transform of belonging probability
%           images
% M       - starting estimates
% regtype - regularisation type (see above)
% ff      - a fudge factor (derived from the one above)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$


if nargin>2,
    M = affreg(varargin{:});
    return;
end;

if nargin==0
    V = spm_select(1,'image');
else
    V = varargin{1};
end;
if ischar(V), V = spm_vol(V); end;

B = strvcat(...
    fullfile(spm('Dir'),'tpm','grey.nii'),...
    fullfile(spm('Dir'),'tpm','white.nii'),...
    fullfile(spm('Dir'),'tpm','csf.nii'));

opts0 = struct('M',eye(4), 'tpm',B,'regtyp','mni','fudge',9,'samp',2);

if nargin < 2
    opts = opts0;
else
    opts = varargin{2};
    fnms = fieldnames(opts0);
    for i=1:length(fnms)
        if ~isfield(opts,fnms{i}), opts.(fnms{i}) = opts0.(fnms{i}); end;
    end;
end;

B  = spm_vol(opts.tpm);
b0 = spm_load_priors(B);

% Load the image
V         = spm_vol(V);
d         = V(1).dim(1:3);
vx        = sqrt(sum(V(1).mat(1:3,1:3).^2));
sk        = max([1 1 1],round(opts.samp*[1 1 1]./vx));
[x1,x2,o] = ndgrid(1:sk(1):d(1),1:sk(2):d(2),1);
x3        = 1:sk(3):d(3);
d         = [size(x1) length(x3)];
g         = zeros(d);
spm_progress_bar('Init',V.dim(3),'Loading volume','Planes loaded');
for i=1:d(3)
    g(:,:,i) = spm_sample_vol(V,x1,x2,o*x3(i),0);
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');

% Convert the image to unsigned bytes
[mn,mx] = spm_minmax(g);
for z=1:length(x3),
    gz         = g(:,:,z);
    buf(z).msk = gz>mn;
    buf(z).nm  = sum(buf(z).msk(:));
    gz         = double(gz(buf(z).msk));
    buf(z).g   = uint8(round(gz*(255/mx)));
end;
clear g

ff   = opts.fudge;
ff   = max(1,ff^3/prod(sk)/abs(det(V.mat(1:3,1:3))));
[M,h] = affreg(buf,{x1,x2,x3},b0,V.mat,B(1).mat,opts.M,opts.regtyp,ff);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [M,h0] = affreg(buf,x,b0,MG,MF,M,regtyp,ff)
% Do the work

x1 = x{1};
x2 = x{2};
x3 = x{3};
[mu,isig] = priors(regtyp);
isig      =  isig*ff;
Alpha0    =  isig;
Beta0     = -isig*mu;

sol  = M2P(M);
sol1 = sol;
ll   = -Inf;
nsmp = sum(cat(1,buf.nm));
lam  = nsmp*1e-5;
I    = diag([1 1 1  [1 1 1]*pi/180 [1 1 1]*0.01  [1 1 1]*0.01]);
pr   = struct('b',[],'db1',[],'db2',[],'db3',[]);
spm_chi2_plot('Init','Registering','Log-likelihood','Iteration');

for iter=1:200
    penalty = (sol1-mu)'*isig*(sol1-mu);
    T       = MF\P2M(sol1)*MG;
    R       = derivs(MF,sol1,MG);
    y1a     = T(1,1)*x1 + T(1,2)*x2 + T(1,4);
    y2a     = T(2,1)*x1 + T(2,2)*x2 + T(2,4);
    y3a     = T(3,1)*x1 + T(3,2)*x2 + T(3,4);
    h0      = zeros(256,length(b0))+eps;
    for i=1:length(x3),
        if ~buf(i).nm, continue; end;
        y1    = y1a(buf(i).msk) + T(1,3)*x3(i);
        y2    = y2a(buf(i).msk) + T(2,3)*x3(i);
        y3    = y3a(buf(i).msk) + T(3,3)*x3(i);
        for k=1:length(b0),
            pr(k).b = spm_sample_priors(b0{k},y1,y2,y3,k==length(b0));
            h0(:,k) = h0(:,k) + spm_hist(buf(i).g,pr(k).b);
        end;
    end;
    h1  = (h0+eps);
    sh  = sum(h1,1);
    ssh = sum(sh);
    krn = spm_smoothkern(2,(-256:256)',0);
    h1  = conv2(h1,krn,'same');
    h1  = h1/ssh;
    sh  = sh/ssh;
    ll1 = sum(sum(h0.*log2(h1./(sum(h1,2)*sh))))/ssh;

    if ll1-penalty/ssh >ll, cr = ':)'; else cr = ':('; end;
    % fprintf('%3d) %-12.6g - %-12.6g = %-12.6g %s\n',...
    %    iter,ll1,penalty/ssh,ll1-penalty/ssh,cr);

    ll1 = ll1 - penalty/ssh;

    spm_chi2_plot('Set',max(ll1,ll));

    if ll1>ll
        sol   = sol1;
        ll    = ll1;
        Alpha = zeros(12);
        Beta  = zeros(12,1);
        h2    = log2(h1./(sum(h1,2)*sh));
        for i=1:length(x3),
            nz    = buf(i).nm;
            if ~nz, continue; end;
            msk   = buf(i).msk;
            gi    = double(buf(i).g)+1;
            y1    = y1a(msk) + T(1,3)*x3(i);
            y2    = y2a(msk) + T(2,3)*x3(i);
            y3    = y3a(msk) + T(3,3)*x3(i);

            %mi   = zeros(nz,1)+eps;
            dmi1  = zeros(nz,1);
            dmi2  = zeros(nz,1);
            dmi3  = zeros(nz,1);
            for k=1:length(b0),
                [pr(k).b, pr(k).db1, pr(k).db2, pr(k).db3] = spm_sample_priors(b0{k},y1,y2,y3,k==length(b0));
                tmp  = 1    - h2(gi,k);
                %mi  = mi   + tmp.*pr(k).b;
                dmi1 = dmi1 + tmp.*pr(k).db1;
                dmi2 = dmi2 + tmp.*pr(k).db2;
                dmi3 = dmi3 + tmp.*pr(k).db3;
            end;
            x1m = x1(msk);
            x2m = x2(msk);
            x3m = x3(i);
            A = [dmi1.*x1m dmi2.*x1m dmi3.*x1m...
                 dmi1.*x2m dmi2.*x2m dmi3.*x2m...
                 dmi1 *x3m dmi2 *x3m dmi3 *x3m...
                 dmi1      dmi2      dmi3];
            Alpha = Alpha + spm_atranspa(A);
            Beta  = Beta  + sum(A,1)';
        end;
        drawnow;
        Alpha = R'*Alpha*R;
        Beta  = R'*Beta;
        %lam  = lam/2;
    else
        %break;
        % The Levenberg-Marquardt strategy may not really be
        % appropriate.  The cost function decreasing by a
        % tiny amount appears to indicate that interpolation
        % errors close to the solution are making things
        % a tiny bit unstable.  Interpolation could possibly
        % be improved by using higher degree linear interpolation
        % on the log-transformed prior probability images, and
        % normalising so that the interpolated values sum to one.
        lam = lam*10;
        if lam>1e5*nsmp, break; end;
    end;

    % LM update
    sol1  = (Alpha+Alpha0+lam*I)\((Alpha+lam*I)*sol - Beta - Beta0);
end;

spm_chi2_plot('Clear');
M = P2M(sol);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function P = M2P(M)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs
J  = M(1:3,1:3);
V  = sqrtm(J*J');
R  = V\J;

lV = logm(V);
lR = -logm(R);
if sum(sum(imag(lR).^2))>1e-6, error('Rotations by pi are still a problem.'); end;
P       = zeros(12,1);
P(1:3)  = M(1:3,4);
P(4:6)  = lR([2 3 6]);
P(7:12) = lV([1 2 3 5 6 9]);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function M = P2M(P)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs

% Translations
D      = P(1:3);
D      = D(:);

% Rotation part
ind    = [2 3 6];
T      = zeros(3);
T(ind) = -P(4:6);
R      = expm(T-T');

% Symmetric part (zooms and shears)
ind    = [1 2 3 5 6 9];
T      = zeros(3);
T(ind) = P(7:12);
V      = expm(T+T'-diag(diag(T)));

M      = [V*R D ; 0 0 0 1];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function R = derivs(MF,P,MG)
% Numerically compute derivatives of Affine transformation matrix w.r.t.
% changes in the parameters.
R = zeros(12,12);
M0 = MF\P2M(P)*MG;
M0 = M0(1:3,:);
for i=1:12
    dp     = 0.000000001;
    P1     = P;
    P1(i)  = P1(i) + dp;
    M1     = MF\P2M(P1)*MG;
    M1     = M1(1:3,:);
    R(:,i) = (M1(:)-M0(:))/dp;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [mu,isig] = priors(typ)
% The parameters for this distribution were derived empirically from 227
% scans, that were matched to the ICBM space.

mu   = zeros(6,1);
isig = zeros(6);
switch deblank(lower(typ))

case 'mni', % For registering with MNI templates...
        mu   = [0.0667 0.0039 0.0008 0.0333 0.0071 0.1071]';
        isig = 1e4 * [
            0.0902   -0.0345   -0.0106   -0.0025   -0.0005   -0.0163
           -0.0345    0.7901    0.3883    0.0041   -0.0103   -0.0116
           -0.0106    0.3883    2.2599    0.0113    0.0396   -0.0060
           -0.0025    0.0041    0.0113    0.0925    0.0471   -0.0440
           -0.0005   -0.0103    0.0396    0.0471    0.2964   -0.0062
           -0.0163   -0.0116   -0.0060   -0.0440   -0.0062    0.1144];

case 'imni', % For registering with MNI templates...
        mu   = -[0.0667 0.0039 0.0008 0.0333 0.0071 0.1071]';
        isig = 1e4 * [
            0.0902   -0.0345   -0.0106   -0.0025   -0.0005   -0.0163
           -0.0345    0.7901    0.3883    0.0041   -0.0103   -0.0116
           -0.0106    0.3883    2.2599    0.0113    0.0396   -0.0060
           -0.0025    0.0041    0.0113    0.0925    0.0471   -0.0440
           -0.0005   -0.0103    0.0396    0.0471    0.2964   -0.0062
           -0.0163   -0.0116   -0.0060   -0.0440   -0.0062    0.1144];

case 'rigid', % Constrained to be almost rigid...
        mu   = zeros(6,1);
        isig = eye(6)*1e8;

case 'subj', % For inter-subject registration...
        mu   = zeros(6,1);
        isig = 1e3 * [
            0.8876    0.0784    0.0784   -0.1749    0.0784   -0.1749
            0.0784    5.3894    0.2655    0.0784    0.2655    0.0784
            0.0784    0.2655    5.3894    0.0784    0.2655    0.0784
           -0.1749    0.0784    0.0784    0.8876    0.0784   -0.1749
            0.0784    0.2655    0.2655    0.0784    5.3894    0.0784
           -0.1749    0.0784    0.0784   -0.1749    0.0784    0.8876];

case 'none', % No regularisation...
        mu   = zeros(6,1);
        isig = zeros(6);

otherwise
        error(['"' typ '" not recognised as type of regularisation.']);
end;
mu   = [zeros(6,1) ; mu];
isig = [zeros(6,12) ; zeros(6,6) isig];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [h0,d1] = reg_unused(M)
% Try to analytically compute the first and second derivatives of a
% penalty function w.r.t. changes in parameters.  It works for first
% derivatives, but I couldn't make it work for the second derivs - so
% I gave up and tried a new strategy.

T   = M(1:3,1:3);
[U,S,V] = svd(T);
s   = diag(S);
h0  = sum(log(s).^2);
d1s = 2*log(s)./s;
%d2s = 2./s.^2-2*log(s)./s.^2;
d1  = zeros(12,1);

for j=1:3
    for i1=1:9
        T1     = zeros(3,3);
        T1(i1) = 1;
        t1     = U(:,j)'*T1*V(:,j);
        d1(i1) = d1(i1) + d1s(j)*t1;
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function M = P2M_unused(P)
% SVD parameterisation of affine transform, based on matrix-logs.

% Translations
D      = P(1:3);
D      = D(:);

% Rotation U
ind    = [2 3 6];
T      = zeros(3);
T(ind) = P(4:6);
U      = expm(T-T');

% Diagonal zooming matrix
S      = expm(diag(P(7:9)));

% Rotation V'
T(ind) = P(10:12);
V      = expm(T'-T);

M      = [U*S*V' D ; 0 0 0 1];
return;
%_______________________________________________________________________
%_______________________________________________________________________

