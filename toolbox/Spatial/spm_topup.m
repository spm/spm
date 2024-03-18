function VDM = spm_topup(data, acqorder, FWHM, reg, rinterp, jac, pref, outdir)
% Correct susceptibility distortions using topup
% FORMAT VDM = spm_topup(vol1, vol2, FWHM, reg, save)
% data       - path to first image (blip-up/PA) followed for the second image (blip-down/AP)
% acqorder   - indicates in which order were acquired the images with
%              different polarities
%               - 0 Blip up first (PA - AP)
%               - 1 Blip down first (AP - PA)
% fwhm       - Gaussian kernel spatial scales (default: [8 4 2 1 0])
% reg        - regularisation settings (default: [0 10 100])
%            See spm_field for details:
%               - [1] Penalty on absolute values.
%               - [2] Penalty on the `membrane energy' of the deformation. 
%                 This penalises the sum of squares of the gradients of the
%                 values.
%               - [3] Penalty on the `bending energy'. This penalises
%                  the sum of squares of the 2nd derivatives of the parameters.
% rinterp    - Order of B-spline by which the images are sampled. A higher 
%              degree provides the better interpolation but it is slower.
% jac        - Option to include jacobian scaling in the process to take 
%              into account the changes of intensities due to stretching 
%              and compression.
% pref       - string to be prepended to the VDM files.
% outdir     - output directory.
%
% VDM        - voxel displacement map.
%
% This is a re-implementation of the topup approach in FSL. Perhaps it
% should be renamed to "putup" or something similar to avoid confusion.
%
% Reference:
%
% J.L.R. Andersson, S. Skare, J. Ashburner. How to correct susceptibility
% distortions in spin-echo echo-planar images: application to diffusion
% tensor imaging. Neuroimage, 20(2):870-888, 2003.
% https://doi.org/10.1016/s1053-8119(03)00336-7
%__________________________________________________________________________

% John Ashburner, Nicole Labra Avila
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

%-Optional input parameters
%--------------------------------------------------------------------------
if nargin < 2
    acqorder = 0;         % order of Blip up/ down acquisition
end
if nargin < 3
    FWHM = [8 4 2 1 0];   % spatial scales
end
if nargin < 4
    reg = [0 10 100];     % regularisation
end
if nargin < 5
    rinterp = 2;          % Degree of B-spline
end
if nargin < 6
    jac = 1;               % Refine topip
end
if nargin < 7
    pref = 'vdm5_';       % Prefix for files
end
if nargin < 8
    outdir = '';          % output directory
end

%-Check the amount of Blip-up/Blip-down pairs entered. If there is more
% than one pair, the average image of every group (Blip-ups or Blip-downs)
% are computed with spm_realign.m and Topup is executed using the average
% images.
%--------------------------------------------------------------------------
% Default parameters for spm_realign
flags1.quality = 0.95;
flags1.fwhm    = 2;
flags1.sep     = 2;
flags1.rtm     = 1;
flags1.PW      = '';
flags1.interp  = 2;
flags1.wrap    = [0 0 0];

% Default parameters for spm_reslice
flags2.which   = [2 1];
flags2.interp  = 4;
flags2.wrap    = [0 0 0];
flags2.mask    = 1;
flags2.prefix  = 'r';

vol1 = cell(size(data));
vol2 = cell(size(data));

for i = 1:numel(data)
   vol1{i} = char(data(i).volbup{:});
   vol2{i} = char(data(i).volbdown{:});
end
dvol1 = size(vol1{1});
dvol2 = size(vol2{1});

% Check if the input volumes for blip-up and blip-down have the same amount
% of elements. Send an error and stop the program if not. If yes, check how
% many elements they have. If the amount of images per phase encoding
% direction is > 1, realign the images and obtain a mean image per group
% (blip-up and blip-downs)

%auxvol = cellstr(size(vol1));
if dvol1(1) == dvol2(1)
   if dvol1 > 1     
      if acqorder == 0
         vol1 = cellstr(flipud((vol1{:})));
      elseif acqorder == 1
           vol2 = cellstr(flipud((vol2{:})));
      end 
      spm_realign(vol1,flags1);
      spm_reslice(vol1,flags1)
      spm_realign(vol2,flags2);
      spm_reslice(vol2,flags2)
      basename = spm_file(vol1{1}(1,:),'basename');
      basename = spm_file(basename,'prefix',"mean",'ext','.nii');
      vol1 = [spm_file(vol1{1}(1,:),'fpath') '/' basename];
      basename = spm_file(vol2{1}(1,:),'basename');
      basename = spm_file(basename,'prefix',"mean",'ext','.nii');
      vol2 = [spm_file(vol2{1}(1,:),'fpath') '/' basename];
   else
        vol1 = char(vol1);
        vol2 = char(vol2);
    end
 
else
   error("Number of blip-up volumes is not consistent with number of " + ...
       "blip-down volumes");
end 

%-Load images and estimate noise level
%--------------------------------------------------------------------------
ord  = [0 rinterp 0  0 0 0];
P    = strvcat(fullfile(vol1), fullfile(vol2));
Nii  = nifti(P);
vx   = sqrt(sum(Nii(1).mat(1:3,1:3).^2)); % Voxel sizes
sd   = spm_noise_estimate(P);             % Rice mixture model to estimate image noise
sig2 = sum(sd.^2);                        % Variance of difference is sum of variances

spm_field('bound',1);                     % Set boundary conditions for spm_field

f1_0  = single(Nii(1).dat(:,:,:));
f2_0  = single(Nii(2).dat(:,:,:));
%f2_0 = f2_0*(mean(f1_0(:))/mean(f2_0(:)));

d   = size(Nii(1).dat);
u   = zeros(d,'single'); % Starting estimates of displacements

% Set the windows figure
FG = spm_figure('GetWin','Graphics');

%-Topup
%==========================================================================
for fwhm = FWHM % Loop over spatial scales

    % Smooth the images according to the current spatial scale
    %----------------------------------------------------------------------
    fprintf('%4g:', fwhm);
    f1  = spm_diffeo('bsplinc', f1_0,ord);
    f2  = spm_diffeo('bsplinc', f2_0,ord);
    spm_smooth(f1,f1,[1 1 1]*fwhm); % Note the side effects
    spm_smooth(f2,f2,[1 1 1]*fwhm); % Note the side effects

    if jac==1
        [u,wf1,wf2] = topup_jacobians(u, f1,f2, sig2, vx, reg, ord, 2e-3, 6, FG, fwhm);
    else
        [u,wf1,wf2] = topup_basic(u, f1,f2, sig2, vx, reg, ord, 2e-3, 8, FG, fwhm);
    end
end


%-Save distortion-corrected blip up/down and down/up images) and VDM
%==========================================================================

% wf1 (blip up/down image)
%--------------------------------------------------------------------------
basename = spm_file(vol1,'basename');
pr       = strcat('w', pref);
oname    = spm_file(basename,'prefix',pr,'ext','.nii');
oname    = fullfile(outdir,oname);
Nio      = nifti;
Nio.dat  = file_array(oname,d,'float32');
Nio.mat  = Nii(1).mat;
Nio.mat0 = Nii(1).mat;
create(Nio);
Nio.dat(:,:,:) = wf1;

% wf2 (blip down/up image)
%--------------------------------------------------------------------------
basename = spm_file(vol2,'basename');
pr       = strcat('w', pref);
oname    = spm_file(basename,'prefix',pr,'ext','.nii');
oname    = fullfile(outdir,oname);
Nio      = nifti;
Nio.dat  = file_array(oname,d,'float32');
Nio.mat  = Nii(2).mat;
Nio.mat0 = Nii(2).mat0;
create(Nio);
Nio.dat(:,:,:) = wf2;

% Write VDM file
%--------------------------------------------------------------------------
basename       = spm_file(vol1,'basename');
pr             = strcat(pref,'_pos_');
oname          = spm_file(basename,'prefix',pr,'ext','.nii');
oname          = fullfile(outdir,oname);
Nio            = nifti;
Nio.dat        = file_array(oname,size(u),'float32');
Nio.mat        = Nii(1).mat;
create(Nio);
Nio.dat(:,:,:) = u;

VDM = Nio;

end
% ===========================================================================

% ===========================================================================
function [u,wf1,wf2] = topup_basic(u, f1,f2, sig2, vx, reg, ord, tol, nit, FG, fwhm)
% Basic topup, without accounting for jacobian modulation for stretching and compression.
if isempty(u)
    u   = zeros(size(f1));
end
if nargin<10|| isempty(FG),   FG   = 0;     end
if nargin<9 || isempty(nit),  nit  = 10;    end
if nargin<8 || isempty(tol),  tol  = 2e-3; end
if nargin<7 || isempty(ord),  ord  = [1 1 1  0 0 0]; end
if nargin<6 || isempty(reg),  reg  = [0 10 100]; end
if nargin<5 || isempty(vx),   vx   = [1 1 1]; end

d   = size(f1);
id  = identity2(d);
phi = id; 

n_acceptable = 0;
E = 0;

% Iterate until convergence (10 iterations max)
%----------------------------------------------------------------------
for it = 1:nit

    % Sample the blip up/down image and its gradients in y direction
    phi(:,:,:,2)  = id(:,:,:,2) + u;
    [wf1,~,gf1,~] = spm_diffeo('bsplins',f1,phi,ord);

    % Sample the blip down/up image and its gradients in y direction
    phi(:,:,:,2)  = id(:,:,:,2) - u;
    [wf2,~,gf2,~] = spm_diffeo('bsplins',f2,phi,ord);

    % Regularisation term is \tfrac{1}{2} u^T L u. Compute L u.
    Lu  = spm_field('vel2mom', u, [vx reg]);

    % Compute cost function. Note the slightly ad hoc treatment of
    % missing data (where one of the phis points outside the FOV).
    % Note that there isn't a clear underlying generative model
    % of the data (so "log-likelihood" is in scare quotes).
    res = 0.5*((wf1 - wf2).^2/sig2 + Lu.*u);
    msk = isfinite(res);
    Eprev = E;
    E   = sum(res(msk(:)))/sum(msk(:));
    fprintf(' %7.4f', E)

    % Gradient of "log-likelihood" term
    g   = (gf1 + gf2).*(wf1 - wf2)/sig2;

    % Diagonal of Hessian of "log-likelihood" term
    h   = (gf1 + gf2).^2/sig2;

    % Mask out missing data
    msk = ~isfinite(g);
    g(msk) = 0;
    h(msk) = 0;

    g   = g + Lu; % Gradient of "log-likelihood" and log-prior terms

    % Gauss-Newton update
    % u \gets u - (diag(h) + L)\(g + L u)
    u   = u - spm_field(h, g, [vx reg 3 3]);

    % Check convergence
    if abs(E - Eprev) < tol                 % Stopping criterion
        n_acceptable = n_acceptable + 1;    % Closer to stopping
    else
        n_acceptable = 0;                   % Restart the countdown
    end
    if n_acceptable >= 3, break; end        % Done

    do_display(wf1,wf2,u,[],g, FG, fwhm, it);
end
fprintf('\n');

if nargout>=2
    phi(:,:,:,2) = id(:,:,:,2) + single(u);
    wf1 = spm_diffeo('bsplins',f1,phi,ord);
    phi(:,:,:,2) = id(:,:,:,2) - single(u);
    wf2 = spm_diffeo('bsplins',f2,phi,ord);
end

end
% ===========================================================================

% ===========================================================================
function [u,wf1,wf2] = topup_jacobians(u, f1,f2, sig2, vx, reg, ord, tol, nit, FG, fwhm)
% It include in the process the changes of intensities due to
% stretching and compression.
if isempty(u)
    u   = zeros(size(f1));
end
if nargin<10|| isempty(FG),   FG   = 0;     end
if nargin<9 || isempty(nit),  nit  = 10;    end
if nargin<8 || isempty(tol),  tol  = 2e-3; end
if nargin<7 || isempty(ord),  ord  = [1 1 1  0 0 0]; end
if nargin<6 || isempty(reg),  reg  = [0 10 100]; end
if nargin<5 || isempty(vx),   vx   = [1 1 1]; end

d   = size(f1);
L   = regop(d, vx, reg);
G   = vx(2)*kron(kron(idop(d(3)),diffop(d(2),vx(2))),idop(d(1))); % Difference operator
G   = (G-G')/2; % Gradient operator (c.f. grad in MATLAB)
u   = double(u);
id  = identity2(d);
phi = id;
filter(d,vx,reg);
E   = Inf;

% Iterate until convergence
%----------------------------------------------------------------------
for it = 1:nit

    phi(:,:,:,2)  = id(:,:,:,2) + single(u);
    [wf1,~,gf1,~] = spm_diffeo('bsplins',f1,phi,ord);

    phi(:,:,:,2)  = id(:,:,:,2) - single(u);
    [wf2,~,gf2,~] = spm_diffeo('bsplins',f2,phi,ord);

    % Deal with potential missing data
    msk = find(~isfinite(wf1) | ~isfinite(wf2));
    wf1(msk) = 0;
    wf2(msk) = 0;
    gf1(msk) = 0;
    gf2(msk) = 0;
    nvox     = numel(u) - numel(msk);

    Gu    = reshape(full(G*u(:)),d);
    Lu    = reshape(full(L*u(:)),d); % L*u

    b     = wf1.*(1+Gu) - wf2.*(1-Gu);
    Eprev = E;
    E     = b(:)'*b(:)/(2*sig2) + 0.5*(u(:)'*Lu(:));

    fprintf(' %7.4f', E/nvox)
    if (Eprev-E)/nvox < tol, break; end

    % Construct gradient (g) and Hessian (H) for the Gauss-Newton update. Note that
    % the Hessian is only an approximation, but it is positive definite.
    At = spdiags(double(gf1(:).*(1+Gu(:))+gf2(:).*(1-Gu(:))),0,numel(u),numel(u)) ...
         - double(wf1(:)+wf2(:))'.*G;         % A'
    g   = full(At*double(b(:)))/sig2 + Lu(:); % A'*b + L*u

    % H(i,i) (the diagonal) matches the true Hessian
    % H(i,i+1) & H(i,i-1) do not match
    % H(i,i+2) & H(i,i-2) match the true Hessian
    % H(i,(i+3):end) &  H(i,1:(i-3)) match and are all zero
    H   = (At*At')/sig2;
    h   = reshape(single(full(diag(H))),size(u));
    H   = H + L;          % A'*A + L 

    do_display(wf1,wf2,u,Gu,g,FG,fwhm,it);

    % Gauss-Newton update via a preconditioned conjugate gradient.  It
    % would be nicer if SPM's multigrid could handle this form of Hessian
    % - but sadly it can't.  Conjugate gradient is slow and often does not
    % converge.
    [du,flg,relres,numit] = cgs(H,g,1e-1,1000,@precon,[],[]);
    u(:) = u(:) - full(du(:));
    %fprintf(' (%d,%7.4f,%4d)', flg, relres, numit)
end

fprintf('\n');

if nargout>=2
    jd  = reshape(full(G*u(:)),d);
    phi(:,:,:,2) = id(:,:,:,2) + single(u);
    wf1 = spm_diffeo('bsplins',f1,phi,ord).*(1+jd);
    phi(:,:,:,2) = id(:,:,:,2) - single(u);
    wf2 = spm_diffeo('bsplins',f2,phi,ord).*(1-jd);
end
u = single(u);
end
% ===========================================================================

% ===========================================================================
function do_display(wf1,wf2,u,jd,g,FG,fwhm,it)
% Display topup iterations in Graphics window
if FG==0, return; end
pl = ceil(size(wf1,3)*0.5);

pic1 = wf1(:,:,pl);
pic2 = wf2(:,:,pl);
picu = u(:,:,pl);
if ~isempty(jd)
    jds  = jd(:,:,pl);
    pic1 = pic1.*(1+jds);
    pic2 = pic2.*(1-jds);
end

sgtitle(['Topup! ----> fwhm : ' num2str(fwhm) ' it: ' num2str(it)]);
ax = subplot(3,2,1,'Parent',FG);
imagesc(ax, pic1');
axis(ax,'image','xy','off');

ax = subplot(3,2,2,'Parent',FG);
imagesc(ax, pic2');
axis(ax,'image','xy','off');

ax = subplot(3,2,3,'Parent',FG);
imagesc(ax,picu');
axis(ax,'image','xy','off');
colorbar(ax);
title(ax,'Displacement');

ax = subplot(3,2,4,'Parent',FG);
imagesc(ax, (pic1-pic2)');
axis(ax,'image','xy','off');
colorbar(ax);
title(ax,'Difference');

if ~isempty(g)
    g   = reshape(g,size(u));
    pic = g(:,:,pl);
    ax  = subplot(3,2,5,'Parent',FG);
    imagesc(ax, pic');
    axis(ax,'image','xy','off');
    colorbar(ax);
    title(ax,'1st derivatives');
end

drawnow
end
% ===========================================================================

% %==========================================================================
function L = regop(d, vx, reg, bnd)
% Generate the regulaser used by spm_field
    if nargin<4,      bnd = [1 1 1]; end
    if numel(bnd)==1, bnd = [1 1 1]*bnd; end

    d  = [d(:)'  1 1];
    vx = [vx(:)' 1 1];
    Ix = idop(d(1));
    Iy = idop(d(2));
    Iz = idop(d(3));
    Dx = diffop(d(1),vx(1),bnd(1));
    Dy = diffop(d(2),vx(2),bnd(2));
    Dz = diffop(d(3),vx(3),bnd(3));
    L  = [kron(kron(Iz,Iy),Dx);
          kron(kron(Iz,Dy),Ix);
          kron(kron(Dz,Iy),Ix)];
    L  = L'*L;
    L  = reg(1)*kron(kron(Iz,Iy),Ix) + reg(2)*L + reg(3)*(L'*L);
end

% ===========================================================================
% %==========================================================================
function D = diffop(d, vx, bc)
% Sparse matrix representation of difference operation
    if nargin<3, bc = 0; end
    if nargin<2, vx = 1; end
    if d>=2
        D = (sparse(1:d,1:d,-1,d,d) + sparse(1:d,[2:d 1],1,d,d))/vx;
        switch bc
        case 0
        case 1
            D(end,:) = 0;
        end
    else
        D = sparse([],[],[],d,d);
    end
end

% ===========================================================================
% %==========================================================================
function I = idop(d)
% Sparse identity transform
    I  = sparse(1:d,1:d,1,d,d);
end

% ===========================================================================
% %==========================================================================
function id = identity2(d)
% Identity transform (alternative to concatenating outputs from ndgrid).
    d  = [d(:)' 1 1];
    d  = d(1:3);
    ix = ones(d(1),1,'single');
    iy = ones(d(2),1,'single');
    iz = ones(d(3),1,'single');
    rx = (single(1):single(d(1)))';
    ry = (single(1):single(d(2)))';
    rz = (single(1):single(d(3)))';
    id = reshape([kron(iz,kron(iy,rx))...
                  kron(iz,kron(ry,ix))...
                  kron(rz,kron(iy,ix))],[d 3]);
end

% ===========================================================================
% %==========================================================================
function Fi = filter(d,vx,prm,c)
% Generate a Fourier transform of the Greens function
% for use as a preconditioner in the conjugate gradient
% solver. Note that the preconditioner involves DCTs, rather
% than FFTs, in order to allow more desirable boundary conditions.
    persistent F
    if nargin>=2
      if nargin<4, c = 0.1; end
        d2 = d;
        d2(d2>1) = d2(d2>1)*2;
        l = spm_diffeo('kernel',d2(1:3),[vx prm 0 0])/vx(2);
        l(1,1,1) = l(1,1,1) + c;
        F = real(fftn(double(l)));
        F = 1./F(1:d(1),1:d(2),1:d(3));
    end
    Fi = F;
end

% ===========================================================================
% %==========================================================================
function v = precon(u)
% Preconditioner for the conjugate gradient solver. c.f. the use of
% the "Hilbert gradient" in Beg's LDDMM paper.
    F   = filter;
    fun = spm_dctdst;
    v   = fun.idct2n(fun.dct2n(reshape(u,size(F))).*F);
    v   = v(:);
end

% ===========================================================================
% %==========================================================================

