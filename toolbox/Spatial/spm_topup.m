% %==========================================================================
% %==========================================================================
function VDM = spm_topup(vol1, vol2, FWHM, reg, rinterp, rt, pref, outdir)
% Correct susceptibility distortions using topup
% FORMAT VDM = spm_topup(vol1, vol2, FWHM, reg, save)
% vol1       - path to first image (blip up)
% vol2       - path to second image (blip down)
% fwhm       - Gaussian kernel spatial scales (default: [8 4 2 1 0.1])
% reg        - regularisation settings (default: [0 10 100])
%            See spm_field for details:
%               - [1] Penalty on absolute values.
%               - [2] Penalty on the `membrane energy'. This penalises
%                  the sum of squares of the gradients of the values.
%               - [3] Penalty on the `bending energy'. This penalises
%                  the sum of squares of the 2nd derivatives.
% rinterp    - Degree of B-spline 
% rt         - Option to apply a supplementary refine over topup to include in the 
%              process the changes of intensities due to stretching and compression.
% pref       - string to be prepended to the VDM files.
% outdir     - output directory.
%
% VDM        - voxel displacement map.
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
if nargin < 3
    FWHM = [8 4 2 1 0.1]; % spatial scales
end
if nargin < 4
    reg = [0 10 100];     % regularisation
end
if nargin < 5
    rinterp = [1 1 1];    % Degree of B-spline 
end
if nargin < 6
    rt = 1;               % Refine topip 
end
if nargin < 7
    pref = 'vdm5_';       % Prefix for files
end

if nargin < 8
    outdir = '';          % output directory
end


%-Load images and estimate noise level
%--------------------------------------------------------------------------
P   = strvcat(fullfile(vol1), fullfile(vol2));
Nii = nifti(P);
vx  = sqrt(sum(Nii(1).mat(1:3,1:3).^2)); % Voxel sizes
reg = [vx reg];                          % Regularisation settings (spm_field)
sd  = spm_noise_estimate(P);             % Rice mixture model to estimate image noise
sig2 = sum(sd.^2);                       % Variance of difference is sum of variances

spm_field('bound',1);                    % Set boundary conditions for spm_field

% Set up identity transform (id) and displacements (u)
%--------------------------------------------------------------------------
d   = size(Nii(1).dat);
id  = zeros([d 3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));

u   = zeros(d,'single'); % Starting estimates of displacements
phi = id;                % x and z components are identity.
                         % y component will have displacement added/subtracted

% Set the windows figure                        
FG = spm_figure('GetWin','Graphics');   

%-Topup
%==========================================================================
for fwhm = FWHM % Loop over spatial scales
    fprintf('%3g:', fwhm);

    % Smooth the images according to the current spatial scale
    %----------------------------------------------------------------------
    f1  = single(Nii(1).dat(:,:,:));
    f2  = single(Nii(2).dat(:,:,:));
    spm_smooth(f1,f1,[1 1 1]*fwhm); % Note the side effects
    spm_smooth(f2,f2,[1 1 1]*fwhm); % Note the side effects

    n_acceptable = 0;
    E = 0;

    % Iterate until convergence (10 iterations max)
    %----------------------------------------------------------------------
    for it = 1:10

        % Sample the blip up/down image and its gradients in y direction
        phi(:,:,:,2) = id(:,:,:,2) + u;
        [wf1,~,d1,~] = spm_diffeo('bsplins',f1,phi,[rinterp [0 0 0]]);

        % Sample the blip down/up image and its gradients in y direction
        phi(:,:,:,2) = id(:,:,:,2) - u;
        [wf2,~,d2,~] = spm_diffeo('bsplins',f2,phi,[rinterp [0 0 0]]);

        % Regularisation term is \tfrac{1}{2} u^T L u. Compute L u. 
        gu  = spm_field('vel2mom', u, reg);

        % Compute cost function. Note the slightly ad hoc treatment of
        % missing data (where one of the phis points outside the FOV).
        % Note that there isn't a clear underlying generative model
        % of the data (so "log-likelihood" is in scare quotes).
        res = 0.5*((wf1 - wf2).^2/sig2 + gu.*u);
        msk = isfinite(res);
        Eprev = E;
        E   = sum(res(msk(:)))/sum(msk(:));
        fprintf(' %7.4f', E)

        % Gradient of "log-likelihood" term
        g   = (d1 + d2).*(wf1 - wf2)/sig2;
        
        % Diagonal of Hessian of "log-likelihood" term
        h   = (d1 + d2).^2/sig2;

        % Mask out missing data
        msk = ~isfinite(g);
        g(msk) = 0;
        h(msk) = 0;

        g   = g + gu; % Gradient of "log-likelihood" and log-prior terms

        % Gauss-Newton update
        % u \gets u - (diag(h) + L)\(g + L u)
        u   = u - spm_field(h, g, [reg 3 3]);

        % Check convergence
        if abs(E - Eprev) < 1e-3                % Stopping criterion
            n_acceptable = n_acceptable + 1;    % Closer to stopping
        else
            n_acceptable = 0;                   % Restart the countdown
        end
        if n_acceptable >= 3, break; end        % Done

        % Display topup iterations in Graphics window
        %------------------------------------------------------------------
        if true
            pl = ceil(size(f1,3)*0.2);
      
            ax = subplot(2,2,1,'Parent',FG);
            imagesc(ax,[f1(:,:,pl)' f2(:,:,pl)']);
            axis(ax,'image','xy','off');
            title(ax,'Original');

            ax = subplot(2,2,2,'Parent',FG);
            imagesc(ax,[wf1(:,:,pl)' wf2(:,:,pl)']);
            axis(ax,'image','xy','off');
            title(ax,'Warped');

            ax = subplot(2,2,3,'Parent',FG);
            imagesc(ax,u(:,:,pl)');
            axis(ax,'image','xy','off');
            colorbar(ax);
            title(ax,'VDM');

            ax = subplot(2,2,4,'Parent',FG);
            imagesc(ax,g(:,:,pl)');
            axis(ax,'image','xy','off');
            colorbar(ax);
            title(ax,'dE/du');

            drawnow

        end
    end
    fprintf('\n');
end

% Refine Topup
if rt == 1
   [u,wf1,wf2] = refine_topup(u, f1,f2, sig2, vx, reg(4:end));
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

% Write VDM file
%--------------------------------------------------------------------------
vi             = u*(-1);
basename       = spm_file(vol2,'basename');
pr             = strcat(pref,'_neg_');
oname          = spm_file(basename,'prefix', pr,'ext','.nii');
oname          = fullfile(outdir,oname);
Nio            = nifti;  
Nio.dat        = file_array(oname,size(vi),'float32');
Nio.mat        = Nii(1).mat;  
create(Nio);          
Nio.dat(:,:,:) = vi; 

end

% ===========================================================================
% %==========================================================================output directory
function [u,wf1,wf2] = refine_topup(u, f1,f2, sig2, vx, reg, ord, tol, nit)
% Refine topup. It include in the process the changes of intensities due to 
% stretching and compression.
fprintf('Refine topup')
if isempty(u)
    u   = zeros(size(f1));
end
if nargin<9 || isempty(nit), nit = 6;     end
if nargin<8 || isempty(tol), tol = 0.001; end
if nargin<7 || isempty(ord), ord = [1 1 1  0 0 0]; end
if nargin<6 || isempty(reg), reg = [0 10 100]; end
if nargin<5 || isempty(vx),  vx  = [1 1 1]; end

d   = size(f1);
L   = regop(d, vx, reg);
D   = vx(2)*kron(kron(idop(d(3)),diffop(d(2),vx(2))),idop(d(1)));
D   = (D-D')/2;
u   = double(u);
id  = identity2(d);
phi = id;
filter(d,vx,reg);

if false
    % Testing the operator
    spm_field('bound',1);
    gu1 = full(L*u(:));
    gu2 = spm_field('vel2mom',single(u),[vx reg]);
    disp(norm(gu1(:)-gu2(:))./norm(gu2(:)))
end

E   = Inf;

% Iterate until convergence
%----------------------------------------------------------------------
for it = 1:nit

    phi(:,:,:,2) = id(:,:,:,2) + single(u);
    [wf1,~,d1,~] = spm_diffeo('bsplins',f1,phi,ord);

    phi(:,:,:,2) = id(:,:,:,2) - single(u);
    [wf2,~,d2,~] = spm_diffeo('bsplins',f2,phi,ord);

    % Deal with potential missing data
    msk = find(~isfinite(wf1) | ~isfinite(wf2));
    wf1(msk) = 0;
    wf2(msk) = 0;
    d1(msk)  = 0;
    d2(msk)  = 0;
    nvox     = numel(u) - numel(msk);

    jd    = reshape(full(D*u(:)),d);
    gu    = reshape(full(L*u(:)),d); % L*u
    b     = wf1.*(1-jd) - wf2.*(1+jd);
    Eprev = E;
    E     = b(:)'*b(:)/(2*sig2) + 0.5*(u(:)'*gu(:));
    fprintf(' %7.4f', E/nvox)
    if (Eprev-E)/nvox < tol, break; end

    % Generate a very approximate solution using spm_field, 
    % which can be plugged in as a starting estimate for the
    % conjugate gradient solver.  The multigrid implementation
    % can not yet handle off-diagnoals in the likelihood part
    % of the Hessian matrix.
    g  = b.*(d1.*(1-jd) + d2.*(1+jd))/sig2 + single(gu);
    h  =    (d1.*(1-jd) + d2.*(1+jd)).^2/sig2;
    spm_field('bound',1); % Free boundary (needed for top and bottom slices)
    du = spm_field(h,g,[vx reg 2 2]);

    % Construct gradient (g) and Hessian (H) for the Gauss-Newton update.
    At = spdiags(double(d1(:)+d2(:)),0,numel(u),numel(u)) +...
         double(wf1(:)+wf2(:))'.*D;  % A'
    g   = full(At*double(b(:)))/sig2 + gu(:); % A'*P*b + L*u
    H   = (At*At')/sig2 + L;          % A'*P*A + L

    % Gauss-Newton update via a preconditioned conjugate gradient.  It
    % would be nicer if SPM's multigrid could handle this form of Hessian
    % - but sadly it can't.  Conjugate gradient is slow and often does not
    % converge.
    [du,flg,relres,numit] = cgs(H,g,1e-2,1000,@precon,[],double(du(:)));
    u(:) = u(:) - full(du(:));
    %fprintf(' (%d,%7.4f,%4d)', flg, relres, numit)
end
fprintf('\n');
if nargout>=2
    jd  = reshape(full(D*u(:)),d);
    phi(:,:,:,2) = id(:,:,:,2) + single(u);
    wf1 = spm_diffeo('bsplins',f1,phi,ord).*(1-jd);
    phi(:,:,:,2) = id(:,:,:,2) - single(u);
    wf2 = spm_diffeo('bsplins',f2,phi,ord).*(1+jd);
end
u = single(u);
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
function Fi = filter(d,vx,prm)
% Generate a Fourier transform of the Greens function
% for use as a preconditioner in the conjugate gradient
% solver
    persistent F
    if nargin>=2
        l = spm_diffeo('kernel',d(1:3),[vx prm 0 0])/vx(2);
        l(1,1,1) = l(1,1,1) + 128; % Ad hoc
        F = real(1./fftn(double(l)));
    end
    Fi = F;
end

% ===========================================================================
% %==========================================================================
function v = precon(u)
% Preconditioner for the conjugate gradient solver.  This
% is slightly ad hoc, but it seems to help.
    F = filter;
    v = real(ifftn(fftn(reshape(u,size(F))).*F));
    v = v(:);
end

% ===========================================================================
% %==========================================================================

