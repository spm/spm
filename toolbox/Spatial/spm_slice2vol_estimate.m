function [Q,mu,Mmu,slice_o] = spm_slice2vol_estimate(Nii,slice_code,mm_sd, rad_sd)
% Slice-to-volume alignment estimation
% FORMAT [Q,mu] = spm_slice2vol_estimate(Nii)
%
% Nii - NIfTI data structure encoding volumes to align
%       Must all have the same dimensions
% Q   - A 3D array of slicewise motion parameters
% mu  - Population average
%
% This function has not been thoroughly evaluated yet, but it should serve
% as a useful starting point for coding up several different applications.
% It worked reasonably well for an fMRI time series where the subject moved
% much more than is typical.
%
% Some possible extensions include:
% * Use a higher degree B-spline for pushing and pulling operations, and
%   properly consider slice profiles.
% * Consider a TV regulariser in the template update. See:
%       Brudfors M. Generative Models for Preprocessing of Hospital Brain
%       Scans (Doctoral dissertation, UCL (University College London)).
% * Use a more robust objective function than L2, which may better handle
%   outliers. Alternatively, use a voxel-specific variance (or attempt to
%   model more of the covariance).
% * Consider combining with an unwarping approach to handle EPI
%   distortions. See:
%       Andersson JL, Skare S, Ashburner J. How to correct susceptibility
%       distortions in spin-echo echo-planar images: application to
%       diffusion tensor imaging. Neuroimage. 2003 Oct 1;20(2):870-88.
% * Make it more Bayesian to better handle parameter uncertainty.
% * More nerdy folk may want to improve on the log-Euclidean
%   regularisation.
% * etc
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


% Regularisation settings
if nargin<2, slice_code = 1;            end % Slice ordering
if nargin<3, mm_sd      = 0.001;        end % Translation standard deviation
if nargin<4, rad_sd     = 0.001*pi/180; end % Rotation standard deviation


dat = init_data(Nii);
d   = dat.d;
B   = bases;
Q   = zeros(size(B,3),d(3),d(4));
pad = 3;
dmu = [d(1:3)+2*pad];
mu  = zeros(dmu,'single')+eps;
Mmu = dat.mat*[eye(3) [-pad -pad -pad]'; 0 0 0 1];

% Slice order vector (slice_o)
switch slice_code
case 0 % Unknown - just assume increasing
     slice_o = 1:d(3);
case 1 % sequential increasing
    slice_o = 1:d(3);
case 2 % sequential decreasing
    slice_o = d(3):-1:1;
case 3 % alternating increasing
    slice_o = [1:2:d(3), 2:2:d(3)];
case 4 % alternating decreasing
    if mod(d(3),2) % Odd number of slices
        slice_o = [d(3):-2:1, (d(3)-1):-2:1];
    else
        slice_o = [(d(3)-1):-2:1, d(3):-2:1];
    end
end

% Iterative groupwise slice-to-volume alignment
for iter=1:8
    fprintf('average [%2d]: ', iter); drawnow
    mu  = update_mean(dat,Q,mu,Mmu);
    for subit=1:2
        fprintf(' motion [%2d]: ', iter); drawnow
        Q   = update_motion(dat,Q,mu,Mmu,slice_o,mm_sd,rad_sd);
    end
end

% Undo the padding
mu  = mu((1+pad):(end-pad),(1+pad):(end-pad),(1+pad):(end-pad));
Mmu = Mmu/[eye(3) [-pad -pad -pad]'; 0 0 0 1];


%==========================================================================
function [Q,slice_o] = update_motion(dat,Q,mu,Mmu, slice_o, mm_sd, rad_sd)
%==========================================================================
d   = dat.d;
Mat = dat.mat;
B   = bases;
r1  = 1:d(1); % Subregion to match
r2  = 1:d(2); % Subregion to match

[i1,i2] = ndgrid(single(r1),single(r2));
d(1)    = size(i1,1);
d(2)    = size(i1,2);
K       = size(Q,1);

% Regularisation for each parameter
Reg     = ([mm_sd*[1 1 1] rad_sd*[1 1 1]].^(-2));

% Construct gradients (gT) and Hessians (HT) of the
% likelihood terms for each slice
HT      = zeros(K,K,d(3),d(4));
gT      = zeros(K,d(3),d(4));
phi     = zeros([d(1:3),3],'single');
s2      = 0;
s0      = 0;
for n=1:d(4)
    if ~rem(n,10), fprintf('.'); end
    f            = extract(dat,r1,r2,':',n);
    for i3=1:d(3)
        % Affine mapping from indices in data to indices in template
        Ms            = Mmu\spm_dexpm(Q(:,i3,n),B)*Mat;

        % Spatial transform (indices in template mu)
        phi(:,:,i3,1) = Ms(1,1)*i1 + Ms(1,2)*i2 + (Ms(1,3)*i3+Ms(1,4));
        phi(:,:,i3,2) = Ms(2,1)*i1 + Ms(2,2)*i2 + (Ms(2,3)*i3+Ms(2,4));
        phi(:,:,i3,3) = Ms(3,1)*i1 + Ms(3,2)*i2 + (Ms(3,3)*i3+Ms(3,4));
    end

    % Moved template and its spatial gradients
    [mun,gx,gy,gz] = spm_diffeo('bsplins',mu,phi,[1 1 1 0 0 0]);
    a       = mun - f;

    % Missing values in the images are assumed to be 0
    msk     = (f==0) | ~isfinite(mun) | (mun==0);
    a(msk)  = 0;
    gx(msk) = 0;
    gy(msk) = 0;
    gz(msk) = 0;

    % Increment sufficient statistics for variance estimate
    s2      = s2 + double(sum(a(:).^2));
    s0      = s0 + double(sum(~msk(:)));

    for i3=1:d(3) % Loop over slices

        % Vectorise data for slice
        ai  = vec(a(:,:,i3));
        gxi = vec(gx(:,:,i3));
        gyi = vec(gy(:,:,i3));
        gzi = vec(gz(:,:,i3));

        % Derivatives of mu w.r.t. changes in affine
        % transform elements
        G   = [i1(:).*gxi i1(:).*gyi i1(:).*gzi...
               i2(:).*gxi i2(:).*gyi i2(:).*gzi...
               i3(:).*gxi i3(:).*gyi i3(:).*gzi...
                      gxi        gyi        gzi];

        % Derivative of affine transform elements
        % w.r.t. parameters (for use with chain rule)
        q       = Q(:,i3,n);
        [~,dM3] = spm_dexpm(q,B);
        dM      = zeros(12,size(B,3));
        for m=1:size(B,3)
            tmp     = Mmu\dM3(:,:,m)*Mat;
            dM(:,m) = reshape(tmp(1:3,:),12,1);
        end

        G   = G*dM;       % Chain rule
        g   = G'*ai;      % dW/dq
        H   = G'*G;       % d^2E/dq^2 (Gauss-Newton approx)
        g   = H*q - g;    % note: q <- H\g, not: q <- q - H\g
        gT(:,i3,n)   = g;
        HT(:,:,i3,n) = H;
    end
end

% Update motion parameters as a MAP solution regularised by
% assuming position is correlated with time.  Note that this
% uses a log-Euclidean approximation, rather than properly
% considers the manifold on which the transformation matrices
% fall.
v       = s2/s0; % Maximum likelihood estimate of variance

% Construct gradients (g) and Hessian (H) for the entire set
% of parameters
Ks      = numel(gT);          % Total number of parameters
H       = sparse([],[],[],Ks,Ks,numel(HT)+Ks*4);
g       = zeros(Ks,1);
[ii,jj] = ndgrid(1:K,1:K);
index   = [repmat(slice_o(:),[d(4) 1]) kron((1:d(4))',ones(d(3),1))];
for i=1:size(index,1)

    % Insert the likelihood parts into H
    i3  = index(i,1);
    n   = index(i,2);
    o   = K*(i3-1+d(3)*(n-1)); % Current offset
    h   = HT(:,:,i3,n);
    H   = H + sparse(ii(:)+o,jj(:)+o,h(:)/v,Ks,Ks);

    % Likelihood parts of g
    g((1:K)+o) = gT(:,i3,n)/v;

    if i>1
        % Regularisation part of H using parameters for
        % current and previous (in time) slice
        i3d = index(i-1,1);
        nd  = index(i-1,2);
        od  = K*(i3d-1+d(3)*(nd-1)); % Previous offset
        H   = H + sparse((1:K)+o ,(1:K)+o , Reg,Ks,Ks);
        H   = H + sparse((1:K)+od,(1:K)+o ,-Reg,Ks,Ks);
    end
    if i<size(index,1)
        % Regularisation part of H using parameters for
        % current and next (in time) slice
        i3d = index(i+1,1);
        nd  = index(i+1,2);
        od  = K*(i3d-1+d(3)*(nd-1)); % Next offset
        H   = H + sparse((1:K)+o ,(1:K)+o , Reg,Ks,Ks);
        H   = H + sparse((1:K)+od,(1:K)+o ,-Reg,Ks,Ks);
    end
end

% Gauss-Newton update (note q <- H\g, rather than q <- q - H\g)
Q = reshape(full(H\g),[K d(3) d(4)]);

fprintf('  %g\n', s2/s0); % Show variance estimate

% % Visualise
% t  = (1:d(3)*d(4))'/d(3);
% subplot(4,1,3);
% plot(t,reshape(Q(1:3,slice_o,:),3,size(Q,2)*size(Q,3))','-');
% subplot(4,1,4);
% plot(t,reshape(Q(4:6,slice_o,:),3,size(Q,2)*size(Q,3))','-');
% drawnow


%==========================================================================
function x = vec(X)
%==========================================================================
% Vectorise and set NaNs to 0
x = X(:);
x(~isfinite(x)) = 0;


%==========================================================================
function mu = update_mean(dat,Q,mu,Mmu)
%==========================================================================
d   = dat.d;
Mat = dat.mat;
B   = bases;
[i1,i2] = ndgrid(single(1:d(1)),single(1:d(2)));
phi = zeros([d(1:3),3],'single');
G   = zeros(size(mu),'single');     % Gradient
H   = zeros(size(mu),'single')+eps; % Hessian (diagonal)
s2  = 0;
s0  = 0;
for n=1:d(4)
    if ~rem(n,10), fprintf('.'); end
    f       = extract(dat,':',':',':',n);

    for i3=1:d(3)
        % Affine mapping from indices in data to indices in template
        Ms            = Mmu\spm_dexpm(Q(:,i3,n),B)*Mat;

        % Spatial transform (indices in template mu)
        phi(:,:,i3,1) = Ms(1,1)*i1 + Ms(1,2)*i2 + (Ms(1,3)*i3+Ms(1,4));
        phi(:,:,i3,2) = Ms(2,1)*i1 + Ms(2,2)*i2 + (Ms(2,3)*i3+Ms(2,4));
        phi(:,:,i3,3) = Ms(3,1)*i1 + Ms(3,2)*i2 + (Ms(3,3)*i3+Ms(3,4));
    end

    % \Phi \mu (or \mu(\phi) )
    mun     = spm_diffeo('pull',mu,phi);
    a       = f-mun;             % Residuals
    msk     = (f==0) | ~isfinite(mun) | (mun==0); % Exclude missing
    a(msk)  = 0;
    s2      = s2 + sum(a(:).^2); % Variance sufficient stats
    s0      = s0 + sum(~msk(:)); % Variance sufficient stats
    a(msk)  = NaN;

    % Note: Push, rather than pull (adjoint of resampling)
    % g <- \Phi^T a, h <- \Phi^T 1
    [gn,hn] = spm_diffeo('push',a,phi,size(mu));
    G       = G + gn;
    H       = H + hn;
end

% Gauss-Newton update
mu = mu + H.\G;

fprintf('  %g\n', s2/s0); % Show variance estimate

% subplot(2,2,1);
% imagesc(mu(:,:,ceil(size(mu,3)/2))');
% axis image xy off;
% subplot(2,2,2);
% imagesc(squeeze(mu(:,ceil(size(mu,2)/2),:))');
% axis xy off;
% drawnow


%==========================================================================
function B = bases
%==========================================================================
% SE(3) - Special Euclidean basis functions
B        = zeros(4,4,6);
B(1,4,1) =  1;
B(2,4,2) =  1;
B(3,4,3) =  1;
B(1,2,4) =  1;
B(2,1,4) = -1;
B(1,3,5) =  1;
B(3,1,5) = -1;
B(2,3,6) =  1;
B(3,2,6) = -1;


%==========================================================================
function dat = init_data(Nii)
%==========================================================================
% Create a data structure to make life easier
% with the extract function (below)
ind = zeros(0,2);
if isa(Nii,'char'), Nii = nifti(Nii); end
for i=1:numel(Nii)
    d   = [Nii(i).dat.dim 1 1];
    d   = d(1:4);
    if i==1
        d1 = d;
    else
        if ~all(d1(1:3) == d(1:3))
            error('Incompatible dimensions');
        end
    end
    ind = [ind; ones(d(4),1)*i (1:d(4))'];
end
dat = struct('Nii',Nii,'ind',ind,'d',[d(1:3) size(ind,1)],'mat',Nii(1).mat);


%==========================================================================
function f = extract(dat,r1,r2,r3,n)
%==========================================================================
% Extract data, which may be spread over a number
% of 3D or 4D NIfTI files
f = single(dat.Nii(dat.ind(n,1)).dat(r1,r2,r3,dat.ind(n,2)));
