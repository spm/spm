function that = spm_swarp(this,def,M)
% Warp surface.
% FORMAT that = spm_swarp(this,def)
% this - a gifti object
% def  - a deformation (nifti object or filename)
% that - the warped gifti object
%
% FORMAT that = spm_swarp(this,def,M)
% this - a gifti object
% def  - a deformation field (nx*ny*nz*1*3)
% M    - mapping from voxels to world, for deformation field
% that - the warped gifti object
%
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_swarp.m 2859 2009-03-11 16:52:19Z guillaume $

if ~isa(this,'gifti'), this = gifti(this); end

if nargin<2, that = this; return; end

if ischar(def) || isa(def,'nifti'),
    if ischar(def),
        def = nifti(def);
    end
    y   = def(1).dat(:,:,:,:,:);
    M   = def(1).mat;
else
    y = def;
    if nargin<3, M = eye(4); end
end

v   = this.vertices;
iM  = inv(M);
v   = iM(1:3,1:4)*[v'; ones(1,size(v,1))];
xyz = {double(v(1,:)'),double(v(2,:)'),double(v(3,:)')};
v   = [spm_bsplins(y(:,:,:,1,1),xyz{:},[1 1 1 0 0 0]),...
       spm_bsplins(y(:,:,:,1,2),xyz{:},[1 1 1 0 0 0]),...
       spm_bsplins(y(:,:,:,1,3),xyz{:},[1 1 1 0 0 0])];

% Much of the surface data is likely to fall outside the FOV of the deformation
% field.  For this reason, the following code attempts to extrapolate the surfaces
% outside this FOV by replacing NaNs in the vertice coordinates by some smooth mesh.

if isfield(this, 'faces')
    f = this.faces;
    m = size(v,1);

    % Gradients
    b = double([v(:,1); v(:,2); v(:,3)]);
    w = isfinite(b);
    b(~w) = 0;

    % Hessian
    A = spdiags(w,0,m*3,m*3);
    f = double(f);
    W = sparse(f(:,1),f(:,1),1,m,m) + sparse(f(:,2),f(:,2),1,m,m) + sparse(f(:,3),f(:,3),1,m,m)...
        - sparse(f(:,1),f(:,2),1,m,m) - sparse(f(:,2),f(:,3),1,m,m) - sparse(f(:,3),f(:,1),1,m,m);
    W = (W'*W + 0.25*W)*1e-2;
    Z = sparse([],[],[],m,m);
    A = A + [W Z Z; Z W Z; Z Z W];

    % Solution to simple linear model.
    v = full(reshape(A\b,m,3));
end

% Generate the gifti structure for the warped data.
that = subsasgn(this,substruct('.','vertices'),v);

