function y = spm_invert_def(y,M1,d0,M0,args)
% Invert a deformation field
% FORMAT iY = spm_invert_def(Y,M,idim,iM,args)
% Y    - input deformation field (floating point n1*n2*n3*3)
% M    - voxel-to-world mapping of deformation field
% idim - dimensions of inverse deformation
% iM   - voxel-to-world mapping of inverse deformation
% args - additional arguments (optional)
%        * If the first element is non-zero, then assume the
%          original mapping is to voxel indices, rather than to
%          mm coordinates.
%        * If the second element is non-zero, then assume the
%          output mapping is to voxel indices.
%
% iY   - inverse deformation field
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_invert_def.m 3211 2009-06-19 12:34:26Z john $


d1 = size(y);
d1 = d1(1:3);
d0 = d0(1:3);

if nargin<5 || numel(args)<1 || ~args(1),
    % y maps from individual to mm coordinates in the template.  Change
    % this so that it maps to voxel indices in the template instead.
    y      = affind(y,inv(M1));
end

% Generate a grid of mm coordinates for each voxel within the
% individual's image.
x      = affind(rgrid(d1),M0);

% These mm coordinates are pushed to their new locations according
% to the original mapping (y).  Note that the resulting y is scaled
% at each point by the number of voxels that are mapped (w).
[y,w]  = dartel3('push',x,y,d0);

% Generate another grid of mm indices at each voxel in the template.
x      = affind(rgrid(d0),M1);

% Fit a the closest affine transform through the mapping.
M      = spm_get_closest_affine(x,y,w);

% Multiply the mm indices by this affine transform.  This transform
% is subtracted from the nonlinear warp so that a smooth displacement
% field can be determined.  The transform is then re-added to the
% displacement field.
x      = affind(x,M);
vx     = sqrt(sum(M1(1:3,1:3).^2));
for m=1:3,
    % Essentially, divide (y-x*w) by w but avoid divisions by zero by
    % biasing the result to be spatially smooth.
    y(:,:,:,m) = optimNn(w,y(:,:,:,m)-x(:,:,:,m).*w,[2  vx  0.001 1e-6 0  2 1]) + x(:,:,:,m);
end

if nargin>=5 && numel(args)>=2 && args(2),
    y = affind(y,inv(M0));
end

%=======================================================================

%=======================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================

%=======================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3),
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%=======================================================================

%=======================================================================

