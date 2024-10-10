function varargout = spm_shp_warps(varargin)
%__________________________________________________________________________
% Collection of tools for manipulating non-linear transformations (warps).
%
% FORMAT out = warps(('warp'), in, y, (vs_in), (itrp), (bnd))
% FORMAT y   = warps('compose', y_1, (vs_1), ..., y_n, (vs_n), (itrp))
% FORMAT y   = warps('identity', lat_dim, (lat_vs))
% FORMAT y   = warps('translation', T, lat_dim, (lat_vs))
% FORMAT y   = warps('linear', L, lat_dim, (lat_vs))
% FORMAT y   = warps('affine', A, lat_dim, (lat_vs))
% FORMAT y   = warps('mm2vox', y, vs)
% FORMAT y   = warps('transform', A, y)
%
% FORMAT help warps>function
% Returns the help file of the selected function.
%__________________________________________________________________________

% Yael Balbastre
% Copyright (C) 2017-2024 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help warps
    error('Not enough argument. Type ''help warps'' for help.');
end
if ~ischar(varargin{1})
    [varargout{1:nargout}] = warp(varargin{:});
else
    id = varargin{1};
    varargin = varargin(2:end);
end
switch lower(id)
    case 'warp'
        [varargout{1:nargout}] = warp(varargin{:});
    case 'compose'
        [varargout{1:nargout}] = compose(varargin{:});
    case 'identity'
        [varargout{1:nargout}] = identity(varargin{:});
    case 'translation'
        [varargout{1:nargout}] = translation(varargin{:});
    case 'linear'
        [varargout{1:nargout}] = linear(varargin{:});
    case 'affine'
        [varargout{1:nargout}] = affine(varargin{:});
    case 'mm2vox'
        [varargout{1:nargout}] = mm2vox(varargin{:});
    case 'transform'
        [varargout{1:nargout}] = transform(varargin{:});
    otherwise
        help warps
        error('Unknown function %s. Type ''help warps'' for help.', id)
end

% === Functions ===========================================================

% -------------------------------------------------------------------------
function y = identity(lat_dim, lat_vs)
% FORMAT y = warps('identity', lat_dim, (lat_vs))
% lat_dim - Dimensions of the lattice on which to compute the map
% lat_vs  - Voxel size of the lattice [default: 1 1 1]
%
% Generate the identity warp on a given lattice.
%__________________________________________________________________________
if nargin < 2, lat_vs = [1 1 1]; end

lat_dim = [lat_dim 1 1 1];
lat_dim = lat_dim(1:3);
lat_vs  = [lat_vs 1 1 1];
lat_vs  = lat_vs(1:3);

y = zeros([lat_dim 3], 'single');
[y(:,:,:,1), y(:,:,:,2), y(:,:,:,3)] = ...
    ndgrid(lat_vs(1) * single(1:lat_dim(1)), ...
           lat_vs(2) * single(1:lat_dim(2)), ...
           lat_vs(3) * single(1:lat_dim(3)));
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function y = translation(T, lat_dim, lat_vs)
% FORMAT y = warps('translation', T, lat_dim, (lat_vs))
% T       - [3 double] Translation
% lat_dim - Dimensions of the lattice on which to compute the map
% lat_vs  - Voxel size of the lattice [default: 1 1 1]
%
% Generate a translation warp on a given lattice.
%__________________________________________________________________________
if nargin < 3, lat_vs = [1 1 1]; end

lat_dim = [lat_dim 1 1 1];
lat_dim = lat_dim(1:3);
lat_vs  = [lat_vs 1 1 1];
lat_vs  = lat_vs(1:3);
T       = [T 0 0 0];
T       = T(1:3);

y = zeros([lat_dim 3], 'single');
[y(:,:,:,1), y(:,:,:,2), y(:,:,:,3)] = ...
    ndgrid(lat_vs(1) * single(1:lat_dim(1)) + T(1), ...
           lat_vs(2) * single(1:lat_dim(2)) + T(2), ...
           lat_vs(3) * single(1:lat_dim(3)) + T(3));
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function y = linear(L, lat_dim, lat_vs)
% FORMAT y = warps('linear', L, lat_dim, (lat_vs))
% L       - [3x3 double] Linear transform
% lat_dim - Dimensions of the lattice on which to compute the map
% lat_vs  - Voxel size of the lattice [default: 1 1 1]
%
% Generate a linear warp on a given lattice.
%__________________________________________________________________________
if nargin < 3, lat_vs = [1 1 1]; end

lat_dim = [lat_dim 1 1 1];
lat_dim = lat_dim(1:3);
lat_vs  = [lat_vs 1 1 1];
lat_vs  = lat_vs(1:3);
dL = size(L);
L2 = L;
L  = eye(3);
L(1:dL(1), 1:dL(2)) = L2;

y = L * reshape(identity(lat_dim, lat_vs), [], 3)';
y = reshape(y', lat_dim, 3);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function y = affine(A, lat_dim, lat_vs)
% FORMAT y = warps('affine', A, lat_dim, (lat_vs))
% A       - [4x4 double] Affine transform
% lat_dim - Dimensions of the lattice on which to compute the map
% lat_vs  - Voxel size of the lattice [default: 1 1 1]
%
% Generate an affine warp on a given lattice.
%__________________________________________________________________________
if nargin < 3, lat_vs = [1 1 1]; end

lat_dim = [lat_dim 1 1 1];
lat_dim = lat_dim(1:3);
lat_vs  = [lat_vs 1 1 1];
lat_vs  = lat_vs(1:3);
dA = size(A);
A2 = A;
A  = eye(4);
A(1:dA(1), 1:dA(2)) = A2;

y = A(1:3,1:3) * reshape(identity(lat_dim, lat_vs), [], 3)';
y = bsxfun(@plus, y, A(1:3,4));
y = reshape(y', [lat_dim 3]);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function y = compose(varargin)
% FORMAT y = warps('compose', y_1, (vs_1), ..., y_n, (vs_n), (itrp))
% y_i  - A warp OR an affine matrix.
% vs_i - Voxel size of the lattice (for warps only) [default: 1 1 1]
% itrp - Interpolation degree [default: 1]
% y    - A warp with the same lattice and voxel size as y_n.
%
% NB:
% - The right-most transform should always be a warp.
% - To specify the output lattice, add an identity warp as the
%   right-most argument.
%
% The argument parsing scheme is the following:
% - If the number of (ns) dimensions is 1       : Voxel size
% - If the number of (ns) dimensions is 2       : Affine matrix
% - If the number of (ns) dimensions is >= 3    : Warp
% - If the last argument is scalar              : Interpolation order
%
% Compose a series of transformations.
%__________________________________________________________________________

if nargin == 0
    y = [];
    return
end

if ~isempty(varargin) && isscalar(varargin{end})
    itrp     = varargin{end};
    varargin = varargin(1:end-1);
else
    itrp = 1;
end

if isempty(varargin)
    y = [];
    return
end

% --- Initialise
if isvector(varargin{end})
    vs = varargin{end};
    y = varargin{end-1};
    varargin = varargin(1:end-2);
else
    vs = [1 1 1];
    y = varargin{end};
    varargin = varargin(1:end-1);
end

% --- Loop
while ~isempty(varargin)
    if isvector(varargin{end})
        cur_vs = varargin{end};
        cur_y = varargin{end-1};
        varargin = varargin(1:end-2);
    else
        cur_vs = [1 1 1];
        cur_y = varargin{end};
        varargin = varargin(1:end-1);
    end
    if length(size(cur_y)) == 2
        % Affine o Warp
        y = transform(cur_y, y);
    else
        % Warp o Warp
        % > Warp the deformation field with circulant boundary
        cur_lat = size(cur_y);
        cur_lat = cur_lat(1:3);
        id = identity(cur_lat, cur_vs);
        y = warp(cur_y - id, y, cur_vs, itrp, 1) + y;
    end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function out = warp(in, y, vs_in, itrp, bnd)
% FORMAT out = warps(('warp'), in, y, (vs_in), (itrp), (bnd))
% in    - Input image (or function R^3 -> R^d).
% y     - Non-linear warp.
% vs_in - Voxel size of the input image lattice [default: 1 1 1].
% itrp  - Interpolation order [default: 1 1 1].
% bnd   - Boundary conditions (0/1 = mirror/circulant) [default: 1 1 1]
%
% Warps an image with a non-linear transform, i.e., computes in(y).
% The input image can be non-scalar.
%__________________________________________________________________________
if nargin < 5, bnd   = [1 1 1]; end
if nargin < 4, itrp  = [1 1 1]; end
if nargin < 3, vs_in = [1 1 1]; end

if numel(itrp) < 3
    itrp = padarray(itrp, [0 3 - numel(itrp)], 'replicate', 'post');
end
if numel(bnd) < 3
    bnd = padarray(bnd, [0 3 - numel(bnd)], 'replicate', 'post');
end
if numel(vs_in) < 3
    vs_in = padarray(vs_in, [0 3 - numel(vs_in)], 'replicate', 'post');
end

dim_in = size(in);
in = reshape(in, dim_in(1), dim_in(2), dim_in(3), []);
dim_out = size(y);
dim_out = dim_out(1:3);

out = zeros([dim_out size(in, 4)], 'like', in);
for k=1:size(in, 4)
    out(:,:,:,k) = warp_scalar(in(:,:,:,k), y, vs_in, itrp, bnd);
end
out = reshape(out, [dim_out dim_in(4:end)]);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function y = mm2vox(y, vs)
% FORMAT y = warps('mm2vox', y, vs)
% y  - Non-linear warp
% vs - Voxel size of the target lattice
%
% /!\ vs is not the voxel size of the warp lattice !
%
% Transform millimetric warps into voxel warps.
% Target coordinates are divided by the voxel size of the target image.
%__________________________________________________________________________
dim = size(y);
y = reshape(y, [], 3);
y = bsxfun(@times, y', vs(:));
y = reshape(y', dim);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function y = transform(A, y)
% FORMAT y = warps('transform', A, y)
% A - Affine matrix
% y - Non-linear warp
%
% Affine transform a warp.
% Note that you can obtain the same result by doing warp('compose', A, y).
%__________________________________________________________________________
dim = size(y);
dim = dim(1:3);

y = A(1:3,1:3) * reshape(y, [], 3)';
y = bsxfun(@plus, y, A(1:3,4));
y = reshape(y', [dim 3]);
% -------------------------------------------------------------------------

% === Helpers =============================================================

% -------------------------------------------------------------------------
function out = warp_scalar(in, y, vs_in, itrp, bnd)
% FORMAT out = warp_scalar(in, y, (vs_in), (itrp))
% in    - Input _scalar_ image (or function R^3 -> R).
% y     - Non-linear warp.
% vs_in - Voxel size of the input image lattice..
% itrp  - Interpolation order.
% bnd   - Boundary conditions (0/1 = mirror/circulant).
%
% Warps an image with a non-linear transform, i.e., computes in(y).
% The input image must be scalar.
%__________________________________________________________________________
% Interpolate input image with circulant boundaries
in_coeff = spm_diffeo('bsplinc', single(in), [itrp bnd]);
% Convert coordinates from mm to voxels
y = mm2vox(y, vs_in);
% Interpolate image on output grid
out = spm_diffeo('bsplins', in_coeff, single(y), [itrp bnd]);
% -------------------------------------------------------------------------

