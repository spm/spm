function f = spm_pac_transform_mesh(f, varargin)
% FORMAT wf = transform_mesh(f, iy, [M])
%   f    - Input mesh (gifti object)
%   iy   - Transformation [Mx My Mz 3]
%   M    - Voxel-to-world matrix of the transformation {default: eye(4)}
%   wf   - Warped mesh (gifti object)
%
% FORMAT wf = transform_mesh(f, T, [iy], [M])
%   T    - World-to-world transformation to apply to the mesh first.
%
% . The transformation (iy) should be expressed in millimetres, that is,
%   each voxel contains a millimetric coordinate in world space.
% . The voxel-to-world matrix (M) should map voxels to mm.
% . The affine transformation matrix (T) should map mm to mm.

%% Yael Balbastre 2024

if nargin < 2
    error('At least two arguments required.');
end

% -------------------------------------------------------------------------
% Read arguments and select mode
if numel(size(varargin{1})) == 2 ...
        && size(varargin{1},1) == 4 ...
        && size(varargin{1},2) == 4
    T = varargin{1};
    if numel(varargin) > 1
        iy = varargin{2};
        if numel(varargin) > 2
            M = varargin{3};
        else
            M = eye(4);
        end
    else
        iy = [];
    end
else
    T = [];
    iy = varargin{1};
    if numel(varargin) > 1
        M = varargin{2};
    else
        M = eye(4);
    end
end


% -------------------------------------------------------------------------
% Apply first affine (worl to world) matrix
if ~isempty(T)
    f = spm_mesh_transform(f, T);
end

% -------------------------------------------------------------------------
% Apply non-linear transformation
if ~isempty(iy)
    dimy = size(iy);
    if numel(dimy) == 4
        dimy = [dimy(1:3) 1 dimy(4)];
        iy   = reshape(iy, dimy);
    end
    f  = spm_swarp(f, double(iy), M);
end