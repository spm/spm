function M = spm_mesh_isosurface(V, t, s)
% Compute isosurface geometry from volume data
% FORMAT M = spm_mesh_isosurface(V, t, s)
% V        - volume data
%            spm_vol struct, nifti object or 3D array
% t        - isosurface value
% s        - Gaussian filter width (FWHM) in {edges} [Default: 0]
%
% M        - patch structure
%
% This is merely a wrapper around isosurface.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_isosurface.m 7677 2019-10-24 09:55:24Z guillaume $


if ischar(V) || isstruct(V)
    V = spm_vol(V);
    V = struct('dat',spm_read_vols(V),'mat',V.mat);
elseif isa(V,'nifti')
    V = struct('dat',full(V.dat),'mat',V.mat);
elseif isnumeric(V) || islogical(V)
    V = struct('dat',V,'mat',eye(4));
else
    error('Invalid volume data type.');
end

if nargin < 3, s = 0; end

if any(s)
    spm_smooth(V.dat, V.dat, s);
end

for i=1:numel(t)
    
    [faces,vertices] = isosurface(V.dat, t(i));
    
    if isempty(vertices)
        faces    = zeros(0,3);
        vertices = zeros(0,3);
    end
    
    % Swap around x and y because isosurface uses meshgrid and not ndgrid
    mat      = V(1).mat(1:3,:) * [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
    vertices = (mat * [vertices'; ones(1,size(vertices,1))])';
    
    M(i) = struct('faces',faces, 'vertices',vertices);
    
end
