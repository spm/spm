function P = spm_mesh_project(M, dat, method, varargin)
% Project volumetric data onto a mesh
% FORMAT P = spm_mesh_project(M, dat, method)
% M        - a patch structure, a handle to a patch 
%            or a [nx3] vertices array
% dat      - a structure with fields dim, mat, XYZ and t (see spm_render.m)
% method   - interpolation method {'nn'}
% varargin - other parameters required by the interpolation method
%
% P        - a [nx1] curvature vector
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_project.m 3206 2009-06-16 13:06:44Z guillaume $

if ishandle(M)
    V = get(M,'Vertices');
elseif isstruct(M)
    V = M.vertices;
else
    V = M;
end

if nargin < 3, method = 'nn'; end
if ~strcmpi(method,'nn')
    error('Only Nearest Neighbours interpolation is available.');
end

Y      = zeros(dat.dim(1:3)');
OFF    = dat.XYZ(1,:) + dat.dim(1)*(dat.XYZ(2,:)-1 + dat.dim(2)*(dat.XYZ(3,:)-1));
Y(OFF) = dat.t .* (dat.t > 0);
XYZ    = double(inv(dat.mat)*[V';ones(1,size(V,1))]);
P      = spm_sample_vol(Y,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
