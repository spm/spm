function P = spm_mesh_project(M, dat, method, varargin)
% Project volumetric data onto a mesh
% FORMAT P = spm_mesh_project(M, dat, method)
% M        - a patch structure, a handle to a patch 
%            or a [nx3] vertices array
% dat      - a structure array [1xm] with fields dim, mat, XYZ and t 
%            (see spm_render.m)
% method   - interpolation method {'nn'}
% varargin - other parameters required by the interpolation method
%
% P        - a [mxn] curvature vector
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_project.m 3277 2009-07-15 11:47:40Z guillaume $

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

P = zeros(length(dat),size(V,1));
for i=1:length(dat)
    Y      = zeros(dat(i).dim(1:3)');
    OFF    = dat(i).XYZ(1,:) + dat(i).dim(1)*(dat(i).XYZ(2,:)-1 + dat(i).dim(2)*(dat(i).XYZ(3,:)-1));
    Y(OFF) = dat(i).t .* (dat(i).t > 0);
    XYZ    = double(inv(dat(i).mat)*[V';ones(1,size(V,1))]);
    P(i,:) = spm_sample_vol(Y,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
end
