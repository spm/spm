function P = spm_mesh_project(M, dat, method, varargin)
% Project volumetric data onto a mesh
% FORMAT P = spm_mesh_project(M, dat, method)
% M        - a patch structure, a handle to a patch 
%            or a [nx3] vertices array
% dat      - a structure array [1xm] with fields dim, mat, XYZ and t 
%            (see spm_render.m)
%            or a structure array [1xm] with fields mat and dat
%            or a char array/cellstr of image filenames
% method   - interpolation method {'nn'}
% varargin - other parameters required by the interpolation method
%
% P        - a [mxn] curvature vector
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_project.m 6987 2017-01-13 16:07:00Z volkmar $

if ishandle(M)
    V = get(M,'Vertices');
elseif isstruct(M) || isa(M,'gifti')
    V = M.vertices;
else
    V = M;
end

if nargin < 3, method = 'nn'; end
if ~strcmpi(method,'nn')
    error('Only Nearest Neighbours interpolation is available.');
end

if ischar(dat), dat = cellstr(dat); end
P = zeros(length(dat),size(V,1));
for i=1:numel(dat)
    if iscellstr(dat)
        v      = spm_vol(dat{i});
        Y      = spm_read_vols(v);
        mat    = v.mat;
    elseif isfield(dat,'dat')
        Y      = dat(i).dat;
        mat    = dat(i).mat;
    else
        if isfield(dat,'Z') %-xSPM structure
           dat = struct('dim',dat.DIM,'XYZ',dat.XYZ,'t',dat.Z,'mat',dat.M);
        end
        dim    = dat(i).dim(1:3);
        Y      = zeros(dim(:)');
        OFF    = dat(i).XYZ(1,:) + ...
                 dat(i).dim(1)*(dat(i).XYZ(2,:)-1 + ...
                 dat(i).dim(2)*(dat(i).XYZ(3,:)-1));
        Y(OFF) = dat(i).t .* (dat(i).t > 0);
        mat    = dat(i).mat;
    end
    XYZ        = double(inv(mat)*[V';ones(1,size(V,1))]);
    P(i,:)     = spm_sample_vol(Y,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
end
