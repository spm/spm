function M = spm_mesh(varargin)
% Load mesh file(s) into memory as patch structure
% FORMAT M = spm_mesh(meshfilename1,meshfilename2,...)
%
% M        - patch structure array (.faces and .vertices) 
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh.m 8155 2021-09-26 16:29:44Z guillaume $


M = [];
if ~nargin
    [f, sts] = spm_select([1,Inf],'.*','Select mesh file');
    if ~sts, return; end
    f = cellstr(f);
else
    f = varargin;
end

for i=1:numel(f)
    M = [M, export(gifti(f{i}),'patch')];
end
