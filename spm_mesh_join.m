function [M, I] = spm_mesh_join(varargin)
% Join a list of surface meshes into a single one
% FORMAT [M, I] = spm_mesh_join(Ms)
% Ms            - a patch structure array or list of scalar patch structures
%
% M             - a scalar patch structure
% I             - a column vector of face indices
%
% See also spm_mesh_split
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


Ms = varargin{1};
for i=2:numel(varargin)
    fn = fieldnames(varargin{i});
    for j=1:numel(fn)
        Ms(i).(fn{j}) = varargin{i}.(fn{j});
    end
end
if isfield(Ms,'mat')
    for i=1:numel(Ms)
        if isempty(Ms(i).mat)
            Ms(i).mat = eye(4);
        end
    end
end

fn = fieldnames(Ms(1));
M  = cell2struct(cell(numel(fn),1),fn);
I  = zeros(numel(Ms),1);

for i=1:numel(Ms)
    if i==1
        M.faces   = Ms(i).faces;
    else
        M.faces   = [M.faces; Ms(i).faces+size(M.vertices,1)];
    end
    I(i)          = size(Ms(i).faces,1);
    M.vertices    = [M.vertices; Ms(i).vertices];
    try, M.cdata  = [M.cdata; Ms(i).cdata]; end
    try
        if i==1, M.mat = Ms(i).mat; end
    end
    if isfield(M,'mat')
        if sum(sum(M.mat - Ms(i).mat)) > 10*eps
            error('Meshes have different orientation.');
        end
    end
end
if nargout > 1, I = repelem((1:numel(I))',I); end
