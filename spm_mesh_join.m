function M = spm_mesh_join(Ms)
% Join a list of surface meshes into a single one
% FORMAT M = spm_mesh_join(Ms)
% Ms       - a patch structure array
%
% M        - a patch structure
%
% See also spm_mesh_split
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_join.m 6912 2016-10-26 14:50:53Z guillaume $


fn = fieldnames(Ms(1));
M  = cell2struct(cell(numel(fn),1),fn);

for i=1:numel(Ms)
    if i==1
        M.faces   = Ms(i).faces;
    else
        M.faces   = [M.faces; Ms(i).faces+size(M.vertices,1)];
    end
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
