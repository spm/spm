function E = spm_mesh_edges(M)
% Return edges of a surface mesh
% FORMAT E = spm_mesh_edges(M)
% M        - a [nx3] faces array or a patch handle/structure
%
% E        - a [mx2] edges array 
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_edges.m 4018 2010-07-27 18:22:42Z guillaume $


%-Parse input arguments
%--------------------------------------------------------------------------
if ishandle(M)
    M = get(M,'Faces');
elseif ~isnumeric(M)
    M = M.faces;
end

%-Compute edges
%--------------------------------------------------------------------------
M = sort(M,2);
E = unique([M(:,[1 2]);M(:,[2 3]);M(:,[1 3])],'rows');
