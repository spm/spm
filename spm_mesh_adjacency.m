function A = spm_mesh_adjacency(F)
% Compute the adjacency matrix of a triangle mesh
% FORMAT A = spm_mesh_adjacency(F)
% F        - a [fx3] faces array or a patch structure
% 
% A        - adjacency matrix as a sparse [vxv] array
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_adjacency.m 4035 2010-08-05 18:54:32Z guillaume $

if ~isnumeric(F) && isfield(F,'vertices')
    N = size(F.vertices,1);
    F = double(F.faces);
    A = sparse([F(:,1); F(:,1); F(:,2); F(:,2); F(:,3); F(:,3)], ...
           [F(:,2); F(:,3); F(:,1); F(:,3); F(:,1); F(:,2)], 1, N, N);
else
    if isstruct(F), F = F.faces; end
    F = double(F);
    A = sparse([F(:,1); F(:,1); F(:,2); F(:,2); F(:,3); F(:,3)], ...
           [F(:,2); F(:,3); F(:,1); F(:,3); F(:,1); F(:,2)], 1);
end
       
A = double(A > 0);
