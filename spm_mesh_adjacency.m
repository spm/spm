function A = spm_mesh_adjacency(F)
% Compute the adjacency matrix of a triangle mesh
% FORMAT A = spm_mesh_adjacency(F)
% F        - a [nx3] faces array or a patch structure
% 
% A        - adjacency matrix as a sparse [nxn] array
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_adjacency.m 3135 2009-05-19 14:49:42Z guillaume $

if ~isnumeric(F)
    F = F.faces;
end

F = double(F);

A = sparse([F(:,1); F(:,1); F(:,2); F(:,2); F(:,3); F(:,3)], ...
           [F(:,2); F(:,3); F(:,1); F(:,3); F(:,1); F(:,2)], 1);
A = double(A > 0);