function C = spm_mesh_label(M)
% Label connected components of a surface mesh
% FORMAT C = spm_mesh_label(M)
% M        - a [nx3] faces array or a patch structure
%
% C        - a [nx1] vector containing labels for the connected components
%            in M
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_label.m 3152 2009-05-27 10:54:49Z guillaume $

if ishandle(M)
    M = get(M,'Faces');
elseif ~isnumeric(M)
    M = M.faces;
end

A = spm_mesh_adjacency(M);
[p,q,r] = dmperm(A);

C = zeros(size(M,1),1);
for i=1:length(r)-1
    C(any(ismember(M,r(i):r(i+1)-1),2)) = i;
end
