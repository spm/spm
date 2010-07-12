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
% $Id: spm_mesh_label.m 3989 2010-07-12 18:39:16Z guillaume $

if ishandle(M)
    M = get(M,'Faces');
elseif ~isnumeric(M)
    M = M.faces;
end

A = spm_mesh_adjacency(M);
A = A + speye(size(A));
[p,q,r] = dmperm(A);

%-Label faces
C = zeros(size(M,1),1);
for i=1:length(r)-1
    C(any(ismember(M,p(r(i):r(i+1)-1)),2)) = i;
end

%-Label vertices
% C = zeros(size(A,1),1);
% for i=1:length(r)-1
%     C(p(r(i):r(i+1)-1)) = i;
% end
