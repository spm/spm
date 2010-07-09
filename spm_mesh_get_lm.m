function L = spm_mesh_get_lm(M,T)
% Identification of local maxima on a textured surface mesh
% FORMAT L = spm_mesh_get_lm(M,T)
% M        - a [nx3] faces array or a patch structure or a [nxn] adjacency
%            matrix
% T        - a [nx1] texture vector
%
% L        - indices of vertices that are local maxima
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_get_lm.m 3982 2010-07-09 13:36:53Z guillaume $

%-Obtain the adjacency matrix
%--------------------------------------------------------------------------
if ~isnumeric(M) || size(M,1) ~= size(M,2)
    A = spm_mesh_adjacency(M);
else
    A = M;
end

%-Restrict it to vertices that actually have data (~= NaN)
%--------------------------------------------------------------------------
out      = isnan(T(:)');
A(:,out) = 0; % A(out,:) is not necessary for usage below

%-Loop over vertices and identify local maxima
%--------------------------------------------------------------------------
L = [];

for i=find(~out)
    v = T(find(A(i,:)));
    if all(v<T(i)), L = [L i]; end
end
