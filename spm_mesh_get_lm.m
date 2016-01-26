function L = spm_mesh_get_lm(M,T)
% Identification of local maxima on a textured surface mesh
% FORMAT L = spm_mesh_get_lm(M,T)
% M        - a [nx3] faces array or a patch structure or a [nxn] adjacency
%            matrix
% T        - a [nx1] texture vector
%
% L        - indices of vertices that are local maxima
%__________________________________________________________________________
% Copyright (C) 2010-2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_get_lm.m 6694 2016-01-26 17:09:11Z guillaume $


%-Get adjacency matrix
%--------------------------------------------------------------------------
if ~isnumeric(M) || size(M,1) ~= size(M,2)
    A = spm_mesh_adjacency(M);
else
    A = M;
end

%-Get neighbours, restricted to vertices that actually have data (~= NaN)
%--------------------------------------------------------------------------
out      = isnan(T(:)');
A(:,out) = 0;
A(out,:) = 0;
N        = spm_mesh_neighbours(A);

%-Identify local maxima
%--------------------------------------------------------------------------
T        = [T; -Inf];
N(N<1)   = numel(T);
L        = all(T(N) < repmat(T(1:end-1),1,size(N,2)),2); % bsxfun?
L        = find(L');
