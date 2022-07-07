function L = spm_mesh_get_lm(M,T)
% Identification of local maxima on a textured surface mesh
% FORMAT L = spm_mesh_get_lm(M,T)
% M        - a [nx3] faces array or a patch structure or a [nxn] adjacency
%            matrix
% T        - a [nx1] texture vector
%
% L        - indices of vertices that are local maxima
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


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
if all(out), L = []; return; end             % empty domain
A(:,out) = 0;
A(out,:) = 0;
N        = spm_mesh_neighbours(A);
if isempty(N), L = find(~out); return; end   % only singletons
%S       = ~any(N(~out,:),2);                % some singletons

%-Identify local maxima
%--------------------------------------------------------------------------
T        = [T; -Inf];
N(N<1)   = numel(T);
L        = all(bsxfun(@gt,T(1:end-1),T(N)),2);
L        = find(L');
