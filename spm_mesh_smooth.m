function T = spm_mesh_smooth(M, T, S)
% Perform Gaussian smoothing on data lying on a surface mesh
% FORMAT GL = spm_mesh_smooth(M)
% M        - a patch structure or a handle to a patch
% GL       - graph Laplacian
%
% FORMAT T = spm_mesh_smooth(M, T, S)
% FORMAT T = spm_mesh_smooth(GL, T, S)
% T        - [vx1] data vector
% S        - smoothing parameter (number of iterations)
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Guillaume Flandin
% $Id: spm_mesh_smooth.m 4079 2010-10-07 11:41:54Z guillaume $

if isstruct(M) || numel(M) == 1
    A  = spm_mesh_distmtx(M,0);
    N  = size(A,1);
    GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/16;
else
    GL = M;
end

if nargin == 1, T = GL; return; end

for i=1:S
    T = GL * T;
end

return;
