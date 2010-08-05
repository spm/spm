function MS = spm_mesh_split(M, C)
% Split a surface mesh into its connected components 
% FUNCTION [MS] = spm_mesh_split(M, C)
% M        - a [nx3] faces array or a patch structure
% C        - a [nx1] vector containing labels for the connected components
%            or a logical vector indicating vertices to keep
%
% MS       - a patch structure array
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_split.m 4035 2010-08-05 18:54:32Z guillaume $

if isnumeric(M), M = struct('faces',M); end

if nargin < 2
    [C, N] = spm_mesh_label(M,'vertices');
else
    if islogical(C)
        N  = nnz(C);
        %[C, N] = spm_mesh_clusters(M, C);
    else
        N  = histc(C,1:max(C));
    end
end

for i=1:numel(N)
    j                    = C == i;
    F                    = M.faces;
    B                    = NaN(size(C));
    B(j)                 = 1:N(i);
    F                    = B(F);
    F(any(isnan(F),2),:) = [];
    MS(i).faces          = F;
    try, MS(i).mat       = M.mat; end
    try, MS(i).cdata     = M.cdata(j,:); end
    try, MS(i).vertices  = M.vertices(j,:); end
end
