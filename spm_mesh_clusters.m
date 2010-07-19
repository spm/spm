function [C, N] = spm_mesh_clusters(M,T)
% Label connected components of surface mesh data
% FORMAT [C, N] = spm_mesh_clusters(M,T)
% M        - a [mx3] faces array or a patch structure
% T        - a [nx1] data vector (using NaNs or logicals), n = #vertices
%
% C        - a [nx1] vector of cluster indices
% N        - a [px1] size of connected components {in vertices}
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_clusters.m 4003 2010-07-19 18:22:38Z guillaume $


%-Input parameters
%--------------------------------------------------------------------------
if isnumeric(M)
    F   = M;
else
    F   = M.faces;
end

%-Compute a reduced mesh corresponding to the data
%--------------------------------------------------------------------------
B       = NaN(size(T));
if islogical(T)
    i   = find(T);
else
    i   = find(~isnan(T));
end
B(i)    = 1:length(i);
F       = B(F);
F(any(isnan(F),2),:) = [];

%-Label connected components
%--------------------------------------------------------------------------
[CC,N]  = spm_mesh_label(F, 'vertices');
C       = NaN(numel(T),1);
C(i)    = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
[N,ni]  = sort(N(:), 1, 'descend');
[ni,ni] = sort(ni);
C(i)    = ni(C(i));
