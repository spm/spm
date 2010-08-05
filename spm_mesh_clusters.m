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
% $Id: spm_mesh_clusters.m 4035 2010-08-05 18:54:32Z guillaume $


%-Input parameters
%--------------------------------------------------------------------------
if ~islogical(T)
    T   = ~isnan(T);
end

%-Compute a reduced mesh corresponding to the data
%--------------------------------------------------------------------------
F       = spm_mesh_split(M, T);

%-Label connected components
%--------------------------------------------------------------------------
% will it find two connected vertices that do not form a triangle?
% -> maybe use spm_mesh_neighbours
[CC,N]  = spm_mesh_label(F, 'vertices');
C       = NaN(numel(T),1);
C(T)    = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
[N,ni]  = sort(N(:), 1, 'descend');
[ni,ni] = sort(ni);
C(T)    = ni(C(T));
