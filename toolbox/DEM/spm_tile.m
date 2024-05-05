function [G,M,H] = spm_tile(L,d)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT [G,M,H] = spm_tile(L,d)
%--------------------------------------------------------------------------
% L  - locations
% d  - size (diameter) of tile in pixels
%
% G  - cell array of outcome indices
% M  - and their locations
% H  - cell array of weights (c.f., receptive field)
%
% This routine identifies overlapping groups of pixels, returning their
% mean locations and a cell array of weights (based upon radial Gaussian
% basis functions) for each group. In other words, the grouping is based
% upon the location of pixels; in the spirit of a receptive field afforded
% by a sensory epithelium.Effectively, this leverages the conditional
% independencies that inherit from local interactions; of the kind found in
% metric spaces that preclude action at a distance.
%--------------------------------------------------------------------------

% centroid locations
%--------------------------------------------------------------------------
Ni    = round((max(L(:,1)) - min(L(:,1)) + 1)/d);
Nj    = round((max(L(:,2)) - min(L(:,2)) + 1)/d);
x     = linspace(min(L(:,1)) + d/2 - 1,max(L(:,1)) - d/2,Ni);
y     = linspace(min(L(:,2)) + d/2 - 1,max(L(:,2)) - d/2,Nj);

% groups of pixels within 2 x d of centroids
%--------------------------------------------------------------------------
h     = cell(Ni,Nj);
g     = cell(Ni,Nj);
for i = 1:Ni
    for j = 1:Nj
        D      = sum( minus(L,[x(i),y(j)]).^2,2 );
        ij     = find(sqrt(D) < 2*d);
        h{i,j} = exp(-D/(2*(d/2)^2));
        g{i,j} = ij;
    end
end
G    = g(:);

% weighting of groups (applying sum to one contraint)
%--------------------------------------------------------------------------
Ng    = numel(G);
h     = spm_cat(h(:)');
h     = spm_dir_norm(h');
H     = cell(Ng,1);
for g = 1:Ng
    H{g} = h(g,G{g});
end

% and mean location of groups
%--------------------------------------------------------------------------
for g = 1:Ng
    M(g,:) = mean(L(G{g},:));
end

return
