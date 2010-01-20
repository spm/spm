function I = spm_grid(I)
% Superimposes a Talairach and Tournoux grid
% FORMAT O = spm_grid(I);
% I - image matrix
% O - image matrix with grid added
%___________________________________________________________________________
%
% spm_grid adds a grid to the input argument.  The grid is scaled
% to 10% of the iuput's maximum.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_grid.m 3691 2010-01-20 17:08:30Z guillaume $

GRID = 0.1;

load('Grid.mat','i','j');
[x y] = size(I);
i     = round(1 + (i - 1)*(x - 1)/64);
j     = round(1 + (j - 1)*(y - 1)/86);
G     = full(sparse(i,j,max(I(:))*GRID*ones(length(i),1)));
I     = max(I,G);
