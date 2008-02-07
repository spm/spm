function I = spm_grid(I)
% Superimposes a Talairach and Tournoux grid
% FORMAT O = spm_grid(I);
% I - image matrix
% O - image matrix with grid added
%___________________________________________________________________________
%
% spm_grid adds a grid to the input argument.  The grid is scaled
% to defaults.grid times the imput's maximum, where defaults.grid is a user
% specified global variable
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_grid.m 1143 2008-02-07 19:33:33Z spm $


%---------------------------------------------------------------------------
global defaults
if ~isempty(defaults) & isfield(defaults,'grid'),
    GRID = defaults.grid/4;
else,
    return;
end;

if GRID
    load('Grid.mat','i','j');
    [x y] = size(I);
    i     = round(1 + (i - 1)*(x - 1)/64);
    j     = round(1 + (j - 1)*(y - 1)/86);
    G     = full(sparse(i,j,max(I(:))*GRID*ones(length(i),1)));
    I     = max(I,G);
end
