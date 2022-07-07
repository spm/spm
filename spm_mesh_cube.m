function M = spm_mesh_cube
% Triangle mesh of a unit cube
% FORMAT M = spm_mesh_cube
% M        - patch structure
%__________________________________________________________________________
%
% Return a triangle mesh of a unit cube (sides of 1 unit long).
% See https://www.wikipedia.org/wiki/Unit_cube
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


M.vertices = [...
    1 0 1
    1 0 0
    1 1 0
    1 1 1
    0 0 1
    0 0 0
    0 1 0
    0 1 1];

M.faces = [...
    5 1 4
    5 4 8
    1 2 3
    1 3 4
    2 6 7
    2 7 3
    6 5 8
    6 8 7
    8 4 3
    8 3 7
    1 6 2
    1 5 6];
