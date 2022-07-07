function X = spm_mesh_euler(M)
% Compute the Euler characteristic of a triangle mesh
% M        - patch structure
%
% X        - Euler characteristic
%__________________________________________________________________________
%
% The Euler characteristic is defined according to the formula:
%
%                    X  =  V - E + F  =  2 - 2g - b
%
% where g is the genus and b the number of boundary components.
% See https://www.wikipedia.org/wiki/Euler_characteristic
%     https://www.wikipedia.org/wiki/Genus_(mathematics)
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


X = size(M.vertices,1) - size(spm_mesh_edges(M),1) + size(M.faces,1);
