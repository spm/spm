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
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_euler.m 8081 2021-03-16 16:12:27Z guillaume $


X = size(M.vertices,1) - size(spm_mesh_edges(M),1) + size(M.faces,1);
