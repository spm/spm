function ksi = spm_mesh_euler(M)
% Compute the Euler characteristic of a triangle mesh
% M        - patch structure
%
% ksi      - Euler characteristic
%__________________________________________________________________________
%
% The Euler characteristic is defined according to the formula:
%                           \ksi = V - E + F
% See https://www.wikipedia.org/wiki/Euler_characteristic
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_euler.m 7397 2018-08-15 11:04:26Z guillaume $


ksi = size(M.vertices,1) - size(spm_mesh_edges(M),1) + size(M.faces,1);
