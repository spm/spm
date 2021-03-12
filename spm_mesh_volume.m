function V = spm_mesh_volume(M)
% Compute the volume of a closed surface mesh
% FORMAT V = spm_mesh_volume(M)
% M        - a patch structure
% 
% V        - volume
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_volume.m 8078 2021-03-12 22:19:57Z guillaume $


V = spm_mesh_utils('volume',M);
