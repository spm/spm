function V = spm_mesh_volume(M)
% Compute the volume of a closed surface mesh
% FORMAT V = spm_mesh_volume(M)
% M        - a patch structure
% 
% V        - volume
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


V = spm_mesh_utils('volume',M);
