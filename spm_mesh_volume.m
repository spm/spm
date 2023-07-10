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

% v1 = M.vertices(M.faces(:,1),:);
% v2 = M.vertices(M.faces(:,2),:);
% v3 = M.vertices(M.faces(:,3),:);
% 
% B = (v1 + v2 + v3) / 3;
% 
% N = cross(v2-v1, v3-v1, 2);
% 
% V = sum(dot(B, N, 2)) / 6;
