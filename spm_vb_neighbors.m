function [vxyz] = spm_vb_neighbors (xyz,vol)
% Create list of neighbors of voxels to be analysed
% FORMAT [vxyz] = spm_vb_neighbors (xyz,vol)
%
% xyz      [Nvoxels x 3] list of voxel positions which are to be analysed
% vol      vol=1 for volumetric neighbors, vol=0 for within-slice neighbors 
%          (default vol=0)
%
% vxyz     [Nvoxels x 4] list of neighbouring voxels
%          or [Nvoxels x 6] list of neighbouring voxels for vol=1
%
%          vxyz(j,:)=[N1 N2 N3 0] means that there are only 3 neighbors
%          of voxel j, and their numbers (ie. where they appear in the xyz list) 
%          are N1, N2 and N3
%
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_neighbors.m 112 2005-05-04 18:20:52Z john $

if nargin<2
    vol=0;
end

[voxels,Ndims]=size(xyz);

if vol
    vxyz=zeros(voxels,6);
else
    vxyz=zeros(voxels,4);
end    

for i=1:voxels,
    Ni=0;
    for j=1:voxels,
        dx=abs(xyz(i,1)-xyz(j,1));
        dy=abs(xyz(i,2)-xyz(j,2));
        if vol
            dz=abs(xyz(i,3)-xyz(j,3));
        else
            dz=0;
        end
        if (dx+dy+dz==1)
            % j is a neighbor of i
            Ni=Ni+1;
            vxyz(i,Ni)=j;
        end
    end
end
    
