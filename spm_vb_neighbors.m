function [vxyz] = spm_vb_neighbors (xyz)
% Create list of neighbors of voxels to be analysed
% FORMAT [vxyz] = spm_vb_neighbors (xyz)
%
% xyz      [Nvoxels x 3] list of voxel positions which are to be analysed
%
% vxyz     [Nvoxels x 4] list of neighbouring voxels
%          vxyz(j,:)=[N1 N2 N3 0] means that there are only 3 neighbors
%          of voxel j, and their numbers (ie. where they appear in the xyz list) 
%          are N1, N2 and N3
%
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id$


[voxels,Ndims]=size(xyz);
% if voxels<Ndims
%     disp('Error in neighbours: xyz wrong way round');
%     return
% end

vxyz=zeros(voxels,4);

for i=1:voxels,
    Ni=0;
    for j=1:voxels,
        dx=abs(xyz(i,1)-xyz(j,1));
        dy=abs(xyz(i,2)-xyz(j,2));
        if (dx+dy==1)
            % j is a neighbor of i
            Ni=Ni+1;
            vxyz(i,Ni)=j;
        end
    end
end
    
