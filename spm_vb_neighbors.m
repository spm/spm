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
% %W% Will Penny and Nelson Trujillo-Barreto %E%

[voxels,Ndims]=size(xyz);
% if voxels<Ndims
%     disp('Error in neighbours: xyz wrong way round');
%     return
% end

x=kron(xyz(:,1),ones(voxels,1));
xx=repmat(xyz(:,1),voxels,1);
yz=nonzeros(abs(x-xx));
hx=min(yz);

y=kron(xyz(:,2),ones(voxels,1));
yy=repmat(xyz(:,2),voxels,1);
xz=nonzeros(abs(y-yy));
hy=min(yz);

vxyz=zeros(voxels,4);

for j=1:voxels,

    pxyz=repmat(xyz(j,:),voxels,1);
    dist=sqrt(sum((pxyz-xyz).^2,2));
    vx=find((dist<=hx)&(dist~=0));
    vy=find((dist<=hy)&(dist~=0));
    vxyzt=union(vx,vy);
    vxyz(j,1:length(vxyzt))=vxyzt';
    
end;

    