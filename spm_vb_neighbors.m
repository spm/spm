function vxyz = spm_vb_neighbors(xyz,vol)
% Create list of neighbors of voxels to be analysed
% FORMAT vxyz = spm_vb_neighbors (xyz,vol)
%
% xyz    - [Nvoxels x 3] list of voxel positions which are to be analysed
% vol    - vol=1 for volumetric neighbors, vol=0 for within-slice neighbors 
%          (default vol=0)
%
% vxyz   - [Nvoxels x 4] list of neighbouring voxels
%          or [Nvoxels x 6] list of neighbouring voxels for vol=1
%
%          vxyz(j,:)=[N1 N2 N3 0] means that there are only 3 neighbors
%          of voxel j, and their numbers (ie. where they appear in the xyz
%          list) are N1, N2 and N3
%
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny, Nelson Trujillo-Barreto and Lee Harrison
% $Id: spm_vb_neighbors.m 6079 2014-06-30 18:25:37Z spm $

if nargin<2
    vol = 0;
end

nearestneighbor = 1;
N               = size(xyz,1);
DIM             = max(xyz,[],1);
if vol,
    I               = sub2ind([DIM(1),DIM(2),DIM(3)],xyz(:,1),xyz(:,2),xyz(:,3));
    is              = zeros(DIM(1),DIM(2),DIM(3));
    is(I)           = 1:N; % index voxels in cluster
    [mm,nn,oo]      = deal(DIM(1)+2,DIM(2)+2,DIM(3)+2); % pad is
    iS              = zeros(mm,nn,oo);
    for i = 1:DIM(3),
        iS(2:DIM(1)+1,2:DIM(2)+1,i+1) = ones(DIM(1),DIM(2));
    end
    Ip          = find(iS);
    iS(Ip)      = is; % padded is
    [xx,yy,zz]  = ndgrid(-1:1,-1:1,-1:1);
    du          = [xx(:) yy(:) zz(:)];
    d           = [mm, 1, mm*nn]';
    vxyz        = zeros(N,6);
else
    I           = sub2ind([DIM(1),DIM(2)],xyz(:,1),xyz(:,2));
    is          = zeros(DIM(1),DIM(2));
    is(I)       = 1:N; % index voxels in cluster
    [mm,nn]     = deal(DIM(1)+2,DIM(2)+2);
    iS          = zeros(mm,nn);
    iS(2:DIM(1)+1,2:DIM(2)+1) = ones(DIM(1),DIM(2));
    Ip          = find(iS);
    iS(Ip)      = is; % padded is
    [xx,yy]     = ndgrid(-1:1,-1:1);
    du          = [xx(:) yy(:)];
    d           = [1, mm]';
    vxyz        = zeros(N,4);
end
    
% Indices of interior points
p           = find(iS);

% stencil
ind         =   find(sum(abs(du),2)~=0); % remove (0,0,0)
du          =   du(ind,:);
if nearestneighbor
    ind     =   find(sum(abs(du),2)==1); % remove off-diagonals
    du      =   du(ind,:);
end
Nq          =   size(du,1);

% edge set, from node k to node n
for j = 1:Nq,
    duj         = du(j,:)'; % vector of displacement on graph
    in          = d'*duj; % displacemnet index
    Q           = iS(p+in); % neighbouring nodes, zero if no neighbour
    q           = find(Q); % location of neighbours 
    kj          = iS(p(q)); % all nodes in p with neighbours
    nj          = Q(q); % all neighbours of p in direction duj
    vxyz(kj,j)  = nj;
end
