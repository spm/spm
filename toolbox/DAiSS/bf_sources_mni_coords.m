function res = bf_sources_mni_coords(BF, S)
% Generate beamforming grid
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Mark Woolrich
% $Id: bf_sources_mni_coords.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0 
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'Pos coords';        
    pos.val = {};
     
    mni_coords = cfg_branch;
    mni_coords.tag = 'mni_coords';
    mni_coords.name = 'Mni Coords';
    mni_coords.val = {pos};    
    
    res=mni_coords;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

% transform MNI coords in MNI space into space where we are doing the
% beamforming
M = inv(BF.data.transforms.toMNI);

grid.pos   = S.pos;

res = ft_transform_geometry(M, grid);

% establish index of nearest bilateral grid point
% for potential use in lateral beamformer.
res.bilateral_index=zeros(size(S.pos,1),1);
for jj=1:size(S.pos,1),
    mnic=S.pos(jj,:);
    mnic(1)=-mnic(1);
    res.bilateral_index(jj)=nearest_vec(S.pos,mnic);
end;

%pos=S.pos;figure;scatter3(pos(:,1),pos(:,2),pos(:,3),'.');
%hold on; jj=100; kk=res.bilateral_index(jj);
%scatter3(pos(jj,1),pos(jj,2),pos(jj,3),'og');
%scatter3(pos(kk,1),pos(kk,2),pos(kk,3),'or');
