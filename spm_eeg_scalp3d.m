function scalp3d(data)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% James Kilner
% $Id: spm_eeg_scalp3d.m 193 2005-07-01 13:37:59Z james $


[xs,ys,zs]=sphere(50);
xs=xs.*0.90;
ys=ys.*0.90;
zs=zs.*0.90;

P = fullfile(spm('dir'), 'EEGtemplates');
if length(data)==128
	load(fullfile(P, '3d_bdf'));
end
if length(data)==275
	load(fullfile(P, '3d_ctf'));
end

spm_defaults


map=griddata3(pos(:,1),pos(:,2),pos(:,3),data,xs,ys,zs,'nearest');


load(fullfile(P, 'params.mat'));

map=map(:);
map(aind)='';
C=map;
%[col_scalp]=spm_eeg_marry_scalp_head_quick(V,C,V3,elec_verts,dv,id);
[col_scalp]=spm_eeg_marry_scalp_head_quick(V,C,V3,elec_verts,dv,id);


figure
load(fullfile(P, 'head2'));

p = patch(datas);
set(p, 'FaceColor', [0.6 0.6 0.6] , 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
axis tight
%reducepatch(p,.15)
elec_pos=spm_eeg_elec_pos3d;
hold on

p2=patch(datas,'FaceVertexCData',col_scalp);

set(p2, 'FaceColor','interp','EdgeColor', 'none');
daspect([1 1 1])
view(3)
axis tight
plot3(elec_pos(:,1),elec_pos(:,2),elec_pos(:,3),'w.')

camlight
axis('off')
colorbar
