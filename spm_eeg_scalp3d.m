function scalp3d(data)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% James Kilner
% $Id$


[xs,ys,zs]=sphere(50);
xs=xs.*0.90;
ys=ys.*0.90;
zs=zs.*0.90;

P = fullfile(spm('dir'), 'EEGtemplates');

load(fullfile(P, '3d_bdf'));

spm_defaults
%% Plot first scalp map
%data=[ones(32,1);ones(32,1).*2;ones(32,1).*3;ones(32,1).*4];

map=griddata3(pos(:,1),pos(:,2),pos(:,3),data,xs,ys,zs,'nearest');


% [xs1,ys1,zs1]=sphere(50);
% [rzi,rzj]=find(zs1<-.55);
% xs1(rzi,rzj)=0;
% ys1(rzi,rzj)=0;
% zs1(rzi,rzj)=0;
% 
% [rzi,rzj]=find(zs1<-.2);
% [ryi,ryj]=find(ys1(rzi,rzj)>-0.4);
% 
% xs1(rzi(ryi),rzj(ryj))=0;
% ys1(rzi(ryi),rzj(ryj))=0;
% zs1(rzi(ryi),rzj(ryj))=0;

load(fullfile(P, 'params.mat'));
%[F,V,C]=make_faces(xs1,ys1,zs1,map);
map=map(:);
map(aind)='';
C=map;
[col_scalp]=spm_eeg_marry_scalp_head_quick(V,C,V3,elec_verts,dv,id);
%[V3,col_scalp,elec_pos]=marry_scalp_head(V,C,pos);

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
