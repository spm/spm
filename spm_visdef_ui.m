function spm_visdef_ui
% Visualise deformations
%____________________________________________________________________________
% John Ashburner $Id$

P = spm_get(1,{'*y_*.img','noexpand'},'Select deformation field');
P = [repmat([deblank(P) ','],3,1) num2str([1 2 3]')];

bb = [-78 -112 0 ; 78 76 0];
directions = 'XYZ';
for d=1:3,
	str = sprintf('%d %d', bb(1,d), bb(2,d));
	bb(:,d) = sort(spm_input(['Bounding Box ' directions(d) ' (mm)'],....
		'+1', 'e',str, 2));
end;
bb = sort(bb);
vx = round((prod(diff(bb)+1)/100)^(1/3));
vx = spm_input('Approximate spacing (mm)', '+1', 'e',sprintf('%g',vx), 1);
vx = abs(vx);

V   = spm_vol(P); dbb = diff(bb);
vx1 = dbb(1)/ceil((dbb(1)+eps)/(vx+eps));
vx2 = dbb(2)/ceil((dbb(2)+eps)/(vx+eps));
vx3 = dbb(3)/ceil((dbb(3)+eps)/(vx+eps));

sc1 = min(10,ceil(2*vx1+eps));
sc2 = min(10,ceil(2*vx2+eps));
sc3 = min(10,ceil(2*vx3+eps));

dd1=vx1/sc1; if ~dd1, dd1=1; end;
dd2=vx2/sc2; if ~dd2, dd2=1; end;
dd3=vx3/sc3; if ~dd3, dd3=1; end;

[z1,z2,z3]=ndgrid(bb(1,1):dd1:bb(2,1),bb(1,2):dd2:bb(2,2),bb(1,3):dd3:bb(2,3));

M = inv(V(1).mat);
x1 = M(1,1)*z1+M(1,2)*z2+M(1,3)*z3+M(1,4);
x2 = M(2,1)*z1+M(2,2)*z2+M(2,3)*z3+M(2,4);
x3 = M(3,1)*z1+M(3,2)*z2+M(3,3)*z3+M(3,4);

d = [size(x1) 1 1];
y1 = spm_sample_vol(V(1),x1,x2,x3,[1 NaN]);y1=reshape(y1,d);
y2 = spm_sample_vol(V(2),x1,x2,x3,[1 NaN]);y2=reshape(y2,d);
y3 = spm_sample_vol(V(3),x1,x2,x3,[1 NaN]);y3=reshape(y3,d);

fig = spm_figure('GetWin','Graphics');
spm_clf
ax  = axes('Parent',fig);

if d(1)>1,
	t1 = y1(:,1:sc1:end,1:sc1:end);
	t2 = y2(:,1:sc1:end,1:sc1:end);
	t3 = y3(:,1:sc1:end,1:sc1:end);
	m  = [size(t1) 1];
	t1 = reshape(t1,m(1),m(2)*m(3));
	t2 = reshape(t2,m(1),m(2)*m(3));
	t3 = reshape(t3,m(1),m(2)*m(3));
	plot3(t1,t2,t3,'k','Parent',ax);
	hold on;
end;

if d(2)>1,
	t1 = y1(1:sc2:end,:,1:sc2:end);
	t2 = y2(1:sc2:end,:,1:sc2:end);
	t3 = y3(1:sc2:end,:,1:sc2:end);
	m  = [size(t1) 1];
	t1 = reshape(shiftdim(t1,1),m(2),m(1)*m(3));
	t2 = reshape(shiftdim(t2,1),m(2),m(1)*m(3));
	t3 = reshape(shiftdim(t3,1),m(2),m(1)*m(3));
	plot3(t1,t2,t3,'k');
	hold on;
end;

if d(3)>1,
	t1 = y1(1:sc3:end,1:sc3:end,:);
	t2 = y2(1:sc3:end,1:sc3:end,:);
	t3 = y3(1:sc3:end,1:sc3:end,:);
	m  = [size(t1) 1];
	t1 = reshape(shiftdim(t1,2),m(3),m(1)*m(2));
	t2 = reshape(shiftdim(t2,2),m(3),m(1)*m(2));
	t3 = reshape(shiftdim(t3,2),m(3),m(1)*m(2));
	plot3(t1,t2,t3,'k');
end;

hold off;
axis image xy off
rotate3d off; drawnow;
rotate3d on

