function spm_orthviews(V,centre,origin,bb,vox,area,ttle)
% Display Orthogonal Views of a Normalized Image
% FORMAT spm_orthviews(V,centre,origin,bb,vox,area,ttle)
% V       - the memory mapped volume.
% centre  - the coordinate in mm which appears in all images.
% origin  - the origin of the volume.
% bb      - bounding box.
% vox     - voxel size to sample.
% area    - the area of the figure to use.
%           area(1) - position x
%           area(2) - position y
%           area(3) - size x
%           area(4) - size y
% ttle    - the title
%__________________________________________________________________________
% %W% John Ashburner %E%

if (nargin < 7)
	ttle = '';
	if (nargin < 6)
		area = [0 0 1 1];
		if (nargin < 5)
			vox = 2;
			if (nargin < 4)
				bb =[ -64 -104 -28 ; 64 68 72];
			end
		end
	end
end


Dims = round(diff(bb)/vox);


TM = [
vox/V(4) 0        0        (origin(1)-1 + bb(1,1)/V(4))
0        vox/V(5) 0        (origin(2)-1 + bb(1,2)/V(5))
0        0        vox/V(6) (origin(3)-1-centre(3)/V(6))
0        0        0        1];

CM = [
vox/V(4) 0        0        (origin(1)-1 + bb(1,1)/V(4))
0        0        vox/V(5) (origin(2)-1-centre(2)/V(5))
0        vox/V(6) 0        (origin(3)-1 + bb(1,3)/V(6))
0        0        0        1];

SM = [
0        0       vox/V(4)  (origin(1)-1-centre(1)/V(4))
0        vox/V(5) 0        (origin(2)-1 + bb(1,2)/V(5))
vox/V(6) 0        0        (origin(3)-1 + bb(1,3)/V(6))
0        0        0        1];

tran  = spm_slice_vol(V(:,1),TM,Dims([1 2]),1);
coron = spm_slice_vol(V(:,1),CM,Dims([1 3]),1);
sag   = spm_slice_vol(V(:,1),SM,Dims([3 2]),1);

for i=2:size(V,2)
	tran  = tran  + spm_slice_vol(V(:,i),TM,Dims([1 2]),1);
	coron = coron + spm_slice_vol(V(:,i),CM,Dims([1 3]),1);
	sag   = sag   + spm_slice_vol(V(:,i),SM,Dims([3 2]),1);
end

un=get(gcf,'Units');
set(gcf,'Units','Pixels');
sz=get(gcf,'Position');
set(gcf,'Units',un);

sz = sz(3:4);
sz(2) = sz(2)-40;
area = [area(1)*sz(1) area(2)*sz(2) area(3)*sz(1) area(4)*sz(2)];

mx = max([max(max(tran)) max(max(sag)) max(max(coron))]);
sx = area(3)/(Dims(1)+Dims(3));
sy = area(4)/(Dims(2)+Dims(3));
s = min([sx sy]);
offx = (area(3)-(Dims(1)+Dims(3))*s)/2 + area(1);
offy = (area(4)-(Dims(2)+Dims(3))*s)/2 + area(2);

ax=axes('Units','pixels', ...
	'Position',[offx offy s*Dims(1) s*Dims(2)],'Units','normalized');
image(rot90(tran*(64/mx)));axis('image');
set(ax,'XGrid','on','YGrid','on','XColor',[1 0 0],'YColor',[1 0 0],...
	'XTickLabels',[],'YTickLabels',[]);

ax=axes('Units','pixels', ...
	'Position',[offx+s*Dims(1) offy s*Dims(3) s*Dims(2)],'Units','normalized');
image(rot90(sag*(64/mx)));axis('image');
set(ax,'XGrid','on','YGrid','on','XColor',[1 0 0],'YColor',[1 0 0],... 
	'XTickLabels',[],'YTickLabels',[]);

ax=axes('Units','pixels', ...
	'Position',[offx offy+s*Dims(2) s*Dims(1) s*Dims(3)],'Units','normalized');
image(rot90(coron*(64/mx)));axis('image');
set(ax,'XGrid','on','YGrid','on','XColor',[1 0 0],'YColor',[1 0 0],...  
	'XTickLabels',[],'YTickLabels',[]);

axes('Units','pixels', ...
	'Position',[offx+s*Dims(1) offy+s*Dims(2) s*Dims(1) s*Dims(3)], ...
	'Visible','off','Xlim',[1 s*Dims(1)],'Ylim',[1 s*Dims(3)],'Units','normalized');

skip = get(gca,'FontSize')+4;

for i=1:size(ttle,1)
	text(10,skip*(size(ttle,1)+2)-skip*i,ttle(i,:));
end
