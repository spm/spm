function spm_surf(P,mode)
% Surface extraction.
% FORMAT spm_surf
%
% This surface extraction is not particularly sophisticated.  It simply
% smooths the data slightly and extracts the surface at a threshold of
% 0.5.
%
% Inputs:
% xxx_seg1.img & xxx_seg2.img - grey and white matter segments created
% using the segmentation routine.  These can be manually cleaned up
% first using e.g., MRIcro.
%
% Outputs:
% A "render_xxx.mat" file can be produced that can be used for
% rendering activations on to.
%
% A "surf_xxx.mat" file can also be written, which is created using
% Matlab's isosurface function.
% This extracted brain surface can be viewed using code something
% like:
%    FV = load(spm_get(1,'surf_*.mat','Select surface data'));
%    fg = spm_figure('GetWin','Graphics');
%    ax = axes('Parent',fg);
%    p  = patch(FV, 'Parent',ax,...
%           'FaceColor', [0.8 0.7 0.7], 'FaceVertexCData', [],...
%           'EdgeColor', 'none',...
%	    'FaceLighting', 'phong',...
%           'SpecularStrength' ,0.7, 'AmbientStrength', 0.1,...
%	    'DiffuseStrength', 0.7, 'SpecularExponent', 10);
%    set(0,'CurrentFigure',fg);
%    set(fg,'CurrentAxes',ax);
%    l  = camlight(-40, 20); 
%    axis image;
%    rotate3d on;
%
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin==0,
	SPMid = spm('FnBanner',mfilename,'%I%');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Surface');
	spm_help('!ContextHelp',mfilename);

	P    = spm_get([1 Inf],'IMAGE','Select images');

	mode = spm_input('Save','+1','m',...
		['Save Rendering|Save Extracted Surface|'...
		 'Save Rendering and Surface|Save Surface as OBJ format'],[1 2 3 4],3);
end;

spm('FigName','Surface: working',Finter,CmdLine);
do_it(P,mode);
spm('FigName','Surface: done',Finter,CmdLine);

return;
%_______________________________________________________________________

%_______________________________________________________________________
function do_it(P,mode)

spm('Pointer','Watch')

V  = spm_vol(P);
br = zeros(V(1).dim(1:3));
for i=1:V(1).dim(3),
	B         = spm_matrix([0 0 i]);
	tmp       = spm_slice_vol(V(1),B,V(1).dim(1:2),1);
	for j=2:length(V),
		M   = V(j).mat\V(1).mat*B;
		tmp = tmp + spm_slice_vol(V(j),M,V(1).dim(1:2),1);
	end;
	br(:,:,i) = tmp;
end;

% Build a 3x3x3 seperable smoothing kernel and smooth
%-----------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;
spm_conv_vol(br,br,kx,ky,kz,-[1 1 1]);

[pth,nam,ext] = fileparts(V(1).fname);

if any(mode==[1 3]),
	% Produce rendering
	%-----------------------------------------------------------------------
	matname = fullfile(pth,['render_' nam '.mat']);
	tmp     = struct('dat',br,'dim',size(br),'mat',V(1).mat);
	renviews(tmp,matname);
end;

if any(mode==[2 3 4]),
	% Produce extracted surface
	%-----------------------------------------------------------------------
	tmp     = struct('dat',br,'dim',size(br),'mat',V(1).mat);

	[faces,vertices] = isosurface(br,0.5);

	% Swap around x and y because isosurface does for some
	% wierd and wonderful reason.
	Mat      = V(1).mat(1:3,:)*[0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
	vertices = (Mat*[vertices' ; ones(1,size(vertices,1))])';
	if any(mode==[2 3]),
		matname = fullfile(pth,['surf_' nam '.mat']);
		save(matname,'faces','vertices');
	if any(mode==[4]),
		fname = fullfile(pth,['surf_' nam '.obj']);
		fid   = fopen(fname,'w');
		%for i=1:size(v,1),
		%	fprintf(fid,'v %f %f %f\n',vertices(i,1),vertices(i,2),vertices(i,3));
		%end;
		fprintf(fid,'v %f %f %f\n',vertices');
		fprintf(fid,'g Cortex\n'); % Group Cortex
		%for i=1:size(f,1),
		%	fprintf(fid,'f %d %d %d\n',faces(i,1),faces(i,2),faces(i,3));
		%end;
		fprintf(fid,'f %d %d %d\n',faces');
		fprintf(fid,'g\n');
		fclose(fid);
	end;

spm('Pointer')

return;
%_______________________________________________________________________

%_______________________________________________________________________
function renviews(V,oname)
% Produce images for rendering activations to
%
% FORMAT renviews(V,oname)
% V     - mapped image to render, or alternatively
%         a structure of:
%         V.dat - 3D array
%         V.dim - size of 3D array
%         V.mat - affine mapping from voxels to millimeters
% oname - the name of the render.mat file.
%_______________________________________________________________________
%
% Produces a matrix file "render_xxx.mat" which contains everything that
% "spm_render" is likely to need.
%
% Ideally, the input image should contain values in the range of zero
% and one, and be smoothed slightly.  A threshold of 0.5 is used to
% distinguish brain from non-brain.
%_______________________________________________________________________

linfun = inline('fprintf([''%-30s%s''],x,[sprintf(''\b'')*ones(1,30)])','x');
linfun('Rendering: ');

linfun('Rendering: Transverse 1..');    rend{1} = make_struct(V,[pi 0 pi/2]);
linfun('Rendering: Transverse 2..');    rend{2} = make_struct(V,[0 0 pi/2]);
linfun('Rendering: Saggital 1..');      rend{3} = make_struct(V,[0 pi/2 pi]);
linfun('Rendering: Saggital 2..');      rend{4} = make_struct(V,[0 pi/2 0]);
linfun('Rendering: Coronal 1..');       rend{5} = make_struct(V,[pi/2 pi/2 0]);
linfun('Rendering: Coronal 2..');       rend{6} = make_struct(V,[pi/2 pi/2 pi]);

linfun('Rendering: Save..');
save(oname,'rend');
linfun('                 ');
disp_renderings(rend);
spm_print;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function str = make_struct(V,thetas)
	[D,M]     = matdim(V.dim(1:3),V.mat,thetas);
	[ren,dep] = make_pic(V,M*V.mat,D);
	str       = struct('M',M,'ren',ren,'dep',dep);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [ren,zbuf]=make_pic(V,M,D)
% A bit of a hack to try and make spm_render_vol produce some slightly
% prettier output.  It kind of works...
if isfield(V,'dat'), vv = V.dat; else, vv = V; end;
[REN, zbuf, X, Y, Z] = spm_render_vol(vv, M, D, [0.5 1]);
fw        = max(sqrt(sum(M(1:3,1:3).^2)));
msk       = find(zbuf==1024);
brn       = ones(size(X));
brn(msk)  = 0;
brn       = spm_conv(brn,fw);
X(msk)    = 0;
Y(msk)    = 0;
Z(msk)    = 0;
msk       = find(brn<0.5);
tmp       = brn;
tmp(msk)  = 100000;
sX        = spm_conv(X,fw)./tmp;
sY        = spm_conv(Y,fw)./tmp;
sZ        = spm_conv(Z,fw)./tmp;
zbuf      = spm_conv(zbuf,fw)./tmp;
zbuf(msk) = 1024;

vec       = [-1 1 3]; % The direction of the lighting.
vec       = vec/norm(vec);
[t,dx,dy,dz] = spm_sample_vol(vv,sX,sY,sZ,3);
IM        = inv(diag([0.5 0.5 1])*M(1:3,1:3))';
ren       = IM(1:3,1:3)*[dx(:)' ; dy(:)' ; dz(:)'];
len       = sqrt(sum(ren.^2,1))+eps;
ren       = [ren(1,:)./len ; ren(2,:)./len ; ren(3,:)./len];
ren       = reshape(vec*ren,[size(dx) 1]);
ren(find(ren<0))  = 0;
ren(msk)  = ren(msk)-0.2;
ren       = ren*0.8+0.2;
mx        = max(ren(:));
ren       = ren/mx;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function disp_renderings(rend)
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
hght = 0.95;
nrow = ceil(length(rend)/2);
ax=axes('Parent',Fgraph,'units','normalized','Position',[0, 0, 1, hght],'Visible','off');
image(0,'Parent',ax);
set(ax,'YTick',[],'XTick',[]);

for i=1:length(rend),
	ren = rend{i}.ren;
	ax=axes('Parent',Fgraph,'units','normalized',...
		'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
		'Visible','off');
	image(ren*64,'Parent',ax);
	set(ax,'DataAspectRatio',[1 1 1], ...
		'PlotBoxAspectRatioMode','auto',...
		'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
end;
drawnow;
return;
%_______________________________________________________________________
function [d,M] = matdim(dim,mat,thetas)
R = spm_matrix([0 0 0 thetas]);
bb = [[1 1 1];dim(1:3)];
c  = [  bb(1,1) bb(1,2) bb(1,3) 1
        bb(1,1) bb(1,2) bb(2,3) 1
        bb(1,1) bb(2,2) bb(1,3) 1
        bb(1,1) bb(2,2) bb(2,3) 1
        bb(2,1) bb(1,2) bb(1,3) 1
        bb(2,1) bb(1,2) bb(2,3) 1
        bb(2,1) bb(2,2) bb(1,3) 1
        bb(2,1) bb(2,2) bb(2,3) 1]';
tc = diag([2 2 1 1])*R*mat*c;
tc = tc(1:3,:)';
mx = max(tc);
mn = min(tc);
M  = spm_matrix(-[mn(1) mn(2)])*diag([2 2 1 1])*R;
d  = ceil(abs(mx(1:2)-mn(1:2)))+1;
return;
