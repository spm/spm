function spm_xbrain
% Brain extraction.
% FORMAT spm_xbrain
%
% A rough and ready routine for extracting the brain from segmented
% images.  It begins by taking the white matter, and eroding it a
% couple of times to get rid of any odd voxels.  The algorithm
% continues on to do conditional dilations for several iterations,
% where the condition is based upon gray or white matter being present.
% The function has not yet been extensively validated and tested
% so it may not work every time.
%
% Inputs:
% xxx_seg1.img & xxx_seg2.img - grey and white matter segments extracted
% using spm_segment.
%
% Outputs:
% The extracted brain is written to "brain_xxx.img" in the same
% directory as xxx_seg1.img.
% A "render_xxx.mat" file can also be produced that can be used for
% rendering activations on to.
% In Matlab 5.3 and over, a "surf_xxx.mat" file can also be written.
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
%    box on;
%    rotate3d on;
%
%
% Note: that the images should be in a right handed co-ordinate system
% (same as the Talairach system) - otherwise the resulting renderings are
% reflections (mirror images) of the true images.
%
%_______________________________________________________________________
% %W% John Ashburner %E%

linfun = inline('fprintf([''%-60s%s''],x,[sprintf(''\b'')*ones(1,60)])','x');

SPMid = spm('FnBanner',mfilename,'%I%');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','XBrain');
spm_help('!ContextHelp','spm_xbrain.m');
P    = spm_get(2,'*_seg?.img','Select gray and white matter images');

v    = sscanf(version,'%f');
v    = v(1);
if v>= 5.3,
	mode = spm_input('Save','+1','m',...
		['Save Extracted Brain|Save Rendering|Save Extracted Surface|'...
		 'Save Extracted Brain and Rendering|'...
		 'Save Extracted Brain and Surface|'...
		 'Save Rendering and Surface|Save All'],[1 2 3 4 5 6 7],4);

else,
	mode = spm_input('Save','+1','m',...
		['Save Extracted Brain|Save Rendering|'...
		 'Save Extracted Brain and Rendering'],[1 2 4],3);
end;

spm('Pointer','Watch')
spm('FigName','Xbrain: working',Finter,CmdLine);


VG=spm_vol(P(1,:));	% Grey matter
VW=spm_vol(P(2,:));	% White matter

% Brain is initially set to white matter...
%-----------------------------------------------------------------------
br=zeros(VW.dim(1:3));
linfun('Initialising');
for i=1:VW.dim(3),
	br(:,:,i) = spm_slice_vol(VW,spm_matrix([0 0 i]),VW.dim(1:2),1);
end;

% Build a 3x3x3 seperable smoothing kernel
%-----------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Erosions and conditional dilations
%-----------------------------------------------------------------------
niter = 32;
spm_progress_bar('Init',niter,'Extracting Brain','Iterations completed');
for j=1:niter,
	if j>2, th=0.15; else, th=0.6; end; % Dilate after two its of erosion. 
	linfun(['Iteration ' num2str(j) ' - thresholding and multiplying']);
	for i=1:VW.dim(3),
		w=spm_slice_vol(VW,spm_matrix([0 0 i]),VW.dim(1:2),1);
		g=spm_slice_vol(VG,spm_matrix([0 0 i]),VW.dim(1:2),1);
		br(:,:,i) = (br(:,:,i)>th).*(w+g);
	end;
	linfun(['Iteration ' num2str(j) ' - convolving']);
	spm_conv_vol(br,br,kx,ky,kz,-[1 1 1]);
	spm_progress_bar('Set',j);
end;
spm_progress_bar('Clear');

% Generate output filename
%-----------------------------------------------------------------------
[pth,nam,ext,ver] = fileparts(deblank(P(1,:)));
ind = findstr(nam,'_seg1');
if ~isempty(ind), nam((0:4)+ind(1))=[]; end;
fname = fullfile(pth,['brain_' nam ext ver]);

if any(mode==[2 4 6 7]),
	% Produce rendering
	%-----------------------------------------------------------------------
	matname = fullfile(pth,['render_' nam '.mat']);
	tmp = struct('dat',br,'dim',size(br),'mat',VG.mat);
	renviews(tmp,matname);
end;

if any(mode==[3 5 6 7]),
	% Produce extracted surface
	%-----------------------------------------------------------------------
	linfun(['Extracting surface - please be patient']);
	matname = fullfile(pth,['surf_' nam '.mat']);
	tmp = struct('dat',br,'dim',size(br),'mat',VG.mat);

	% This is done with an eval statement because Matlab5.2
	% complains otherwise...
	eval('[faces,vertices] = isosurface(br,0.5);');

	% Swap around x and y because isosurface does for some
	% wierd and wonderful reason.
	Mat = VG.mat(1:3,:)*[0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
	vertices = (Mat*[vertices' ; ones(1,size(vertices,1))])';
	save(matname,'faces','vertices');
end;

if any(mode==[1 4 5 7]),
	% Final cleanup
	%-----------------------------------------------------------------------
	linfun(['Final cleanup - thresholding and multiplying']);
	for i=1:VW.dim(3),
		w  = spm_slice_vol(VW,spm_matrix([0 0 i]),VW.dim(1:2),1);
		g  = spm_slice_vol(VG,spm_matrix([0 0 i]),VW.dim(1:2),1);
		th = 0.15;
		br(:,:,i) = ((br(:,:,i)>th).*(w+g))>th;
	end;

	% Write the extracted brain
	%-----------------------------------------------------------------------
	linfun('Writing volume');
	VO=struct('fname',fname,...
		  'dim',[VW.dim(1:3) spm_type('uint8')],...
		  'mat',VW.mat,...
		  'pinfo',[1/255 0 0]',...
		  'descrip','extracted brain');
	spm_write_vol(VO,br);
end;

linfun(' ');
spm('FigName','Xbrain: done',Finter,CmdLine);
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
v      = V.dim(1:3).*sqrt(sum(V.mat(1:3,1:3).^2));
M      = V.mat;
shift0 = inv(spm_matrix(v/2));
shift  = spm_matrix(V.dim(1:3)/2);
zoom   = spm_matrix([0 0 0 0 0 0 sqrt(sum(V.mat(1:3,1:3).^2))]);

MT0 = zoom;

MT1 = MT0 * shift * spm_matrix([0 0 0 0 pi 0]) * inv(shift);

MS0 = spm_matrix(v([3 2 1])/2) * spm_matrix([0 0 0 0 pi/2]) ...
	* shift0 * zoom;

MS1 = MS0 * shift * spm_matrix([0 0 0 0 0 pi]) * inv(shift);

MC0 = spm_matrix(v([3 1 2])/2) * spm_matrix([0 0 0 0 0 pi/2]) ...
	* spm_matrix([0 0 0 pi/2]) * shift0 * zoom;

MC1 = MC0 * shift * spm_matrix([0 0 0 0 0 pi]) * inv(shift);

if 1,  % 0.5mm resolution
	zm  = diag([2 2 1 1]); % vox vox mm 1
	v   = v*2;
	MT0 = zm*MT0; MT1 = zm*MT1;
	MS0 = zm*MS0; MS1 = zm*MS1;
	MC0 = zm*MC0; MC1 = zm*MC1;
end;

if det(M)<0,
	flp=spm_matrix([0 (v(2)+1) 0  0 0 0 1 -1 1]); MT0=flp*MT0; MT1=flp*MT1;
	flp=spm_matrix([0 (v(2)+1) 0  0 0 0 1 -1 1]); MS0=flp*MS0; MS1=flp*MS1;
	flp=spm_matrix([0 (v(1)+1) 0  0 0 0 1 -1 1]); MC0=flp*MC0; MC1=flp*MC1;
end;

linfun('Rendering: Transverse 1..');	rend{1} = make_struct(V,MT0,v([1 2]));
linfun('Rendering: Transverse 2..');	rend{2} = make_struct(V,MT1,v([1 2]));
linfun('Rendering: Saggital 1..');	rend{3} = make_struct(V,MS0,v([3 2]));
linfun('Rendering: Saggital 2..');	rend{4} = make_struct(V,MS1,v([3 2]));
linfun('Rendering: Coronal 1..');	rend{5} = make_struct(V,MC0,v([3 1]));
linfun('Rendering: Coronal 2..');	rend{6} = make_struct(V,MC1,v([3 1]));

linfun('Rendering: Save..');
save(oname,'rend');
linfun(' ');
disp_renderings(rend);
spm_print;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function str = make_struct(V,M,D)
	[ren,dep]=make_pic(V,M,D);
	str=struct('M',M/V.mat,'ren',ren,'dep',dep);
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
ren       = reshape([dx(:) dy(:) dz(:).*brn(:)]*(vec*M(1:3,1:3))',[size(dx) 1]);
tmp       = find(ren<0);
ren(tmp)  = 0;
ren(msk)  = ren(msk)-0.2;
ren       = ren+0.2;
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
