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
% so it may not work every time.  Note that the output images are
% slightly smooth.  This is intentional - so that the images can be
% used for rendering.  A simple modification of the code could would
% easily produce crisper results.
%
% Inputs:
% xxx_seg1.img & xxx_seg2.img - grey and white matter segments extracted
% using spm_segment.
%
% Outputs:
% The extracted brain is written to "brain_xxx.img" in the same
% directory as xxx_seg1.img.  A "render_xxx.mat" file is also produced
% that can be used for rendering activations on to.
%
%_______________________________________________________________________
% %W% John Ashburner %E%

linfun = inline('fprintf([''%-60s%s''],x,[sprintf(''\b'')*ones(1,60)])','x');

SPMid = spm('FnBanner',mfilename,'%W%');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','XBrain');
spm_help('!ContextHelp','spm_xbrain.m');
P=spm_get(2,'*_seg?.img','Select gray and white matter images');
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
th = 0.6;	% Threshold for the erosions
spm_progress_bar('Init',25,'Extracting Brain','Iterations completed');
for j=1:32,
	if j>2, th=0.15; end;	% Dilate after two iterations of erosion. 
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

% Write the extracted brain
%-----------------------------------------------------------------------
linfun('Writing volume');
VO=struct('fname',fname,...
	  'dim',[VW.dim(1:3) spm_type('uint8')],...
	  'mat',VW.mat,...
	  'pinfo',[1/255 0 0]',...
	  'descrip','extracted brain');
spm_write_vol(VO,br);
linfun(' ');

% Produce rendering
%-----------------------------------------------------------------------
matname = fullfile(pth,['render_' nam '.mat']);
spm('FigName','Xbrain: render',Finter,CmdLine);
renviews(fname,matname);

spm('FigName','Xbrain: done',Finter,CmdLine);
spm('Pointer')

return;
%_______________________________________________________________________

%_______________________________________________________________________
function renviews(P0,oname)
% Produce images for rendering activations to
%
% FORMAT renviews(P0,oname)
% P0     - filename of image to render.
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

V      = spm_vol(P0);
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

linfun('Rendering: Transverse 1..');	rend{1} = make_struct(V,MT0,v([1 2]));
linfun('Rendering: Transverse 2..');	rend{2} = make_struct(V,MT1,v([1 2]));
linfun('Rendering: Saggital 1..');	rend{3} = make_struct(V,MS0,v([3 2]));
linfun('Rendering: Saggital 2..');	rend{4} = make_struct(V,MS1,v([3 2]));
linfun('Rendering: Coronal 1..');	rend{5} = make_struct(V,MC0,v([3 1]));
linfun('Rendering: Coronal 2..');	rend{6} = make_struct(V,MC1,v([3 1]));

linfun('Rendering: Save..');
save(oname,'rend');
linfun(' ');
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
[REN, zbuf, X, Y, Z] = spm_render_vol(V, M, D, [0.5 1]);
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
[t,dx,dy,dz] = spm_sample_vol(V,sX,sY,sZ,3);
ren       = reshape([dx(:) dy(:) dz(:).*brn(:)]*(vec*M(1:3,1:3))',[size(dx) 1]);
tmp       = find(ren<0);
ren(tmp)  = 0;
ren(msk)  = ren(msk)-0.2;
ren       = ren+0.2;
mx        = max(ren(:));
ren       = ren/mx;
return;
%_______________________________________________________________________
