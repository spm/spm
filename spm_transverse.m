function spm_transverse(SPM,VOL,hReg,junk1)
% Rendering of regional effects [SPM{T/F}] on transverse sections
% FORMAT spm_transverse(SPM,VOL,hReg)
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .Z     - minimum of n Statistics {filtered on u and k}
% .STAT  - distribution {Z, T, X or F}     
% .u     - height threshold
% .XYZ   - location of voxels {voxel coords}
%
% VOL    - structure containing details of volume analysed
%        - required fields are:
% .iM    - mm  -> voxels matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels}
%
% hReg   - handle of MIP XYZ registry object (see spm_XYZreg for details)
%
% spm_transverse automatically updates its co-ordinates from the
% registry, but clicking on the slices has no effect on the registry.
% i.e., the updating is one way only.
%
% See also: spm_getSPM
%_______________________________________________________________________
%
% spm_transverse is called by the SPM results section and uses
% variables in SPM and VOL to create three transverse sections though a
% background image.  Regional foci from the selected SPM{T/F} are
% rendered on this image.
%
% Although the SPM{.} adopts the neurological convention (left = left)
% the rendered images follow the same convention as the original data.
%_______________________________________________________________________
% %W% Karl Friston - modified by John Ashburner %E%

if nargin == 3 & isstruct(SPM),
	init(SPM,VOL,hReg);
end;
if ~isstruct(SPM) & strcmp(lower(SPM),'setcoords'),
	reposition(VOL);
end;
if ~isstruct(SPM) & strcmp(lower(SPM),'clear'),
	clear_global;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function init(SPM,VOL,hReg)

%-Get figure handles
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');


%-Get the image on which to render
%-----------------------------------------------------------------------
spms   = spm_get(1,'.img','select an image for rendering');
spm('Pointer','Watch');

%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('Clear',Fgraph);

global transv
transv = struct('blob',[],'V',spm_vol(spms),'h',[],'hReg',hReg,'fig',Fgraph);
transv.blob = struct('xyz', round(SPM.XYZ), 't',SPM.Z, 'dim',VOL.DIM(1:3),...
	             'iM',VOL.iM,...
		     'vox', sqrt(sum(VOL.M(1:3,1:3).^2)), 'u', SPM.u);

%-Get current location and convert to pixel co-ordinates
%-----------------------------------------------------------------------
xyzmm  = spm_XYZreg('GetCoords',transv.hReg);
xyz    = round(transv.blob.iM(1:3,:)*[xyzmm; 1]);

% extract data from SPM [at one plane separation]
% and get background slices
%----------------------------------------------------------------------
dim    = ceil(transv.blob.dim(1:3)'.*transv.blob.vox);
A      = transv.blob.iM*transv.V.mat;
hld    = 0;

zoomM  = inv(spm_matrix([0 0 -1  0 0 0  transv.blob.vox([1 2]) 1]));
zoomM1 =     spm_matrix([0 0  0  0 0 0  transv.blob.vox([1 2]) 1]);

Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)) < 0.5);
T2     = full(sparse(transv.blob.xyz(1,Q),transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
T2     = spm_slice_vol(T2,zoomM,dim([1 2]),hld);
D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3);0 0 0 1]*A;
D2     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
maxD   = max(D2(:));

if transv.blob.dim(3) > 1

	Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)+1) < 0.5);
	T1     = full(sparse(transv.blob.xyz(1,Q),...
			transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
	T1     = spm_slice_vol(T1,zoomM,dim([1 2]),hld);
	D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)+1;0 0 0 1]*A;
	D1     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
	maxD   = max([maxD ; D1(:)]);

	Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)-1) < 0.5);
	T3     = full(sparse(transv.blob.xyz(1,Q),...
			transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
	T3     = spm_slice_vol(T3,zoomM,dim([1 2]),hld);
	D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)-1;0 0 0 1]*A;
	D3     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
	maxD   = max([maxD ; D3(:)]);
end

d      = max(T2(:));
D2     = D2/maxD;
if transv.blob.dim(3) > 1,
	D1 = D1/maxD;
	D3 = D3/maxD;
	d = max([d ; T1(:) ; T3(:) ; eps]);
end;

%-Configure {128 level} colormap
%-----------------------------------------------------------------------
cmap   = get(Fgraph,'Colormap');
if size(cmap,1) ~= 128
	figure(Fgraph)
	spm_figure('Colormap','gray-hot')
	cmap   = get(Fgraph,'Colormap');
end

D      = length(cmap)/2;
Q      = T2(:) > transv.blob.u; T2 = T2(Q)/d; D2(Q) = 1 + T2; T2 = D*D2;

if transv.blob.dim(3) > 1
    Q  = T1(:) > transv.blob.u; T1 = T1(Q)/d; D1(Q) = 1 + T1; T1 = D*D1;
    Q  = T3(:) > transv.blob.u; T3 = T3(Q)/d; D3(Q) = 1 + T3; T3 = D*D3;
end

set(Fgraph,'Units','pixels')
siz    = get(Fgraph,'Position');
siz    = siz(3:4);

P = xyz.*transv.blob.vox';

%-Render activation foci on background images
%-----------------------------------------------------------------------
if transv.blob.dim(3) > 1
	zm     = min([(siz(1) - 120)/(dim(1)*3),(siz(2)/2 - 60)/dim(2)]);
	xo     = (siz(1)-(dim(1)*zm*3)-120)/2;
	yo     = (siz(2)/2 - dim(2)*zm - 60)/2;

	transv.h(1) = axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(1)*zm dim(2)*zm],'DeleteFcn','spm_transverse(''clear'');');
	transv.h(2) = image(rot90(spm_grid(T1)));
	axis image; axis off;
	title(sprintf('z = %0.0fmm',(xyzmm(3) - transv.blob.vox(3))));
	transv.h(3) = line([1 1]*P(1),[0 dim(2)],'Color','w');
	transv.h(4) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w');

	transv.h(5) = axes('Units','pixels','Parent',Fgraph,'Position',[40+dim(1)*zm+xo 20+yo dim(1)*zm dim(2)*zm],'DeleteFcn','spm_transverse(''clear'');');
	transv.h(6) = image(rot90(spm_grid(T2)));
	axis image; axis off;
	title(sprintf('z = %0.0fmm',xyzmm(3)));
	transv.h(7) = line([1 1]*P(1),[0 dim(2)],'Color','w');
	transv.h(8) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w');

	transv.h(9) = axes('Units','pixels','Parent',Fgraph,'Position',[60+dim(1)*zm*2+xo 20+yo dim(1)*zm dim(2)*zm],'DeleteFcn','spm_transverse(''clear'');');
	transv.h(10) = image(rot90(spm_grid(T3)));
	axis image; axis off;
	title(sprintf('z = %0.0fmm',(xyzmm(3) + transv.blob.vox(3))));
	transv.h(11) = line([1 1]*P(1),[0 dim(2)],'Color','w');
	transv.h(12) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w');

	% colorbar
	%-----------------------------------------------------------------------
	q      = [80+dim(1)*zm*3+xo 20+yo 20 dim(2)*zm];
	transv.h(13) = axes('Units','pixels','Parent',Fgraph,'Position',q,'Visible','off','DeleteFcn','spm_transverse(''clear'');');
	transv.h(14) = image([0 d/32],[0 d],[1:D]' + D);
	str    = [SPM.STAT ' value'];
	axis xy; title(str,'FontSize',9);
	set(gca,'XTickLabel',[]);

else,
	zm     = min([(siz(1) - 80)/dim(1),(siz(2)/2 - 60)/dim(2)]);
	xo     = (siz(1)-(dim(1)*zm)-80)/2;
	yo     = (siz(2)/2 - dim(2)*zm - 60)/2;

	transv.h(1) = axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(1)*zm dim(2)*zm],'DeleteFcn','spm_transverse(''clear'');');
	transv.h(2) = image(rot90(spm_grid(T2)));
	axis image; axis off;
	title(sprintf('z = %0.0fmm',xyzmm(3)));
	transv.h(3) = line([1 1]*P(1),[0 dim(2)],'Color','w');
	transv.h(4) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w');

	% colorbar
	%-----------------------------------------------------------------------
	q      = [40+dim(1)*zm+xo 20+yo 20 dim(2)*zm];
	transv.h(5) = axes('Units','pixels','Parent',Fgraph,'Position',q,'Visible','off','DeleteFcn','spm_transverse(''clear'');');
	transv.h(6) = image([0 d/32],[0 d],[1:D]' + D);
	str    = [SPM.STAT ' value'];
	axis xy; title(str,'FontSize',9);
	set(gca,'XTickLabel',[]);
end;

spm_XYZreg('Add2Reg',transv.hReg,transv.h(1), 'spm_transverse');

for h=transv.h,
	set(h,'DeleteFcn','spm_transverse(''clear'');');
end;

%-Reset pointer
%-----------------------------------------------------------------------
spm('Pointer','Arrow')
return;
%_______________________________________________________________________
%_______________________________________________________________________
function reposition(xyzmm)
global transv
if ~isstruct(transv), return; end;

spm('Pointer','Watch');


%-Get current location and convert to pixel co-ordinates
%-----------------------------------------------------------------------
% xyzmm  = spm_XYZreg('GetCoords',transv.hReg)
xyz    = round(transv.blob.iM(1:3,:)*[xyzmm; 1]);

% extract data from SPM [at one plane separation]
% and get background slices
%----------------------------------------------------------------------
dim    = ceil(transv.blob.dim(1:3)'.*transv.blob.vox);
A      = transv.blob.iM*transv.V.mat;
hld    = 0;

zoomM  = inv(spm_matrix([0 0 -1  0 0 0  transv.blob.vox([1 2]) 1]));
zoomM1 =     spm_matrix([0 0  0  0 0 0  transv.blob.vox([1 2]) 1]);

Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)) < 0.5);
T2     = full(sparse(transv.blob.xyz(1,Q),transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
T2     = spm_slice_vol(T2,zoomM,dim([1 2]),hld);
D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3);0 0 0 1]*A;
D2     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
maxD   = max(D2(:));

if transv.blob.dim(3) > 1

	Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)+1) < 0.5);
	T1     = full(sparse(transv.blob.xyz(1,Q),...
			transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
	T1     = spm_slice_vol(T1,zoomM,dim([1 2]),hld);
	D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)+1;0 0 0 1]*A;
	D1     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
	maxD   = max([maxD ; D1(:)]);

	Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)-1) < 0.5);
	T3     = full(sparse(transv.blob.xyz(1,Q),...
			transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
	T3     = spm_slice_vol(T3,zoomM,dim([1 2]),hld);
	D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)-1;0 0 0 1]*A;
	D3     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
	maxD   = max([maxD ; D3(:)]);
end

d      = max(T2(:));
D2     = D2/maxD;
if transv.blob.dim(3) > 1,
	D1 = D1/maxD;
	D3 = D3/maxD;
	d = max([d ; T1(:) ; T3(:) ; eps]);
end;

%-Configure {128 level} colormap
%-----------------------------------------------------------------------
cmap   = get(transv.fig,'Colormap');
if size(cmap,1) ~= 128
	figure(transv.fig)
	spm_figure('Colormap','gray-hot')
	cmap   = get(transv.fig,'Colormap');
end

D      = length(cmap)/2;
Q      = T2(:) > transv.blob.u; T2 = T2(Q)/d; D2(Q) = 1 + T2; T2 = D*D2;

if transv.blob.dim(3) > 1
    Q  = T1(:) > transv.blob.u; T1 = T1(Q)/d; D1(Q) = 1 + T1; T1 = D*D1;
    Q  = T3(:) > transv.blob.u; T3 = T3(Q)/d; D3(Q) = 1 + T3; T3 = D*D3;
end

P = xyz.*transv.blob.vox';

%-Render activation foci on background images
%-----------------------------------------------------------------------
if transv.blob.dim(3) > 1

	set(transv.h(2),'Cdata',rot90(spm_grid(T1)));
	set(get(transv.h(1),'Title'),'String',sprintf('z = %0.0fmm',(xyzmm(3) - transv.blob.vox(3))));
	set(transv.h(3),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
	set(transv.h(4),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));

	set(transv.h(6),'Cdata',rot90(spm_grid(T2)));
	set(get(transv.h(5),'Title'),'String',sprintf('z = %0.0fmm',(xyzmm(3))));
	set(transv.h(7),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
	set(transv.h(8),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));

	set(transv.h(10),'Cdata',rot90(spm_grid(T3)));
	set(get(transv.h(9),'Title'),'String',sprintf('z = %0.0fmm',(xyzmm(3) + transv.blob.vox(3))));
	set(transv.h(11),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
	set(transv.h(12),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));

	% colorbar
	%-----------------------------------------------------------------------
	set(transv.h(14), 'Ydata',[0 d], 'Cdata',[1:D]' + D);
	set(transv.h(13),'XTickLabel',[],'Ylim',[0 d]);
else,
	set(transv.h(2),'Cdata',rot90(spm_grid(T2)));
	set(get(transv.h(1),'Title'),'String',sprintf('z = %0.0fmm',xyzmm(3)));
	set(transv.h(3),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
	set(transv.h(4),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));

	% colorbar
	%-----------------------------------------------------------------------
	set(transv.h(6), 'Ydata',[0 d], 'Cdata',[1:D]' + D);
	set(transv.h(5),'XTickLabel',[],'Ylim',[0 d]);
end;


%-Reset pointer
%-----------------------------------------------------------------------
spm('Pointer','Arrow')
return;
%_______________________________________________________________________
%_______________________________________________________________________
function clear_global
global transv
if isstruct(transv),
	for h=transv.h,
		if ishandle(h), set(h,'DeleteFcn',';'); end;
	end;
	for h=transv.h,
		if ishandle(h), delete(h); end;
	end;
	transv = [];
	clear global transv;
end;
return;
