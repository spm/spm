function spm_transverse(SPM,VOL,hReg)
% Rendering of regional effects [SPM{T/F}] on transverse sections
% FORMAT spm_transverse(SPM,VOL,hReg)
%
% SPM  - SPM structure      {'Z' 'n' 'STAT' 'df' 'u' 'k'}
% VOL  - Spatial structure  {'R' 'FWHM' 'S' 'DIM' 'VOX' 'ORG' 'M' 'XYZ' 'QQ'}
% hReg - handle of MIP register
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_transverse is called by spm_results and uses variables in SPM and
% VOL  to create three transverse  sections though a background image.
% Regional foci from the selected SPM{T/F} are rendered on this image.
%
% Although the SPM{Z} adopts the neurological convention (left = left)
% the rendered images follow the same convention as the original data.
%_______________________________________________________________________
% %W% Karl Friston %E%

%-Get figure handles and filenames
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');


% get the image on which to render
%-----------------------------------------------------------------------
spms   = spm_get(1,'.img','select an image for rendering');
set([Finter,Fgraph],'Pointer','Watch');

%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('ClearPane',Fgraph,'RNP');


IM     = inv(VOL.M);
XYZ    = IM(1:3,:)*[VOL.XYZ ; ones(1,size(VOL.XYZ,2))];
L0     = spm_XYZreg('GetCoords',hReg);
L      = round(IM(1:3,:)*[L0 ; 1]);

% extract data from SPM [at one plane separation]
% and get background slices
%----------------------------------------------------------------------
vox    = sqrt(sum(VOL.M(1:3,1:3).^2));
dim    = ceil(VOL.DIM(1:3)'.*vox);
Vs     = spm_vol(spms);
A      = VOL.M\Vs.mat;
hld    = 0;

zoomM  = spm_matrix([0 0 -1  0 0 0  vox([1 2]) 1]);

Q      = find(abs(XYZ(3,:) - L(3)) < 0.5);
T2     = full(sparse(XYZ(1,Q),XYZ(2,Q),SPM.Z(Q),VOL.DIM(1),VOL.DIM(2)));
T2     = spm_slice_vol(T2,inv(zoomM),dim([1 2]),hld);
D      = zoomM*[1 0 0 0;0 1 0 0;0 0 1 -L(3);0 0 0 1]*A;
D2     = spm_slice_vol(Vs,inv(D),dim([1 2]),1);
maxD   = max(D2(:));

if VOL.DIM(3) > 1

	Q      = find(abs(XYZ(3,:) - L(3)+1) < 0.5);
	T1     = full(sparse(XYZ(1,Q),XYZ(2,Q),SPM.Z(Q),VOL.DIM(1),VOL.DIM(2)));
	T1     = spm_slice_vol(T1,inv(zoomM),dim([1 2]),hld);
	D      = zoomM*[1 0 0 0;0 1 0 0;0 0 1 -L(3)+1;0 0 0 1]*A;
	D1     = spm_slice_vol(Vs,inv(D),dim([1 2]),1);
	maxD   = max([maxD ; D1(:)]);

	Q      = find(abs(XYZ(3,:) - L(3)-1) < 0.5);
	T3     = full(sparse(XYZ(1,Q),XYZ(2,Q),SPM.Z(Q),VOL.DIM(1),VOL.DIM(2)));
	T3     = spm_slice_vol(T3,inv(zoomM),dim([1 2]),hld);
	D      = zoomM*[1 0 0 0;0 1 0 0;0 0 1 -L(3)-1;0 0 0 1]*A;
	D3     = spm_slice_vol(Vs,inv(D),dim([1 2]),1);
	maxD   = max([maxD ; D3(:)]);
end

d      = max(T2(:));
D2     = D2/maxD;
if VOL.DIM(3) > 1,
	D1 = D1/maxD;
	D3 = D3/maxD;
	d = max([d ; T1(:) ; T3(:)]);
end;

% configure {128 level} colormap
%-----------------------------------------------------------------------
colormap([gray(64); pink(64)])
D      = length(colormap)/2;
Q      = T2(:) > SPM.u; T2 = T2(Q)/d; D2(Q) = 1 + T2; T2 = D*D2;

if VOL.DIM(3) > 1
    Q  = T1(:) > SPM.u; T1 = T1(Q)/d; D1(Q) = 1 + T1; T1 = D*D1;
    Q  = T3(:) > SPM.u; T3 = T3(Q)/d; D3(Q) = 1 + T3; T3 = D*D3;
end

set(Fgraph,'Units','pixels')
siz    = get(Fgraph,'Position');
siz    = siz(3:4);

% render activation foci on background images
%-----------------------------------------------------------------------
if VOL.DIM(3) > 1
	zm     = min([(siz(1) - 120)/(dim(1)*3),(siz(2)/2 - 60)/dim(2)]);
	xo     = (siz(1)-(dim(1)*zm*3)-120)/2;
	yo     = (siz(2)/2 - dim(2)*zm - 60)/2;

	P = L.*vox';

	axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(1)*zm dim(2)*zm])
	image(rot90(spm_grid(T1)))
	axis image; axis off;
	title(sprintf('z = %0.0fmm',(L0(3) - vox(3))))
	line([1 1]*P(1),[0 dim(2)],'Color','w')
	line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w')

	axes('Units','pixels','Parent',Fgraph,'Position',[40+dim(1)*zm+xo 20+yo dim(1)*zm dim(2)*zm])
	image(rot90(spm_grid(T2)))
	axis image; axis off;
	title(sprintf('z = %0.0fmm',L0(3)))
	line([1 1]*P(1),[0 dim(2)],'Color','w')
	line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w')

	axes('Units','pixels','Parent',Fgraph,'Position',[60+dim(1)*zm*2+xo 20+yo dim(1)*zm dim(2)*zm])
	image(rot90(spm_grid(T3)))
	axis image; axis off;
	title(sprintf('z = %0.0fmm',(L0(3) + vox(3))))
	line([1 1]*P(1),[0 dim(2)],'Color','w')
	line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w')

	% colorbar
	%-----------------------------------------------------------------------
	q      = [80+dim(1)*zm*3+xo 20+yo 20 dim(2)*zm];
	axes('Units','pixels','Parent',Fgraph,'Position',q,'Visible','off')
	image([0 d/32],[0 d],[1:D]' + D)
	str    = [SPM.STAT ' value'];
	axis xy; title(str,'FontSize',9);
	set(gca,'XTickLabel',[])

else,
	zm     = min([(siz(1) - 80)/dim(1),(siz(2)/2 - 60)/dim(2)]);
	xo     = (siz(1)-(dim(1)*zm)-80)/2;
	yo     = (siz(2)/2 - dim(2)*zm - 60)/2;

	axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(1)*zm dim(2)*zm])
	image(rot90(spm_grid(T2)))
	axis image; axis off;
	title(sprintf('z = %0.0fmm',L0(3)))
	line([1 1]*P(1),[0 dim(2)],'Color','w')
	line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w')

	% colorbar
	%-----------------------------------------------------------------------
	q      = [40+dim(1)*zm+xo 20+yo 20 dim(2)*zm];
	axes('Units','pixels','Parent',Fgraph,'Position',q,'Visible','off')
	image([0 d/32],[0 d],[1:D]' + D)
	str    = [SPM.STAT ' value'];
	axis xy; title(str,'FontSize',9);
	set(gca,'XTickLabel',[])
end;

% colorbar
%-----------------------------------------------------------------------
u      = get(gca,'Position');
axes('position', [(u(1) + u(3) + 0.1) u(2) 0.01 u(3)])
image([0 d/32],[0 d],[1:D]' + D)
str    = [SPM.STAT ' value'];

axis xy; title(str,'FontSize',9);
set(gca,'XTickLabel',[])

% reset pointer
%-----------------------------------------------------------------------
set([Finter,Fgraph],'Pointer','Arrow')
