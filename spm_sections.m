function spm_sections(SPM,VOL,hReg)
% rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections(SPM,VOL,hReg)
%
% SPM  - SPM structure      {'Z' 'n' 'STAT' 'df' 'u' 'k'}
% VOL  - Spatial structure  {'R' 'FWHM' 'S' 'DIM' 'VOX' 'ORG' 'M' 'XYZ' 'QQ'}
% hReg - handle of MIP register
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_sections is called by spm_results and uses variables in SPM and
% VOL to create three orthogonal sections though a background image.
% Regional foci from the selected SPM are rendered on this image.
%
% Although the SPM{Z} adopts the neurological convention (left = left)
% the rendered images follow the same convention as the original data.
%
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
L      = spm_XYZreg('GetCoords',hReg);
L      = round(IM(1:3,:)*[L ; 1]);


% extract data from SPM{t} [sagittal {Ts} coronal {Tc} transverse {Tt}]
% and get background slices
%----------------------------------------------------------------------
vox    = sqrt(sum(VOL.M(1:3,1:3).^2));
dim    = ceil(VOL.DIM(1:3)'.*vox);
Vs     = spm_vol(spms);
A      = VOL.M\Vs.mat;
hld    = 1;

Q      = find(abs(XYZ(1,:) - L(1)) < 0.5);
Ts     = full(sparse(XYZ(3,Q),XYZ(2,Q),SPM.Z(Q),VOL.DIM(3),VOL.DIM(2)));
zoomM  = spm_matrix([0 0 -1  0 0 0  vox([3 2]) 1]);
Ts     = spm_slice_vol(Ts,inv(zoomM),dim([3 2]),hld);
D      = zoomM*[0 0 1 0;0 1 0 0;1 0 0 -L(1);0 0 0 1]*A;
Ds     = spm_slice_vol(Vs,inv(D),dim([3 2]),1);

Q      = find(abs(XYZ(2,:) - L(2)) < 0.5);
Tc     = full(sparse(XYZ(3,Q),XYZ(1,Q),SPM.Z(Q),VOL.DIM(3),VOL.DIM(1)));
zoomM  = spm_matrix([0 0 -1  0 0 0  vox([3 1]) 1]);
Tc     = spm_slice_vol(Tc,inv(zoomM),dim([3 1]),hld);
D      = zoomM*[0 0 1 0;1 0 0 0;0 1 0 -L(2);0 0 0 1]*A;
Dc     = spm_slice_vol(Vs,inv(D),dim([3 1]),1);

Q      = find(abs(XYZ(3,:) - L(3)) < 0.5);
Tt     = full(sparse(XYZ(1,Q),XYZ(2,Q),SPM.Z(Q),VOL.DIM(1),VOL.DIM(2)));
zoomM  = spm_matrix([0 0 -1  0 0 0  vox([1 2]) 1]);
Tt     = spm_slice_vol(Tt,inv(zoomM),dim([1 2]),hld);
D      = zoomM*[1 0 0 0;0 1 0 0;0 0 1 -L(3);0 0 0 1]*A;
Dt     = spm_slice_vol(Vs,inv(D),dim([1 2]),1);

P = L.*vox';

colormap([gray(64); pink(64)])

scal   = max([Ds(:) ; Dc(:) ; Dt(:)]);
Dt     = Dt/scal;
Ds     = Ds/scal;
Dc     = Dc/scal;
d      = max([Ts(:); Tc(:); Tt(:)]);
D      = length(colormap)/2;
Q      = Ts(:) > SPM.u; Ts = Ts(Q)/d; Ds(Q) = 1 + Ts; Ts = D*Ds;
Q      = Tc(:) > SPM.u; Tc = Tc(Q)/d; Dc(Q) = 1 + Tc; Tc = D*Dc;
Q      = Tt(:) > SPM.u; Tt = Tt(Q)/d; Dt(Q) = 1 + Tt; Tt = D*Dt;


set(Fgraph,'Units','pixels')
siz    = get(Fgraph,'Position');
siz    = siz(3:4);
zm     = min([(siz(1) - 60)/(dim(2)+dim(1)),(siz(2)/2 - 100)/(dim(1)+dim(3))]);
xo     = (siz(1)-(dim(2)+dim(1))*zm-60)/2;
yo     = (siz(2)/2 - (dim(1)+dim(3))*zm - 100)/2;

% render activation foci on background image
%-----------------------------------------------------------------------
axes('Units','pixels','Parent',Fgraph,'Position',[20+xo dim(1)*zm+60+yo dim(2)*zm dim(3)*zm]);
image(Ts)
axis image; axis('xy'); axis off;
title 'sagittal'
line([0 dim(2)],[1 1]*P(3),'Color','w')
line([1 1]*P(2),[0 dim(3)],'Color','w')

axes('Units','pixels','Parent',Fgraph,'Position',[dim(2)*zm+40+xo dim(1)*zm+60+yo dim(1)*zm dim(3)*zm]);
image(Tc)
axis image; axis('xy'); axis off;
title 'coronal'
line([0 dim(1)],[1 1]*P(3),'Color','w')
line([1 1]*P(1),[0 dim(3)],'Color','w')

axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(2)*zm dim(1)*zm])
image(Tt)
axis image; axis off;
title 'transverse'
line([0 dim(2)],[1 1]*P(1),'Color','w')
line([1 1]*P(2),[0 dim(1)],'Color','w')

% colorbar
%-----------------------------------------------------------------------
q      = [dim(2)*zm+40+xo 20+yo 20 dim(1)*zm];
axes('Units','pixels','Parent',Fgraph,'Position',q,'Visible','off')
image([0 d/32],[0 d],[1:D]' + D)
str    = [SPM.STAT ' value'];

axis xy; title(str,'FontSize',9);
set(gca,'XTickLabel',[])

% reset pointer (and x locations is necessary)
%----------------------------------------------------------------------
set([Finter,Fgraph],'Pointer','Arrow')
