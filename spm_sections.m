function spm_sections(SPM,VOL,hReg)
% Rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections(SPM,VOL,hReg)
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
% .M     - voxels - > mm matrix
% .iM    - mm -> voxels matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels}
%
% hReg   - handle of MIP XYZ registry object (see spm_XYZreg for details)
%
% See also: spm_getSPM
%_______________________________________________________________________
%
% spm_sections is called by the SPM results section and uses variables
% in SPM and VOL to create three orthogonal sections though a
% background image.  Regional foci from the selected SPM are rendered
% on this image.
%
% Although the SPM{Z} adopts the neurological convention (left = left)
% the rendered images follow the same convention as the original data.
%
%_______________________________________________________________________
% %W% Karl Friston %E%


%-Get figure handles
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');


% get the image on which to render
%-----------------------------------------------------------------------
spms   = spm_get(1,'.img','select an image for rendering');
spm('Pointer','Watch');


%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('ClearPane',Fgraph,'RNP');


%-Get current location & convert to pixel co-ordinates
%-----------------------------------------------------------------------
xyz    = round(VOL.iM(1:3,:)*[spm_XYZreg('GetCoords',hReg); 1]);


%-Extract data from SPM [sagittal {Zs} coronal {Zc} transverse {Zt}]
% and get background slices
%----------------------------------------------------------------------
vox    = VOL.VOX';
dim    = ceil(VOL.DIM.*VOL.VOX)';
Vs     = spm_vol(spms);
A      = VOL.M\Vs.mat;
hld    = 1;

Q      = find(abs(SPM.XYZ(1,:) - xyz(1)) < 0.5);
Zs     = full(sparse(SPM.XYZ(3,Q),SPM.XYZ(2,Q),SPM.Z(Q),VOL.DIM(3),VOL.DIM(2)));
zoomM  = spm_matrix([0 0 -1  0 0 0  vox([3 2]) 1]);
Zs     = spm_slice_vol(Zs,inv(zoomM),dim([3 2]),hld);
D      = zoomM*[0 0 1 0;0 1 0 0;1 0 0 -xyz(1);0 0 0 1]*A;
Ds     = spm_slice_vol(Vs,inv(D),dim([3 2]),1);

Q      = find(abs(SPM.XYZ(2,:) - xyz(2)) < 0.5);
Zc     = full(sparse(SPM.XYZ(3,Q),SPM.XYZ(1,Q),SPM.Z(Q),VOL.DIM(3),VOL.DIM(1)));
zoomM  = spm_matrix([0 0 -1  0 0 0  vox([3 1]) 1]);
Zc     = spm_slice_vol(Zc,inv(zoomM),dim([3 1]),hld);
D      = zoomM*[0 0 1 0;1 0 0 0;0 1 0 -xyz(2);0 0 0 1]*A;
Dc     = spm_slice_vol(Vs,inv(D),dim([3 1]),1);

Q      = find(abs(SPM.XYZ(3,:) - xyz(3)) < 0.5);
Zt     = full(sparse(SPM.XYZ(1,Q),SPM.XYZ(2,Q),SPM.Z(Q),VOL.DIM(1),VOL.DIM(2)));
zoomM  = spm_matrix([0 0 -1  0 0 0  vox([1 2]) 1]);
Zt     = spm_slice_vol(Zt,inv(zoomM),dim([1 2]),hld);
D      = zoomM*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3);0 0 0 1]*A;
Dt     = spm_slice_vol(Vs,inv(D),dim([1 2]),1);

P = xyz.*vox';

colormap([gray(64); pink(64)])

scal   = max([Ds(:) ; Dc(:) ; Dt(:)]);
Dt     = Dt/scal;
Ds     = Ds/scal;
Dc     = Dc/scal;
d      = max([Zs(:); Zc(:); Zt(:)]);
D      = length(colormap)/2;
Q      = Zs(:) > SPM.u; Zs = Zs(Q)/d; Ds(Q) = 1 + Zs; Zs = D*Ds;
Q      = Zc(:) > SPM.u; Zc = Zc(Q)/d; Dc(Q) = 1 + Zc; Zc = D*Dc;
Q      = Zt(:) > SPM.u; Zt = Zt(Q)/d; Dt(Q) = 1 + Zt; Zt = D*Dt;


set(Fgraph,'Units','pixels')
siz    = get(Fgraph,'Position');
siz    = siz(3:4);
zm     = min([(siz(1) - 60)/(dim(2)+dim(1)),(siz(2)/2 - 100)/(dim(1)+dim(3))]);
xo     = (siz(1)-(dim(2)+dim(1))*zm-60)/2;
yo     = (siz(2)/2 - (dim(1)+dim(3))*zm - 100)/2;

% render activation foci on background image
%-----------------------------------------------------------------------
axes('Units','pixels','Parent',Fgraph,'Position',[20+xo dim(1)*zm+60+yo dim(2)*zm dim(3)*zm]);
image(Zs)
axis image; axis('xy'); axis off;
title 'sagittal'
line([0 dim(2)],[1 1]*P(3),'Color','w')
line([1 1]*P(2),[0 dim(3)],'Color','w')

axes('Units','pixels','Parent',Fgraph,'Position',[dim(2)*zm+40+xo dim(1)*zm+60+yo dim(1)*zm dim(3)*zm]);
image(Zc)
axis image; axis('xy'); axis off;
title 'coronal'
line([0 dim(1)],[1 1]*P(3),'Color','w')
line([1 1]*P(1),[0 dim(3)],'Color','w')

axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(2)*zm dim(1)*zm])
image(Zt)
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
spm('Pointer','Arrow')
