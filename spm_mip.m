function spm_mip(Z,XYZ,M,DIM)
% SPM maximum intensity projection
% FORMAT spm_mip(Z,XYZ,M,DIM);
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% M       - voxels - > mm matrix
% DIM     - image dimensions {voxels}
%_______________________________________________________________________
%
% If the data are 2 dimensional [DIM(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (Z) and their location (XYZ) in three orthogonal
% views of the brain.  It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988).
%
% This routine loads a mip putline from MIP.mat. This is an image with
% contours and grids defining the space of Talairach & Tournoux (1988).
% mip95 corresponds to the Talairach atlas, mip96 to the MNI templates.
% The outline and grid are superimposed at intensity GRID (global default),
% defaulting to 0.6.
%
% A default colormap of 64 levels is assumed. The pointlist image is
% scaled to fit in the interval [0.25,1]*64 for display. Flat images
% are scaled to 1*64.
%
%_______________________________________________________________________
% %W% Karl Friston et al. %E%


%-Get GRID value
%-----------------------------------------------------------------------
GRID = spm('GetGlobal','GRID');
if isempty(GRID), GRID = 0.6; end

%-Scale & offset point list values to fit in [0.25,1]
%-----------------------------------------------------------------------
Z    = Z(:)';
Z    = Z - min(Z);
mZ   = max(Z);
if mZ	%-Scale
	Z = (1 + 3*Z/mZ)/4;
else	%-Image is flat: set to ones:
	Z = ones(1,length(Z));
end

%-Single slice case
%=======================================================================
if DIM(3) == 1,
	VOX = sqrt(sum(M(1:3,1:3).^2));
	XYZ = round(M\[XYZ; ones(1,size(XYZ,2))]);
	mip = full(sparse(XYZ(1,:),XYZ(2,:),Z,DIM(1),DIM(2)));
	imagesc([1 DIM(1)*VOX(1)],[1 DIM(2)*VOX(2)],-mip');
	axis xy image; 
	set(gca,'FontSize',8,'TickDir','in')
	xlabel('x'), ylabel('y')
	return
end

%-3d case
%=======================================================================
%-Load mip and create maximum intensity projection
%-----------------------------------------------------------------------
load MIP
mip  = mip96*GRID;
c    = [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1 
	1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1] -0.5;
c    = (M(1:3,1:3)*c')';
dim  = [(max(c)-min(c)) size(mip)];
d    = spm_project(Z,round(XYZ),dim);
mip  = max(d,mip);
image(rot90((1 - mip)*64)); axis image; axis off;
