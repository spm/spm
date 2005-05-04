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
% The outline and grid are superimposed at intensity defaults.grid,
% defaulting to 0.6.
%
% A default colormap of 64 levels is assumed. The pointlist image is
% scaled to fit in the interval [0.25,1]*64 for display. Flat images
% are scaled to 1*64.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston et al.
% $Id: spm_mip.m 112 2005-05-04 18:20:52Z john $



%-Get GRID value
%-----------------------------------------------------------------------
global defaults
try
	GRID  = defaults.grid;
catch
	GRID  = 0.6;
end
try
	units = defaults.units;
catch
	units = {'mm' 'mm' 'mm'};
end

%-Scale & offset point list values to fit in [0.25,1]
%-----------------------------------------------------------------------
Z    = Z(:)';
Z    = Z - min(Z);
m    = max(Z);
if m
	Z = (1 + 3*Z/m)/4;
else
	Z = ones(1,length(Z));
end

%-Single slice case
%=======================================================================
if DIM(3) == 1,
        set(gca,'Position',[0.08 0.62 0.48 0.32])
	IJK = round(M\[XYZ; ones(1,size(XYZ,2))]);
	MIP = full(sparse(IJK(1,:),IJK(2,:),Z,DIM(1),DIM(2)));
	MIP = (1 - MIP'/max(MIP(:)))*64;
        x   = M(1:2,1:2)*[1 DIM(1); 1 DIM(2)];
	image([x(1,1) x(1,2)],[x(2,1) x(2,2)],MIP);
        axis xy
        xlabel(units{1})
	axis tight
        set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])
	return
end

%-3d case
%=======================================================================
%-Load mip and create maximum intensity projection
%-----------------------------------------------------------------------
load MIP
mip  = mip96*GRID;
c    = [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1 
	1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1] - 0.5;
c    = (M(1:3,1:3)*c')';
dim  = [(max(c) - min(c)) size(mip)];
d    = spm_project(Z,round(XYZ),dim);
mip  = max(d,mip);
image(rot90((1 - mip)*64)); axis tight; axis off;
