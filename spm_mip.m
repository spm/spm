function spm_mip(Z,XYZ,M,units)
% SPM maximum intensity projection
% FORMAT spm_mip(Z,XYZ,M);
% Z       - vector point list of SPM values for MIP
% XYZ     - matrix of coordinates of points (mip coordinates)
% M       - voxels - > mip matrix or size of voxels (mm)
% units   - defining space     [default {'mm' 'mm' 'mm'}]
%         - Scalar specifies intensity of grid
%_______________________________________________________________________
%
% If the data are 2 dimensional [DIM(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (Z) and their location (XYZ) in three orthogonal
% views of the brain.  It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988); unless the third
% dimension is time.
%
% This routine loads a mip outline from MIP.mat. This is an image with
% contours and grids defining the space of Talairach & Tournoux (1988).
% mip95 corresponds to the Talairach atlas, mip96 to the MNI templates.
% The outline and grid are superimposed at intensity 0.4.
%
% A customised mip outline can be used instead of the default.
%
% A default colormap of 64 levels is assumed. The pointlist image is
% scaled to fit in the interval [0.25,1]*64 for display. Flat images
% are scaled to 1*64.
%
% If M is not specified, it is assumed the XYZ locations are
% in Talairach mm.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston et al.
% $Id: spm_mip.m 4404 2011-07-22 11:36:40Z volkmar $

%-Get units and grid scaling
%--------------------------------------------------------------------------
try, units;                catch, units = {'mm' 'mm' 'mm'}; end
try, M;                    catch, M = 1; end
Grid = 0.4;

% transpose locations if necessary
%--------------------------------------------------------------------------
if size(XYZ,1) ~= 3, XYZ = XYZ';         end
if size(Z,1)   ~= 1, Z   = Z';           end
if size(M,1)   == 1, M   = speye(4,4)*M; end

%-Scale & offset point list values to fit in [0.25,1]
%==========================================================================
Z    = Z - min(Z);
mx   = max(Z);
Scal = 8;
if isempty(mx),
    Z = [];
elseif isfinite(mx) && (numel(Z) ~= 1),
    Z = (1 + Scal*Z/mx)/(Scal + 1);
else
    Z = ones(1,length(Z));
end

%-Display format
%==========================================================================
% load various grids, DXYZ, CXYZ, scale (see spm_mip_ui and spm_project)
load(spm_get_defaults('stats.results.mipmat'));

%-Single slice case
%--------------------------------------------------------------------------
if isempty(units{3}) && ~strcmp(units{2},'mm')

    %-2d case: Time-Frequency or Frequency-Frequency
    %----------------------------------------------------------------------
    mip = 4*grid_trans;
      
elseif isempty(units{3})
    
    %-2d case
    %----------------------------------------------------------------------
    mip = 4*grid_trans + mask_trans;
        
elseif strcmp(units{3},'ms') || strcmp(units{3},'Hz')
    
    %-3d case: Space-time
    %----------------------------------------------------------------------
    mip = 4*grid_time + mask_trans;

else
    %-3d case: Space
    %----------------------------------------------------------------------
    mip = 4*grid_all + mask_all;
end

% Load mip and create maximum intensity projection
%--------------------------------------------------------------------------
mip  = mip/max(mip(:));
c    = [0 0 0 ;
        0 0 1 ;
        0 1 0 ;
        0 1 1 ;
        1 0 0 ;
        1 0 1 ; 
        1 1 0 ; 
        1 1 1 ] - 0.5;
c    = c*M(1:3,1:3);
dim  = [(max(c) - min(c)) size(mip)];
d    = spm_project(Z,round(XYZ),dim,DXYZ,CXYZ);
mip  = max(d,Grid*mip);
image(rot90((1 - mip)*64)); axis tight; axis off;

