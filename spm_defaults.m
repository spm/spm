
% Sets the defaults which are used by SPM
%
% FORMAT spm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
% Individual users can make copies which can be stored in their own
% matlab subdirectories. If ~/matlab is ahead of the SPM directory
% in the MATLABPATH, then the users own personal defaults are used.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% %W%	%E%

global CWD PRINTSTR LOGFILE CMDLINE GRID proj_MultiPage
global PET_UFp PET_DIM PET_VOX PET_TYPE PET_SCALE PET_OFFSET PET_ORIGIN PET_DESCRIP
global fMRI_UFp fMRI_DIM fMRI_VOX fMRI_TYPE fMRI_SCALE fMRI_OFFSET fMRI_ORIGIN fMRI_DESCRIP

% Default command for printing
%-----------------------------------------------------------------------
PRINTSTR = 'print -dpsc2 -append spm97.ps';

% Log user input to SPM. If LOGFILE is '', then don't log.
%-----------------------------------------------------------------------
LOGFILE = '';

% Command Line
% Values can be:
%	0 - button clicking only
%	1 - input filenames etc. from the command line
%-----------------------------------------------------------------------
CMDLINE = 0;

% GRID should be in the range of 0 to 1.
% It determines the intensity of any grids which are superimposed
% on displayed images.
%-----------------------------------------------------------------------
GRID = 0.6;

% proj_MultiPage controls whether multi-page tables of cluster
% statistics are printed out by spm_projections when called from
% spm_projections_ui (the SPM{Z} button)
% This paging facility enables a more exhaustive listing of the 
% suprathreshold clusters in an SPM
%-----------------------------------------------------------------------
proj_MultiPage = 0;	% No multipage tables
%proj_MultiPage = 1;	% Enable paging


% Header defaults
%=======================================================================

% PET header defaults
%-----------------------------------------------------------------------
PET_DIM      = [128 128 63];		% Dimensions [x y z]
PET_VOX      = [2.09 2.09 2.37];		% Voxel size [x y z]
PET_TYPE     = 2;			% Data type
PET_SCALE    = 1.0;			% Scaling coeficient
PET_OFFSET   = 0;			% Offset in bytes
PET_ORIGIN   = [0 0 0];			% Origin in voxels
PET_DESCRIP  = 'SPM-compatible';

% fMRI header defaults
%-----------------------------------------------------------------------
fMRI_DIM     = [64 64 64];		% Dimensions [x y z]
fMRI_VOX     = [3 3 3];			% Voxel size [x y z]
fMRI_TYPE    = 4;			% Data type
fMRI_SCALE   = 1;			% Scaling coeficient
fMRI_OFFSET  = 0;			% Offset in bytes
fMRI_ORIGIN  = [0 0 0];			% Origin in voxels
fMRI_DESCRIP = 'SPM-compatible';

% Stats defaults
%=======================================================================

% UFp - Upper tail F probability threshold used to filter voxels after
% stats
%-----------------------------------------------------------------------
PET_UFp  = 0.05;
fMRI_UFp = 0.001;


% Realignment defaults
%=======================================================================
global sptl_WhchPtn sptl_CrtWht sptl_DjstFMRI sptl_MskOptn

% Which Option?
% This gives the flexibility for coregistering and reslicing images
% seperately, or alternatively reduces flexibility and keeps things
% simple for the user.
%-----------------------------------------------------------------------
%sptl_WhchPtn = 1;	% Combine coregistration and reslicing 
sptl_WhchPtn = -1;	% Allow separate coregistration and reslicing


% Create What? Give flexibility about which images are written resliced,
% or remove an extra question.
%-----------------------------------------------------------------------
%sptl_CrtWht = 1;	% All Images + Mean Image
sptl_CrtWht = -1;	% Full options

% Adjust FMRI?
% Adjust the data (fMRI) to remove movement-related components
% The adjustment procedure is based on a autoregression-moving 
% average-like model of the effect of position on signal and 
% explicitly includes a spin excitation history effect.
%-----------------------------------------------------------------------
%sptl_DjstFMRI =  0;	% Never adjust
%sptl_DjstFMRI =  1;	% Always adjust
sptl_DjstFMRI = -1;	% Optional adjust

% Mask Option.
% To avoid artifactual movement-related variance the realigned
% set of images can be internally masked, within the set (i.e.
% if any image has a zero value at a voxel than all images have
% zero values at that voxel).  Zero values occur when regions
% 'outside' the image are moved 'inside' the image during
% realignment.
%-----------------------------------------------------------------------
%sptl_MskOptn = -1;	% Optional mask
sptl_MskOptn =  1;	% Always mask


% Coregistration defaults
%=======================================================================
global QUICK_COREG

% Hidden option to just do a quick between mode coregistration
% This option misses out the segmenting, and coregistering segments
% steps.
%-----------------------------------------------------------------------
%QUICK_COREG = 1;	% Quick and simple
QUICK_COREG = 0;	% Full


% Spatial Normalisation defaults
%=======================================================================
global sptl_Ornt sptl_CO sptl_NAP sptl_NBss sptl_NItr sptl_BB sptl_Vx sptl_Rglrztn

% Orientation/position of images. Used as a starting estimate for
% affine normalisation.
%-----------------------------------------------------------------------
%sptl_Ornt = [0 0 0  0 0 0  1 1 1 0 0 0]; % Neurological Convention (R is R)
sptl_Ornt = [0 0 0  0 0 0 -1 1 1 0 0 0]; % Radiological Convention (L is R)


% Customisation Option. Include option to customise the normalisation
% options.
% ie. # affine params, # nonlinear basis images & # nonlinear iterations.
%-----------------------------------------------------------------------
%sptl_CO = -1;		% Allow customised
sptl_CO  =  1;		% Disallow Customised

% Number of nonlinear basis functions
%-----------------------------------------------------------------------
%sptl_NBss = [0 0 0];	% None (ie. perform affine normalisation only).
sptl_NBss = [7 8 7];

% Regularization fudge factor:
%	small values	-> less regularization -> more warping
%	large values	-> more regularization -> less warping
%-----------------------------------------------------------------------
sptl_Rglrztn = 0.01;

% Number of iterations of nonlinear spatial normalisation.
%-----------------------------------------------------------------------
sptl_NItr = 12;

% Estimate the spatial normalization parameters from the brain only
% or use the whole head.
%-----------------------------------------------------------------------
%sptl_MskBrn = 0; % Estimate from the whole head
sptl_MskBrn = 1; % Estimate from the brain only

% Bounding Box. The definition of the volume of the normalised image
% which is written (mm relative to AC).
% [[lowX lowY lowZ];[highX highY highZ]]
%-----------------------------------------------------------------------
sptl_BB = [[-78 -112 -50];[78 76 85]];
%sptl_BB = [[-90 -126 -72];[91 91 109]];

% Voxel sizes in mm of the normalised images
%-----------------------------------------------------------------------
sptl_Vx = [2 2 2];	% 2mm x 2mm x 2mm
%sptl_Vx = [1.5 1.5 1.5];
