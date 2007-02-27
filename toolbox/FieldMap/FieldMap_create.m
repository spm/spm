function VDM=FieldMap_create(fm_imgs,epi_img,pm_defs)
%
% Function to create VDM file from fieldmap images.
% 
% This function uses routines from the FieldMap toolbox to:
% 1) Create a single field map from inout filedmap data.
% 2) Convert fieldmap to a voxel displacement map (vdm_* file).
% 3) Match vdm_* to an input EPI which should be the first image 
% that everything else will be realigned/unwarped to.
% 4) The selected EPI is unwarped and written out with the prefix 'u'.
% 
% For details about the FieldMap toolbox, see FieldMap.man. For a 
% description of the components of the structure IP, see FieldMap.m.
% For an introduction to the theoretcial and practical principles behind 
% the toolbox, see principles.man.
%
% Updated for SPM5 - 27/02/07
%_________________________________________________________________
% FieldMap_create.m                           Chloe Hutton 27/02/07
  
if nargin < 3
  error('field map images, epi image and defaults');
end

IP = FieldMap('Initialise'); % Gets default params from pm_defaults

% Define parameters for fieldmap creation
IP.et{1} = pm_defs.SHORT_ECHO_TIME;
IP.et{2} = pm_defs.LONG_ECHO_TIME;
IP.maskbrain = pm_defs.MASKBRAIN;
 
% Set parameters for unwrapping
IP.uflags.iformat = pm_defs.INPUT_DATA_FORMAT;
IP.uflags.method = pm_defs.UNWRAPPING_METHOD;
IP.uflags.fwhm = pm_defs.FWHM;
IP.uflags.pad = pm_defs.PAD;
IP.uflags.ws = pm_defs.WS;
IP.uflags.etd = pm_defs.LONG_ECHO_TIME - pm_defs.SHORT_ECHO_TIME;     

% Set parameters for brain extraction
IP.mflags.template=pm_defs.MFLAGS.TEMPLATE;
IP.mflags.fwhm=pm_defs.MFLAGS.FWHM; 
IP.mflags.nerode=pm_defs.MFLAGS.NERODE;
IP.mflags.ndilate=pm_defs.MFLAGS.NDILATE;
IP.mflags.thresh=pm_defs.MFLAGS.THRESH;

% Set parameters for unwarping 
IP.ajm = pm_defs.DO_JACOBIAN_MODULATION;
IP.blipdir = pm_defs.K_SPACE_TRAVERSAL_BLIP_DIR;
IP.tert = pm_defs.TOTAL_EPI_READOUT_TIME;
IP.epifm = pm_defs.EPI_BASED_FIELDMAPS;

% Clear any old handles etc
IP.fm = [];
IP.vdm = [];
IP.jim = [];
IP.pP = [];
IP.epiP = [];
IP.uepiP = [];
IP.vdmP = [];
ID = cell(4,1);

%----------------------------------------------------------------------
% Load measured field map data - phase and magnitude or real and imaginary
%----------------------------------------------------------------------

if ischar(fm_imgs), fm_imgs = spm_vol(fm_imgs); end
n_fms = length(fm_imgs);
switch n_fms
 case 4  % real, imaginary pairs
  for i = 1:n_fms
     IP.P{i} = spm_vol(fm_imgs(i)); 
  end
 case 2  % precalculated Hz map and magnitude image
  IP.P{1} = spm_vol(fm_imgs(1));
  IP.P{2} = spm_vol(fm_imgs(2));
 otherwise 
   error('Funny number of input fieldmap images')
end

%----------------------------------------------------------------------
% Create field map (in Hz) - this routine calls the unwrapping
%----------------------------------------------------------------------

IP.fm = FieldMap('CreateFieldMap',IP);
  
%----------------------------------------------------------------------
% Write out field map
% Outputs -> fpm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------
  
FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Smoothed phase map');

%----------------------------------------------------------------------
% Convert Hz to voxels and write voxel displacement map 
% Outputs -> vdm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

[IP.vdm, IP.vdmP]=FieldMap('FM2VDM',IP);

%----------------------------------------------------------------------
% Select an EPI to unwarp
%----------------------------------------------------------------------

if ischar(epi_img), epi_img = spm_vol(epi_img); end
IP.epiP = epi_img;

%----------------------------------------------------------------------
% Match voxel displacement map to image
% Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

if isfield(pm_defs, 'match_vdm')
  if pm_defs.match_vdm
    IP.vdmP = FieldMap('MatchVDM',IP);
  end
end

%----------------------------------------------------------------------
% Unwarp EPI
%----------------------------------------------------------------------

IP.uepiP = FieldMap('UnwarpEPI',IP);

%----------------------------------------------------------------------
% Write unwarped EPI 
% Outputs -> uNAME-OF-EPI.img
%----------------------------------------------------------------------
unwarp_info=sprintf('Unwarped EPI:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
IP.uepiP = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info);

VDM=IP.vdmP ;
return









