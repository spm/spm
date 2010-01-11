function [VDM IPcell]=FieldMap_create(fm_imgs,epi_img,pm_defs)
%
% Function to create VDM file from fieldmap images and can be called
% using FieldMap_preprocess.m
% 
% This function uses routines from the FieldMap toolbox to:
% 1) Create a single field map from input fieldmap data.
% 2) Convert fieldmap to a voxel displacement map (vdm_* file).
% 3) Match vdm_* to input EPI(s) which should be the first image 
% that each session will be realigned/unwarped to. Writes out matched vdm
% file with name extension 'session' or a user-specified name.
% 4) Each selected EPI is unwarped and written out with the prefix 'u'.
% 
% For details about the FieldMap toolbox, see FieldMap.man. For a 
% description of the components of the structure IP, see FieldMap.m.
% For an introduction to the theoretcial and practical principles behind 
% the toolbox, see principles.man.
%
% 27/02/07 - Updated for SPM5.
% 20/02/08 - Updated to work with multiple sessions
%_______________________________________________________________________
% Copyright (C) 2006 Wellcome Department of Imaging Neuroscience

% Chloe Hutton 
% $Id: FieldMap_create.m 3670 2010-01-11 14:08:02Z chloe $
%_______________________________________________________________________

  
if nargin < 3
  error('field map images, epi image and defaults');
end

IP = FieldMap('Initialise'); % Gets default params from pm_defaults

% Define parameters for fieldmap creation

if ~isfield(pm_defs,'et')
   IP.et{1} = pm_defs.SHORT_ECHO_TIME;
   IP.et{2} = pm_defs.LONG_ECHO_TIME;
else
   IP.et{1} = pm_defs.et{1};
   IP.et{2} = pm_defs.et{2};
end

if ~isfield(pm_defs,'maskbrain')
   IP.maskbrain = pm_defs.MASKBRAIN;
else
   IP.maskbrain = pm_defs.maskbrain;
end
     
% Set parameters for unwrapping
if ~isfield(pm_defs,'uflags')
   IP.uflags.iformat = pm_defs.INPUT_DATA_FORMAT;
   IP.uflags.method = pm_defs.UNWRAPPING_METHOD;
   IP.uflags.fwhm = pm_defs.FWHM;
   IP.uflags.pad = pm_defs.PAD;
   IP.uflags.ws = pm_defs.WS;
   IP.uflags.etd = pm_defs.LONG_ECHO_TIME - pm_defs.SHORT_ECHO_TIME;     
else
   IP.uflags.iformat = pm_defs.uflags.iformat;
   IP.uflags.method = pm_defs.uflags.method;
   IP.uflags.fwhm = pm_defs.uflags.fwhm;
   IP.uflags.pad = pm_defs.uflags.pad;
   IP.uflags.ws = pm_defs.uflags.ws;
   IP.uflags.etd = pm_defs.uflags.etd;      
end

% Set parameters for brain extraction
if ~isfield(pm_defs,'mflags')
   IP.mflags.template=pm_defs.MFLAGS.TEMPLATE;
   IP.mflags.fwhm=pm_defs.MFLAGS.FWHM; 
   IP.mflags.nerode=pm_defs.MFLAGS.NERODE;
   IP.mflags.ndilate=pm_defs.MFLAGS.NDILATE;
   IP.mflags.thresh=pm_defs.MFLAGS.THRESH;
   IP.mflags.reg=pm_defs.MFLAGS.REG;
else
   IP.mflags.template=pm_defs.mflags.template;
   IP.mflags.fwhm=pm_defs.mflags.fwhm; 
   IP.mflags.nerode=pm_defs.mflags.nerode;
   IP.mflags.ndilate=pm_defs.mflags.ndilate;
   IP.mflags.thresh=pm_defs.mflags.thresh;
   IP.mflags.reg=pm_defs.mflags.reg;
end

% Get FieldMap parameters 
% These may come from FieldMap_preprocess or FieldMap_Run
if ~isfield(pm_defs,'ajm')
   IP.ajm = pm_defs.DO_JACOBIAN_MODULATION;
else
   IP.ajm = pm_defs.ajm;
end
if ~isfield(pm_defs,'blipdir')
   IP.blipdir = pm_defs.K_SPACE_TRAVERSAL_BLIP_DIR;
else
   IP.blipdir = pm_defs.blipdir;
end
if ~isfield(pm_defs,'tert')
   IP.tert = pm_defs.TOTAL_EPI_READOUT_TIME;
else
   IP.tert = pm_defs.tert;
end 
if ~isfield(pm_defs,'epifm')
   IP.epifm = pm_defs.EPI_BASED_FIELDMAPS;
else
   IP.epifm = pm_defs.epifm; 
end

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

if ischar(fm_imgs)
   fm_imgs = spm_vol(fm_imgs); 
end
n_fms = length(fm_imgs);
switch n_fms
case 4  % real, imaginary pairs
   for i = 1:n_fms
      IP.P{i} = spm_vol(fm_imgs(i)); 
   end
case 2  % precalculated phase map and magnitude image
   IP.P{1} = spm_vol(fm_imgs(1));
   IP.P{2} = spm_vol(fm_imgs(2));
case 1  % precalculated and unwrapped Hz map    
   IP.pP = spm_vol(fm_imgs);
   IP.fm.fpm = spm_read_vols(IP.pP);
   IP.fm.jac = pm_diff(IP.fm.fpm,2);
   if isfield(IP,'P') & ~isempty(IP.P{1})
      IP.P = cell(1,4); 
   end
   if isfield(pm_defs,'magfieldmap')
      IP.fmagP=pm_defs.magfieldmap;
   end
otherwise 
   error('Funny number of input fieldmap images')
end

if ~isempty(IP.P{1})
%----------------------------------------------------------------------
% Create field map (in Hz) - this routine calls the unwrapping
%----------------------------------------------------------------------
   IP.fm = FieldMap('CreateFieldMap',IP);
%----------------------------------------------------------------------
% Write out field map
% Outputs -> fpm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------  
   FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Smoothed phase map');
end

%----------------------------------------------------------------------
% Convert Hz to voxels and write voxel displacement map 
% Outputs -> vdm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

[IP.vdm, IP.vdmP]=FieldMap('FM2VDM',IP);

%----------------------------------------------------------------------
% Select an EPI to unwarp
%----------------------------------------------------------------------

if ischar(epi_img)
    nsessions = 1;
    epi_img{1} = epi_img;
elseif iscell(epi_img)
    nsessions=size(epi_img,2);
else
    nsessions=0;
end

%----------------------------------------------------------------------
% Match voxel displacement map to image
% Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------
if nsessions==0
   VDM{1}=IP.vdmP;
   IPcell{1}=IP; 
elseif nsessions==1
    IP.epiP = spm_vol(epi_img{1}(1,:));
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
   VDM{1}=IP.vdmP;
   IPcell{1}=IP;
else
% If multiple sessions, does match to first image of each session
% Copies the written file to 
%----------------------------------------------------------------------
   orig_vdm=IP.vdmP;

   % get session specific fieldmap name
   if isfield(pm_defs,'sessname')
      sessname=pm_defs.sessname;
   else
      sessname='session';
   end

   for sessnum=1:nsessions
      IP.vdmP=orig_vdm; % Make sure we start with original for each session
      Ovdm=IP.vdmP;
      IP.epiP = spm_vol(epi_img{sessnum});
      if isfield(pm_defs, 'match_vdm')
         if pm_defs.match_vdm
            msg=sprintf('\nMatching session %d...\n',sessnum);
            disp(msg);
            IP.vdmP = FieldMap('MatchVDM',IP);
         end
         % Now copy this file to a session specific file
         session_vdm = spm_vol(IP.vdmP.fname);
         vol=spm_read_vols(session_vdm);
         vdm_info=sprintf('Voxel Displacement Map:echo time difference=%2.2fms, EPI readout time=%2.2fms',IP.uflags.etd, IP.tert);    
         newname=sprintf('%s_%s%d.%s',spm_str_manip(IP.vdmP.fname,'r'),sessname,sessnum,spm_str_manip(IP.vdmP.fname,'e'));
         Ovdm=struct('fname',newname,'mat',session_vdm.mat,'dim',session_vdm.dim,'dt',session_vdm.dt,'descrip',vdm_info);
         spm_write_vol(Ovdm,vol);
      end
  
      %----------------------------------------------------------------------
      % Unwarp EPI
      %----------------------------------------------------------------------
      if isfield(pm_defs,'write_unwarped')
         if pm_defs.write_unwarped
            IP.uepiP = FieldMap('UnwarpEPI',IP);

         %----------------------------------------------------------------------
         % Write unwarped EPI 
         % Outputs -> uNAME-OF-EPI.img
         %----------------------------------------------------------------------
            unwarp_info=sprintf('Unwarped EPI:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
            IP.uepiP = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info);
         end
      end
      % Update current values in memory
      IP.vdmP=Ovdm;
      VDM{sessnum}=IP.vdmP ;
      IPcell{sessnum}=IP;
   end  
end
return









