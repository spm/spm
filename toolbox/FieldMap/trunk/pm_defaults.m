%
% Sets the default values for the FieldMap toolbox
%
% FORMAT pm_defaults
%_______________________________________________________________________
%
% This file is intended for site-specific customisations.
%
%_______________________________________________________________________
% pm_defaults Chloe Hutton, Jesper Andersson 20/01/04

global pm_def

% Defaults for creating field map. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%=======================================================================
pm_def.INPUT_DATA_FORMAT = 'RI';      % 'RI' = load two real and 
                                      % imaginary image pairs
                                      % 'PM' = load one or two
                                      % phase and magnitude image
                                      % pairs.
pm_def.SHORT_ECHO_TIME = 25;          % Short echo time in ms
pm_def.LONG_ECHO_TIME = 34.5;         % Long echo time in ms

pm_def.MASK = 0;                      % Do brain masking (1 or 0 - unnecessary for EPI fieldmaps)? 

% Defaults for unwrapping options. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%=======================================================================
pm_def.UNWRAPPING_METHOD = 'Mark3D';  % Unwrapping options are:
                                      % 'Huttonish', 'Mark3D' or 'Mark2D'
pm_def.FWHM = 10;                     % FWHM (mm) of Gaussian filter used to 
                                      % implement weighted smoothing of
                                      % unwrapped maps.
pm_def.PAD = 0;                       % Size of padding kernel if required.
pm_def.WS = 1;                        % Weighted or normal smoothing.
pm_def.BMASK = [];                    % Mask used for unwrapping

% Flags for brain extraction
%=======================================================================
pm_def.MFLAGS.TEMPLATE = fullfile(spm('Dir'),'templates','T1.mnc');
pm_def.MFLAGS.FWHM = 5; % In mm
pm_def.MFLAGS.NERODE = 1;% In voxels
pm_def.MFLAGS.NDILATE = 2; % In voxels
pm_def.MFLAGS.THRESH = 0.5;

% Defaults for converting field map to voxel displacement map.
%=======================================================================
pm_def.EPI_BASED_FIELDMAPS = 1;        % EPI=1, other=0.
pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = 1; % +ve k-space = 1, -ve = -1.
pm_def.TOTAL_EPI_READOUT_TIME = 32;    % Sonata EPI RO time (500E-6*64)

% Defaults for Unwarping.
%=======================================================================
pm_def.DO_JACOBIAN_MODULATION = 0;    % Do jacobian modulation to adjust 
                                      % for compression or stretching
                                      % No = 0, Yes = 1

% FIL specific additions
%=======================================================================

global SCANNER
global SEQUENCE

if findstr(SCANNER,'Sonata') & findstr(SEQUENCE,'Siemens')
   disp('Using Sonata Siemens parameters'); 
   pm_def.INPUT_DATA_FORMAT = 'PM'; 
   pm_def.SHORT_ECHO_TIME = 10.0; 
   pm_def.LONG_ECHO_TIME = 14.76;
   pm_def.EPI_BASED_FIELDMAPS = 0;
   pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = 1;
   pm_def.MASK = 1; 

elseif findstr(SCANNER, 'Sonata') & findstr(SEQUENCE,'Flash');
    disp('Using Sonata Flash parameters');
    pm_def.SHORT_ECHO_TIME = 4; 
    pm_def.LONG_ECHO_TIME = 13;
    pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = 1;
    pm_def.EPI_BASED_FIELDMAPS = 0;

elseif findstr(SCANNER, 'Allegra') 
   pm_def.TOTAL_EPI_READOUT_TIME = 21.0; 

   if findstr(SEQUENCE,'Siemens')
     disp('Using Allegra Siemens parameters');
      pm_def.INPUT_DATA_FORMAT = 'PM';
      pm_def.SHORT_ECHO_TIME = 10.0; % 
      pm_def.LONG_ECHO_TIME = 12.46;
      pm_def.EPI_BASED_FIELDMAPS = 0;
      pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = 1;
      pm_def.WS = 1;
      pm_def.PAD = 0;
      pm_def.MASK = 1; 
   else
      disp('Using Allegra EPI parameters');       
      pm_def.SHORT_ECHO_TIME = 19;
      pm_def.LONG_ECHO_TIME = 29;
      pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = 1; % This is 1 for FLASH and EPI
      pm_def.EPI_BASED_FIELDMAPS = 1;
   end   
end
