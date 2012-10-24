% Sets the default values for the FieldMap toolbox
%
% FORMAT pm_defaults
%
%__________________________________________________________________________
%
% This file is intended for site and/or scanner and/or
% sequence-specific customisations.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton and Jesper Andersson
% $Id: pm_defaults.m 5014 2012-10-24 10:56:28Z guillaume $

global pm_def

% Defaults for creating field map. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%==========================================================================
pm_def.INPUT_DATA_FORMAT = 'PM';      % 'RI' = load two real and 
                                      % imaginary image pairs
                                      % 'PM' = load one or two
                                      % phase and magnitude image
                                      % pairs.
pm_def.SHORT_ECHO_TIME = 10.0;        % Short echo time in ms for Allegra
pm_def.LONG_ECHO_TIME = 12.46;        % Long echo time in ms for Allegra
pm_def.MASKBRAIN = 0;                 % Do brain masking (1 or 0,
                                      % 0 for EPI fieldmaps)

% Defaults for unwrapping options. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%==========================================================================
pm_def.UNWRAPPING_METHOD = 'Mark3D';  % Unwrapping options are:
                                      % 'Huttonish', 'Mark3D' or 'Mark2D'
pm_def.FWHM = 10;                     % FWHM of Gaussian filter used to 
                                      % implement weighted smoothing of
                                      % unwrapped maps.
pm_def.PAD = 0;                       % Size of padding kernel if required.
pm_def.WS = 1;                        % Weighted or normal smoothing.

% Flags for brain extraction
%==========================================================================
% Default T1 template for segmentation
pm_def.MFLAGS.TEMPLATE = fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii'); 
pm_def.MFLAGS.FWHM = 5;     % {mm} used for smoothing to fill holes in brain mask
pm_def.MFLAGS.NERODE = 2;   % {vox} used for erosion to remove scalp from brain mask
pm_def.MFLAGS.NDILATE = 4;  % {vox} used for dilation to condition scalp removal
pm_def.MFLAGS.THRESH = 0.5; % Intensity thresholding for filling holes
pm_def.MFLAGS.REG = 0.02;   % A larger value helps segmentation to converge
pm_def.MFLAGS.GRAPHICS = 0; % Don't display segmentation results

% Defaults for converting field map to voxel displacement map.
%==========================================================================
pm_def.EPI_BASED_FIELDMAPS = 0;         % EPI=1, other=0.
pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = -1; % +ve k-space = 1, -ve = -1.
pm_def.TOTAL_EPI_READOUT_TIME = 21.1;   % Allegra EPI RO time (330E-6*64)

% Defaults for Unwarping.
%==========================================================================
pm_def.DO_JACOBIAN_MODULATION = 0;      % Do jacobian modulation to adjust 
                                        % for compression or stretching
                                        % No = 0, Yes = 1

                                      
% The rest is dead code that I leave in for "documentation"
% purposes.  JeAn 161206

% FIL specific additions 
%==========================================================================
%
% global SCANNER
% global SEQUENCE
% 
% if findstr(SCANNER,'Sonata') & findstr(SEQUENCE,'Siemens')
%    pm_def.TOTAL_EPI_READOUT_TIME = 32;    % Sonata EPI RO time (500E-6*64)
%    pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = -1;
%    if findstr(SEQUENCE,'Siemens')
%       disp('Using Sonata Siemens parameters'); 
%       pm_def.INPUT_DATA_FORMAT = 'PM'; 
%       pm_def.SHORT_ECHO_TIME = 10.0; 
%       pm_def.LONG_ECHO_TIME = 14.76;
%       pm_def.EPI_BASED_FIELDMAPS = 0;
%    else
%       disp('Using Sonata EPI parameters');       
%       pm_def.SHORT_ECHO_TIME = 25;
%       pm_def.LONG_ECHO_TIME = 34.5;
%       pm_def.EPI_BASED_FIELDMAPS = 1;
%    end
% elseif findstr(SCANNER, 'Allegra') 
%    pm_def.TOTAL_EPI_READOUT_TIME = 21.1;  % Allegra EPI RO time (330E-6*64)
%    pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = -1;
%    if findstr(SEQUENCE,'Siemens')
%       disp('Using Allegra Siemens parameters');
%       pm_def.INPUT_DATA_FORMAT = 'PM';
%       pm_def.SHORT_ECHO_TIME = 10.0; % 
%       pm_def.LONG_ECHO_TIME = 12.46;
%       pm_def.EPI_BASED_FIELDMAPS = 0;
%    else
%       disp('Using Allegra EPI parameters');   
%       pm_def.INPUT_DATA_FORMAT = 'RI';
%       pm_def.SHORT_ECHO_TIME = 19;
%       pm_def.LONG_ECHO_TIME = 29;
%       pm_def.EPI_BASED_FIELDMAPS = 1;
%    end
% end
