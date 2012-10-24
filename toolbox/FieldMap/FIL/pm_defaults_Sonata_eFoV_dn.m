% SPM5 UPDATE 14/11/06
% UPDATE 27/01/05
% Sets the default values for the FieldMap toolbox
%
% FORMAT pm_defaults_Sonata_eFoV_dn
%_______________________________________________________________________
%
% This file is intended for use with the Siemens fieldmap sequence
% on the Sonata scanner at the FIL, with the eFoV sequence and PE blips=+1:
% nw_mepi_v3d_OFC
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton and Jesper Andersson
% $Id: pm_defaults_Sonata_eFoV_dn.m 5015 2012-10-24 13:40:07Z guillaume $

global pm_def

% Defaults for creating field map. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%=======================================================================
pm_def.INPUT_DATA_FORMAT = 'PM';      % 'RI' = load two real and 
                                      % imaginary image pairs
                                      % 'PM' = load one or two
                                      % phase and magnitude image
                                      % pairs.
pm_def.SHORT_ECHO_TIME = 10.0;        % Short echo time in ms for Sonata
pm_def.LONG_ECHO_TIME = 14.76;        % Long echo time in ms for Sonata
pm_def.MASKBRAIN = 1;                 % Do brain masking (1 or 0,
                      % 0 for EPI fieldmaps)

% Defaults for unwrapping options. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%=======================================================================
pm_def.UNWRAPPING_METHOD = 'Mark3D';  % Unwrapping options are:
                                      % 'Huttonish', 'Mark3D' or 'Mark2D'
pm_def.FWHM = 10;                     % FWHM of Gaussian filter used to 
                                      % implement weighted smoothing of
                                      % unwrapped maps.
pm_def.PAD = 0;                       % Size of padding kernel if required.
pm_def.WS = 1;                        % Weighted or normal smoothing.

% Flags for brain extraction
%=======================================================================
pm_def.MFLAGS.TEMPLATE = fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii');
pm_def.MFLAGS.FWHM = 5; % In mm
pm_def.MFLAGS.NERODE = 2;% In voxels
pm_def.MFLAGS.NDILATE = 4; % In voxels
pm_def.MFLAGS.THRESH = 0.5;
pm_def.MFLAGS.REG = 0.02; % A larger value helps segmentation to converge
pm_def.MFLAGS.GRAPHICS = 0; % A larger value helps segmentation to converge

% Defaults for converting field map to voxel displacement map.
%=======================================================================
pm_def.EPI_BASED_FIELDMAPS = 0;         % EPI=1, other=0.
pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = +1; % +ve k-space = 1, -ve = -1.
pm_def.TOTAL_EPI_READOUT_TIME = 36;     % Sonata EPI RO time (500E-6*72)

% Defaults for Unwarping.
%=======================================================================
pm_def.DO_JACOBIAN_MODULATION = 0;    % Do jacobian modulation to adjust 
                                      % for compression or stretching
                                      % No = 0, Yes = 1

