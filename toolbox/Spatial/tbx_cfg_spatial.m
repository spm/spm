function spatial = tbx_cfg_spatial
% Configuration file for toolbox 'Spatial Tools'
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% John Ashburner
% $Id: tbx_cfg_spatial.m 8119 2021-07-06 13:51:43Z guillaume $

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Spatial')); end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
spatial        = cfg_choice;
spatial.tag    = 'spatial';
spatial.name   = 'Spatial Tools';
spatial.values = {spm_slice2vol_config};
spatial.help   = {[...
'A selection of work-in-progress tools for various spatial processing tasks.',...
]};
%==========================================================================

