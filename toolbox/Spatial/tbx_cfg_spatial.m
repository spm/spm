function spatial = tbx_cfg_spatial
% MATLABBATCH Configuration file for toolbox 'Spatial Tools'

% John Ashburner
% $Id: tbx_cfg_spatial.m 8053 2021-02-08 11:51:17Z spm $

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

