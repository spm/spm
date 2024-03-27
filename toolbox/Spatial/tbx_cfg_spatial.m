function spatial = tbx_cfg_spatial
% Configuration file for toolbox 'Spatial Tools'
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Spatial')); end

spatial        = cfg_choice;
spatial.tag    = 'spatial';
spatial.name   = 'Spatial Tools';
spatial.values = {spm_slice2vol_config spm_scope_config spm_TVdenoise_config};
spatial.help   = {[...
'A selection of work-in-progress tools for various spatial processing tasks.',...
]};
