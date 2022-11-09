function tsss = tbx_cfg_tsss
% Configuration file for toolbox 'TSSS'
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2014-2022 Wellcome Centre for Human Neuroimaging


if ~isdeployed, addpath(fileparts(mfilename('fullpath'))); end

tsss        = cfg_choice;
tsss.tag    = 'tsss';
tsss.name   = 'TSSS';
tsss.help   = {'Temporal Signal Space Separation (TSSS) toolbox'};
tsss.values = {tsss_config, tsss_config_momentspace};
