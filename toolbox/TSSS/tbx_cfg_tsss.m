function tsss = tbx_cfg_tsss
% Configuration file for toolbox 'TSSS'
%__________________________________________________________________________
% Copyright (C) 2014-2021 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: tbx_cfg_tsss.m 8119 2021-07-06 13:51:43Z guillaume $


if ~isdeployed, addpath(fileparts(mfilename('fullpath'))); end

tsss        = cfg_choice;
tsss.tag    = 'tsss';
tsss.name   = 'TSSS';
tsss.help   = {'Temporal Signal Space Separation (TSSS) toolbox'};
tsss.values = {tsss_config, tsss_config_momentspace};
