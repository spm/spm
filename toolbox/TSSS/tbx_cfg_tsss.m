function tsss = tbx_cfg_tsss
% Configuration file for toolbox 'TSSS'
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: tbx_cfg_tsss.m 7703 2019-11-22 12:06:29Z guillaume $

tbxdir = fileparts(mfilename('fullpath'));

if ~isdeployed, addpath(tbxdir); end

components = {   
    'tsss_config'
    'tsss_config_momentspace'
    };

tsss = cfg_choice;
tsss.tag  = 'tsss';
tsss.name = 'TSSS';
tsss.help = {'Temporal Signal Space Separation (TSSS) toolbox'};

for i = 1:numel(components)
  tsss.values{i} = feval(components{i});
end

