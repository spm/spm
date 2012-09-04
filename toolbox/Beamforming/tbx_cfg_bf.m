function bf = tbx_cfg_bf
% Configuration file for toolbox 'Beamforming'
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: tbx_cfg_bf.m 4897 2012-09-04 16:32:18Z vladimir $

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Beamforming')); end

components = {
    'bf_data';
    'bf_sources'
    'bf_features'
    'bf_inverse'
    'bf_output'
    'bf_write'
    };

bf = cfg_choice;
bf.tag = 'beamforming';
bf.name = 'Beamforming';
bf.help = {'Beamforming toolbox'};

for i = 1:numel(components)
  bf.values{i} = feval(components{i});
end

