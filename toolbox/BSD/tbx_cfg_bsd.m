function bf = tbx_cfg_bsd
% Configuration file for toolbox 'Bayesian Spectral Decomposition'
%_______________________________________________________________________
% Copyright (C) 2024-2025 Wellcome Trust Centre for Neuroimaging

% Johan Medrano

tbxdir = fileparts(mfilename('fullpath'));

if ~isdeployed, addpath(tbxdir); end

components = {
    'spm_bsd_batch';
    };

bf = cfg_choice;
bf.tag = 'bsd';
bf.name = 'Bayesian Spectral Decomposition';
bf.help = {'Bayesian spectral decomposition toolbox'};

for i = 1:numel(components)
  bf.values{i} = feval(components{i});
end


