function [cfg, def, ver] = cfg_mlbatch_appcfg
% Add SPM to the application list of MATLABBATCH
% This file must be on MATLAB search path for cfg_util to detect it.
%__________________________________________________________________________

% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


cfg = spm_cfg;
def = [];
ver = spm('Version');
