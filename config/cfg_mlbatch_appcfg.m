function [cfg, def, ver] = cfg_mlbatch_appcfg
% Add SPM to the application list of MATLABBATCH
% This file must be on MATLAB search path for cfg_util to detect it.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% $Id: cfg_mlbatch_appcfg.m 4964 2012-09-26 10:51:05Z guillaume $

cfg = spm_cfg;
def = [];
ver = spm('Version');
