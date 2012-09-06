function [cfg, def] = cfg_mlbatch_appcfg
% Add SPM to the application list of MATLABBATCH
% This file must be on MATLAB search path for cfg_util to detect it.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% $Id: cfg_mlbatch_appcfg.m 4904 2012-09-06 15:08:56Z guillaume $

cfg = spm_cfg;
def = [];
