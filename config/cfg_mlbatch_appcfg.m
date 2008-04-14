function [cfg, def] = cfg_mlbatch_appcfg

% MATLABBATCH initialisation function for SPM. This file must be on
% MATLABs path for cfg_util to pick up SPM configuration.
%_______________________________________________________________________
% SPM is developed by members and collaborators of the
% Wellcome Trust Centre for Neuroimaging

%-$Id$
%-----------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

rev = '$Rev: 207 $';

cfg  = spm_cfg;
def  = spm_def;