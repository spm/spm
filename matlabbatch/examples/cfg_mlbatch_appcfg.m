function [cfg, def] = cfg_mlbatch_appcfg(varargin)

% function [cfg, def] = cfg_mlbatch_appcfg(varargin)
% Application initialisation for matlabbatch example application. This
% application only has a cfg_ file, no defaults. The main application
% body will just add the path to the example directory to MATLAB path and
% call cfg_ui. This will initialise cfg_util with all cfg_mlbatch_appcfg
% files found on the MATLAB path.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_mlbatch_appcfg.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

% toy example configuration
cfg = cfg_example_master;
def = [];