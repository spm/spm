function [cfg, def] = cfg_mlbatch_appcfg(varargin)

% Add BasicIO to applications list of cfg_util. This file is an example how
% to add your own application configuration to cfg_util. To add an
% application, create a file called cfg_mlbatch_appcfg.m in the application
% folder and add this folder to the MATLAB path. cfg_util will look for
% files with the exact name cfg_mlbatch_appcfg.m and run all of them in
% order of their occurence on the path.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_mlbatch_appcfg.m 3179 2009-06-03 12:41:21Z volkmar $

rev = '$Rev: 3179 $'; %#ok

if ~isdeployed
    % get path to this file
    p = fileparts(mfilename('fullpath'));
    % in this case, the cfg_mlbatch_appcfg file is not in the application
    % folder. Therefore, the path to the application needs to be added
    % explicitly
    addpath(fullfile(p, 'cfg_basicio'));
end
% these two files are now on MATLABs path
cfg = cfg_cfg_basicio;
def = cfg_cfg_basicio_def;
