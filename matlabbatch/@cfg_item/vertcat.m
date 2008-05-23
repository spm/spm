function varargout = vertcat(varargin)

% function varargout = vertcat(varargin)
% Prevent vertcat for cfg_item objects.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: vertcat.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

error('matlabbatch:cfg_item:cat', ['Concatenation of cfg_item objects is ' ...
                    'not allowed.']);