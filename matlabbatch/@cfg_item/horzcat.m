function varargout = horzcat(varargin)

% function varargout = horzcat(varargin)
% Prevent horzcat for cfg_item objects.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: horzcat.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

error('matlabbatch:cfg_item:cat', ['Concatenation of cfg_item objects is ' ...
                    'not allowed.']);