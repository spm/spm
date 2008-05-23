function varargout = cat(varargin)

% function varargout = cat(varargin)
% Prevent cat for cfg_item objects.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cat.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if nargin > 2
    error('matlabbatch:cfg_item:cat', ['Concatenation of cfg_item objects is ' ...
                    'not allowed.']);
else
    varargout{1} = varargin{2};
end;
