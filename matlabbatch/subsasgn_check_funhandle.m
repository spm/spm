function sts = subsasgn_check_funhandle(val)

% function sts = subsasgn_check_funhandle(val)
% Return true if val is either empty, or a function or function handle.
% One could also check for nargin == 1 and nargout == 1.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check_funhandle.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

sts = isempty(val) || isa(val, 'function_handle') || (ischar(val) && any(exist(val) == 2:6));