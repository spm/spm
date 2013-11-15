function varargout = cfg_util_persistent(varargin)
%CFG_UTIL_PERSISTENT - store persistent variables for cfg_util
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_util_persistent.m 5750 2013-11-15 15:02:24Z volkmar $

rev = '$Rev: 5750 $'; 

persistent c0;
persistent jobs;
if nargin == 2 && nargout == 0
    c0   = varargin{1};
    jobs = varargin{2};
elseif nargin == 0 && nargout == 2
    varargout{1} = c0;
    varargout{2} = jobs;
else
    cfg_message('matlabbatch:usage', '%s: Usage error.', mfilename);
end
