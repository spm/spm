function warning_flexible(varargin)
% Function allowing to have better control over the warnings
% that might not be necessary at some point
% _______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: warning_flexible.m 5467 2013-05-05 20:03:50Z vladimir $

warning off backtrace
warning(varargin{:});
warning on backtrace