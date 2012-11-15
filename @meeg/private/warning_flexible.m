function warning_flexible(varargin)
% Function allowing to have better control over the warnings
% that might not be necessary at some point
% _______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: warning_flexible.m 5057 2012-11-15 13:03:35Z vladimir $

warning off backtrace
warning(varargin{:});
warning on backtrace