function warning_flexible(varargin)
% Function allowing to have better control over the warnings
% that might not be necessary at some point
% _______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: warning_flexible.m 5464 2013-05-03 15:20:14Z vladimir $

%warning off backtrace
warning(varargin{:});
%warning on backtrace