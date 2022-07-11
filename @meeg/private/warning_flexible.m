function warning_flexible(varargin)
% Function allowing to have better control over the warnings
% that might not be necessary at some point
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


warning off backtrace
warning(varargin{:});
warning on backtrace
