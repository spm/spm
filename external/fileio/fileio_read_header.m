function varargout = funname(varargin);

% this is a SPM wrapper around a FieldTrip function 

% this part is variable
prefix = 'fileio_';

% this part is fixed
funname   = mfilename;
funhandle = str2func(funname((length(prefix)+1):end));
[varargout{1:nargout}] = funhandle(varargin{:});
