function varargout = funname(varargin);

% this is a SPM wrapper around a FieldTrip function 

% this part is variable
funhandle = @trialfun_general;

% this part is fixed
[varargout{1:nargout}] = funhandle(varargin{:});
