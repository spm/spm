function varargout = funname(varargin);

% this is a SPM wrapper around a FieldTrip function 

% this part is variable
funhandle = @prepare_localspheres;

% this part is fixed
[varargout{1:nargout}] = funhandle(varargin{:});
