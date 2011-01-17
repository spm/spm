function varargout = cfg_callbuiltin(varargin)
if nargout == 0
    try
        varargout{1} = builtin(varargin{:});
    catch
        builtin(varargin{:});
    end
else
    varargout = cell(1,nargout);
    [varargout{:}] = builtin(varargin{:});
end
