function varargout = spm_subfun(varargin)
% Enable calling local functions
% FORMAT [o1,o2,...] = spm_subfun(localfunctions,action,i1,i2,...)
% The funcion is supposed to be inserted into multifunction m-files
% so that it calls localfunctions within the scope of the m-file.
% The output of this is used to match the action string with the
% name of each local function to see which of them to call.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% $Id$
if nargin<=1
    [varargout{1:nargout}] = import(varargin{1:nargin});
else
    [varargout{1:nargout}] = select(varargin{1:nargin});
end
%==========================================================================

%==========================================================================
function varargout = select(funs,opt,varargin)
opt = lower(opt);
s   = import(funs);
if ~isfield(s,opt), error('Unknown function (%s)',opt); end
[varargout{1:nargout}] = feval(s.(opt),varargin{:});
%==========================================================================

%==========================================================================
function varargout = import(funs,varargin)
names = cellfun(@(x)lower(char(x)),funs,'UniformOutput',false);
c     = {names{:}; funs{:}};
varargout{1} = struct(c{:});
%==========================================================================

