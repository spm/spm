function varargout = fname(varargin)
% file_array's fname property
% For getting the value
% dat = fname(obj)
%
% For setting the value
% obj = fname(obj,dat)
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: fname.m 7147 2017-08-03 14:07:01Z spm $



if nargin==2
    varargout{1} = asgn(varargin{:});
elseif nargin==1
    varargout{1} = ref(varargin{:});
else
    error('Wrong number of arguments.');
end


%==========================================================================
% function dat = ref(obj)
%==========================================================================
function dat = ref(obj)
dat = obj.fname;


%==========================================================================
% function obj = asgn(obj,dat)
%==========================================================================
function obj = asgn(obj,dat)
if ischar(dat)
    obj.fname = deblank(dat(:)');
else
    error('"fname" must be a character string.');
end