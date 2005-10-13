function varargout = scl_slope(varargin)
% Format
% For getting the value
% dat = scl_slope(obj)
%
% For setting the value
% obj = scl_slope(obj,dat)
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: scl_slope.m 253 2005-10-13 15:31:34Z guillaume $



if nargin==2,
    varargout{1} = asgn(varargin{:});
elseif nargin==1,
    varargout{1} = ref(varargin{:});
else
    error('Wrong number of arguments.');
end;
return;

function dat = ref(obj)
dat = obj.scl_slope;
return;

function obj = asgn(obj,dat)
if isnumeric(dat), % && numel(dat)<=1,
    obj.scl_slope = double(dat);
else
    error('"scl_slope" must be numeric.');
end;
return;
