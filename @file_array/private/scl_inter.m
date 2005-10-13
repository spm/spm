function varargout = scl_inter(varargin)
% Format
% For getting the value
% dat = scl_inter(obj)
%
% For setting the value
% obj = scl_inter(obj,dat)
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: scl_inter.m 253 2005-10-13 15:31:34Z guillaume $



if nargin==2,
    varargout{1} = asgn(varargin{:});
elseif nargin==1,
    varargout{1} = ref(varargin{:});
else
    error('Wrong number of arguments.');
end;
return;

function dat = ref(obj)
dat = obj.scl_inter;
return;

function obj = asgn(obj,dat)
if isnumeric(dat), % && numel(dat)<=1,
    obj.scl_inter = double(dat);
else
    error('"scl_inter" must be numeric.');
end;
return;

