function res = dtype(obj, value)
% returns datatype of embedded file_array object
% FORMAT dtype(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: dtype.m 4432 2011-08-15 12:43:44Z christophe $

if nargin == 1
    res = obj.data.y.dtype;
else
    obj.data.y.dtype = value;
    obj.data.datatype = value;
    res = obj;
end