function res = dtype(obj)
% returns datatype of embedded file_array object
% FORMAT dtype(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel


res = obj.data.y.dtype;