function res = dtype(obj)
% returns datatype of embedded file_array object
% FORMAT dtype(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: dtype.m 1373 2008-04-11 14:24:03Z spm $

res = obj.data.y.dtype;
