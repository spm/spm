function out = numeric(fa)
% Convert to numeric form
% FORMAT numeric(fa)
% fa - a file_array
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: numeric.m 253 2005-10-13 15:31:34Z guillaume $


[vo{1:ndims(fa)}] = deal(':');
out = subsref(fa,struct('type','()','subs',{vo}));

