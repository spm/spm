function out = numeric(fa)
% Convert to numeric form
% FORMAT numeric(fa)
% fa - a file_array
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: numeric.m 174 2005-05-24 11:03:32Z john $


[vo{1:ndims(fa)}] = deal(':');
out = subsref(fa,struct('type','()','subs',{vo}));

