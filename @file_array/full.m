function out = full(fa)
% Convert to numeric form
% FORMAT full(fa)
% fa - a file_array
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

%
% $Id: full.m 7501 2018-11-30 12:16:58Z guillaume $


[vo{1:ndims(fa)}] = deal(':');
out = subsref(fa,struct('type','()','subs',{vo}));
