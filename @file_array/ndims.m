function out = ndims(fa)
% Number of dimensions
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: ndims.m 7147 2017-08-03 14:07:01Z spm $


out = size(fa);
out = length(out);
