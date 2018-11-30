function out = double(fa)
% Convert to double precision
% FORMAT double(fa)
% fa - a file_array
%__________________________________________________________________________
% Copyright (C) 2005-2018 Wellcome Trust Centre for Neuroimaging

%
% $Id: double.m 7501 2018-11-30 12:16:58Z guillaume $


out = double(full(fa));
