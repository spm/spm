function out = double(fa)
% Convert to double precision
% FORMAT double(fa)
% fa - a file_array
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: double.m 7147 2017-08-03 14:07:01Z spm $


out = double(numeric(fa));
