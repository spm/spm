function out = double(fa)
% Convert to double precision
% FORMAT double(fa)
% fa - a file_array
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: double.m 253 2005-10-13 15:31:34Z guillaume $

out = double(numeric(fa));

