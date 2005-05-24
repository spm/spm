function out = double(fa)
% Convert to double precision
% FORMAT double(fa)
% fa - a file_array
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: double.m 174 2005-05-24 11:03:32Z john $

out = double(numeric(fa));

