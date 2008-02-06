function H = spm_hanning(n)
% returns the n-point Hanning window in a column vector
% FORMAT H = spm_hanning(n);
% n  -  length of hanning function
% H  -  hanning function
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hanning.m 1131 2008-02-06 11:17:09Z spm $

H  = (1 - cos(2*pi*[1:n]'/(n + 1)))/2;
