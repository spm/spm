function H = spm_hanning(n)
% Return the n-point Hanning window in a column vector
% FORMAT H = spm_hanning(n)
% n  -  length of hanning function
% H  -  hanning function
%__________________________________________________________________________
% Copyright (C) 2007-2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hanning.m 7749 2019-12-05 17:05:46Z guillaume $


H  = (1 - cos(2*pi*[1:n]'/(n + 1)))/2;
