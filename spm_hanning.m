function H = spm_hanning(n)
% Return the n-point Hanning window in a column vector
% FORMAT H = spm_hanning(n)
% n  -  length of hanning function
% H  -  hanning function
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2007-2022 Wellcome Centre for Human Neuroimaging


H  = (1 - cos(2*pi*[1:n]'/(n + 1)))/2;
