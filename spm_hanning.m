function H = spm_hanning(n)
% returns the n-point Hanning window in a column vector
% FORMAT H = spm_hanning(n);
% n  -  length of hanning function
% H  -  hanning function
%___________________________________________________________________________


H  = (1 - cos(2*pi*[1:n]'/(n + 1)))/2;
