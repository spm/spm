function [U] = spm_voice_dct(N,K,n)
% logarithmically sampled discrete cosine transform matrix
% FORMAT [U] = spm_voice_dct(N,K,n)
%
%  N        - dimension
%  K        - order
%  n        - log scaling parameter (typically 4)
%
% This routine returns a discrete cosine transform matrix sampled
% logarithmically according to a scaling parameter.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_dct.m 7597 2019-05-23 18:42:38Z karl $


% discrete cosine transform matrix
%--------------------------------------------------------------------------
nu = exp(-n*(0:(N - 1))/N);                  % log spacing
nu = nu - min(nu);                           % in short range lies between
nu = N - (N - 1)*nu/max(nu);                 % 0 and N
U  = spm_dctmtx(N,K,nu);                     % DCT matrix

return



