function [U] = spm_voice_dct(N,K,n,s)
% Logarithmically sampled discrete cosine transform matrix
% FORMAT [U] = spm_voice_dct(N,K,n,[s])
%
%  N         - dimension
%  K         - order
%  n         - log scaling parameter (typically 4)
%  s         - optional linear scaling [default 1]
%
% This routine returns a discrete cosine transform matrix sampled
% logarithmically according to a scaling parameter.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 4, s = 1; end

% discrete cosine transform matrix
%--------------------------------------------------------------------------
nu = exp(-n*(0:(N - 1))/N);                  % log spacing
nu = 1 - nu;                                 % ensure range lies between
nu = s*N*nu/nu(end);                         % 0 and s*N
U  = spm_dctmtx(N,K,nu);                     % DCT matrix
