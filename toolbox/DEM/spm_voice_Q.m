function [Q,U,V] = spm_voice_Q(W,G,Ni,ni)
% inverse discrete cosine transform of formant coefficients
% FORMAT [Q,U,V] = spm_voice_Q(W,G,Ni,ni)
%
% W     - log formant coefficients (weights)
% G(1)  - log formant (pitch) Tu
% G(2)  - log timing  (pitch) Tv
% G(3)  - amplitude   (pitch) Tw
% Ni    - number of formant frequencies
% ni    - number of timing  intervals
%
% Q     - formants (time-frequency representation): Q = U*xY.W*V'
% U     - DCT over frequency
% V     - DCT over intervals
%
% This  auxiliary routine scales and transforms log formant coefficients
% using a pair of discrete cosine transforms with logarithmic scaling
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_Q.m 7597 2019-05-23 18:42:38Z karl $

% defaults and (logarithmic) scaling
%--------------------------------------------------------------------------
if nargin < 3, Ni = 256; end
if nargin < 4, ni = 64;  end
if nargin < 2
    Tu = 4;
    Tv = 1;
else
    Tu = exp(G(1));
    Tv = exp(G(2));
    Tw =     G(3);
end


% sizes
%--------------------------------------------------------------------------
[Nu,Nv] = size(W);

%  inverse transform
%--------------------------------------------------------------------------
U  = spm_voice_dct(Ni,Nu,Tu);                % DCT over formants
V  = spm_voice_dct(ni,Nv,Tv);                % DCT over intervals
Q  = U*W*V';                                 % log formants
A  = 1 - exp(-(1:Ni)*8/Ni)*Tw;               % amplitude modulation
Q  = bsxfun(@times,Q,A(:));                  % balance
Q  = Q/std(Q(:));                            % normalise




