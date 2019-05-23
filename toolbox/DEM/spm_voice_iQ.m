function [W] = spm_voice_iQ(Q,G,Nu,Nv)
% discrete cosine transform of formant coefficients
% FORMAT [W] = spm_voice_iQ(Q,G,Nu,Nv)
%
% Q     - log formant frequencies
% G(1)  - log formant (pitch) Tu
% G(2)  - log timing  (pitch) Tv
% Nu    - number of formant coefficients
% Nv    - number of timing  coefficients
%
% W     - log formant coeficients (weights)
%
% This  auxiliary routine scales and transforms log formant coefficients
% using a pair of discrete cosine transforms with logarithmic scaling
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_iQ.m 7597 2019-05-23 18:42:38Z karl $


% defaults and (logarithmic) scaling
%--------------------------------------------------------------------------
if nargin < 3, Nu = 32; end
if nargin < 4, Nv = 8;  end
if nargin < 2
    Tu = 4;
    Tv = 1;
else
    Tu = exp(G(1));
    Tv = exp(G(2));
end

% sizes 
%--------------------------------------------------------------------------
[Ni,ni] = size(Q);

%  inverse transform
%--------------------------------------------------------------------------
Q  = Q/std(Q(:));                            % normalise
U  = spm_voice_dct(Ni,Nu,Tu);                % DCT over formants
V  = spm_voice_dct(ni,Nv,Tv);                % DCT over intervals
W  = U\Q/V';                                 % coeficients


