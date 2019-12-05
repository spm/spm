function [Q,U,V] = spm_voice_Q(W,G,Ni,ni)
% Inverse discrete cosine transform of formant coefficients
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
% using a pair of discrete cosine transforms with logarithmic scaling.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_Q.m 7750 2019-12-05 17:54:29Z spm $


% defaults and (logarithmic) scaling
%--------------------------------------------------------------------------
try, Tu = VOX.Tu; catch, Tu  = 4; end    % log scaling (formants)
try, Tv = VOX.Tv; catch, Tv  = 1; end    % log scaling (interval)
try, Tf = VOX.Tf; catch, Tf  = 1; end    % lin scaling (formants)
try, Tw = VOX.Tv; catch, Tw  = 1; end    % lin scaling (amplitude)

if nargin < 3, Ni = 256; end
if nargin < 4, ni = 64;  end
if nargin > 1
    Tu = Tu*exp(G(1));
    Tv = Tv*exp(G(2));
    Tf = Tf*exp(G(3));
    Tw =     Tw*G(4);

end

% sizes
%--------------------------------------------------------------------------
[Nu,Nv] = size(W);

%  inverse transform
%--------------------------------------------------------------------------
U  = spm_voice_dct(Ni,Nu,Tu,Tf);             % DCT over formants
V  = spm_voice_dct(ni,Nv,Tv);                % DCT over intervals
Q  = U*W*V';                                 % log formants
A  = 1 - exp(-(1:Ni)*8/Ni)*Tw;               % amplitude modulation
Q  = bsxfun(@times,Q,A(:));                  % balance
Q  = Q/std(Q(:));                            % normalise
