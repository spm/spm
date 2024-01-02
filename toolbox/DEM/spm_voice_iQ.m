function [W] = spm_voice_iQ(Q)
% Discrete cosine transform of formant coefficients
% FORMAT [W] = spm_voice_iQ(Q)
%
% Q     - log formant frequencies
% G(1)  - log formant (pitch) Tu
% G(2)  - log timing  (pitch) Tv

%
% W     - (Nu x Nv) log formant coefficients (weights)
%
%   Nu  - number of formant coefficients
%   Nv  - number of timing  coefficients
%   Tu  - log formant (pitch)
%   Tv  - log timing  (pitch) 
%
% This  auxiliary routine scales and transforms log formant coefficients
% using a pair of discrete cosine transforms with logarithmic scaling.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% defaults and (logarithmic) scaling
%--------------------------------------------------------------------------
try, Nu = VOX.Nu; catch, Nu  = 32;         end    % DCT order   (formants)
try, Nv = VOX.Nv; catch, Nv  = 8;          end    % DCT order   (interval)
try, Tu = VOX.Tu; catch, Tu  = 4;          end    % log scaling (formants)
try, Tv = VOX.Tv; catch, Tv  = 1;          end    % log scaling (interval)

% sizes 
%--------------------------------------------------------------------------
[Ni,ni] = size(Q);

%  inverse transform
%--------------------------------------------------------------------------
U  = spm_voice_dct(Ni,Nu,Tu);                % DCT over formants
V  = spm_voice_dct(ni,Nv,Tv);                % DCT over intervals
W  = U\Q/V';                                 % coefficients
W  = W/std(W(:));                            % normalise
