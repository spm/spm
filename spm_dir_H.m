function H = spm_dir_H(a)
% Entropy of a Dirichlet distribution
% FORMAT H = spm_dir_H(a)
% a    - Dirichlet parameters of a joint distribution
% H    - Entropy
%
% The entropy here pertains to the Dirichlet distribution.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% entropy of a (conditioned) Dirichlet distribution 
%--------------------------------------------------------------------------
a   = a + exp(-16);                                 % preclude overflow
a0  = sum(a);
k   = size(a,1);
H   = spm_betaln(a) + (a0 - k).*psi(a0) - sum((a - 1).*psi(a));
