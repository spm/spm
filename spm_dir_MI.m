function E = spm_dir_MI(a)
% Expected information gain (i.e., mutual information)
% FORMAT E = spm_dir_MI(a)
%
% a    - Dirichlet parameters of a joint distribution
% C    - log preferences
%
% E    - expected free energy (information gain minus cost)
%
% The mutual information here pertains to the Dirichlet distribution. See
% spm_MDP_MI for the mutual information of the expected categorical
% distribution.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging


% deal with tensors
%--------------------------------------------------------------------------
a     = a(:,:);

% expected information gain (MI)
%--------------------------------------------------------------------------
E     = zeros(size(a,2),1);
for i = 1:numel(E)
    E(i) = spm_H(a(:,i));
end
E     = sum(a)*E/sum(a,'all');                 % H(Y|X)
E     = spm_H(sum(a,2)) - E;                   % M = IH(Y) - H(Y|X)

return


function I  = spm_H(a)

% entropy of a (conditioned) Dirichlet distribution 
%--------------------------------------------------------------------------
a   = a + 1/32;                                % preclude overflow
a0  = sum(a);
K   = numel(a);
I   = spm_betaln(a) + (a0 - K)*psi(a0) - sum((a - 1).*psi(a));
