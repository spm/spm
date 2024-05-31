function [E] = spm_dir_MI(a,c,h)
% Expected information gain (i.e., mutual information)
% FORMAT [E] = spm_dir_MI(a,c,h)
%
% a    - Dirichlet parameters of a joint distribution
% c    - prior preferences (outcomes)
% h    - prior preferences (states)
%
% E    - expected free energy (information gain minus cost)
%
% The mutual information here pertains to the Dirichlet distribution. See
% spm_MDP_MI for the mutual information of the expected categorical
% distribution.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% deal cells of (multimodal) tensors (omitting gradients)
%==========================================================================
if iscell(a)
    E     = 0;
    for g = 1:numel(a)
        if nargin > 2
            E = E + spm_dir_MI(a{g},c{g},h);
        elseif nargin > 1
            E = E + spm_dir_MI(a{g},c{g});
        else
            E = E + spm_dir_MI(a{g});
        end
    end
    return
end

% deal with tensors
%--------------------------------------------------------------------------
a     = a(:,:);

% mutual information: H(X) + H(Y) - H(X,Y);
%--------------------------------------------------------------------------
E     = spm_H(sum(a,2)) + spm_H(sum(a,1)) - spm_H(a(:));


% deal with costs
%==========================================================================
if nargin > 1
    A  = a/sum(a,'all');
end

% expected (negative) cost : outcomes
%--------------------------------------------------------------------------
if nargin > 1
    if numel(c)
        c = c(:)/sum(c,'all');
        C = spm_log(c);
        E = E + C'*sum(A,2);
    end
end

% expected (negative) cost : latent states
%--------------------------------------------------------------------------
if nargin > 2
    h = spm_cat(h(:));
    if numel(h)
        h = h(:)/sum(h,'all');
        H = spm_log(h);
        E = E + sum(A,1)*H;
    end
end

return


function I  = spm_H(a)

% differential entropy of a Dirichlet distribution 
%--------------------------------------------------------------------------
a0  = sum(a);
I   = psi(a0 + 1) - sum(a.*(psi(a + 1)))/a0;

return

% NB: continuous entropy of a  Dirichlet distribution 
%--------------------------------------------------------------------------
% a0  = sum(a);
% K   = numel(a);
% I   = spm_betaln(a) + (a0 - K)*psi(a0) - sum((a - 1).*psi(a));



