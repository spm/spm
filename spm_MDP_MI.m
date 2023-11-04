function [E,dEda,dEdA] = spm_MDP_MI(a,c)
% Expected information gain (i.e., mutual information)
% FORMAT [E,dEda,dEdA] = spm_MDP_MI(a,c)
%
% a    - Dirichlet parameters of a joint distribution
% c    - prior preferences
%
% E    - expected free energy (information gain minus cost)
% dEda - derivative with respect to Dirichlet parameters (a)
% dEdA - derivative with respect to joint density: A = a/sum(a(:))
%
% The mutual information here pertains to the expected distribution. See
% spm_dir_MI for the mutual information of a Dirichlet distribution per se
%
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% deal cells of (multimodal) tensors (omitting gradients)
%==========================================================================
if iscell(a)
    E     = 0;
    for g = 1:numel(a)
        if nargin > 1
            E = E + spm_MDP_MI(a{g},c{g});
        else
            E = E + spm_MDP_MI(a{g});
        end
    end
    return
end


% deal with tensors
%==========================================================================
a     = a(:,:);

% expected information gain (and negative cost)
%--------------------------------------------------------------------------
s     = sum(a,'all');
A     = a/s;
E     = spm_MI(A);

% expected (negative) cost
%--------------------------------------------------------------------------
if nargin > 1
    C = spm_log(c);
    E = E + C'*sum(A,2);
end

if nargout < 2, return, end

% dEdA = spm_log(A) - 1 - spm_log(sum(A,2))*ones(1,m) - ones(n,1)*spm_log(sum(A,1));
%--------------------------------------------------------------------------
dEdA   = spm_log(A./(sum(A,2)*sum(A,1))) - 1;

% expected (negative) cost
%--------------------------------------------------------------------------
if nargin > 1
    dEdA = plus(dEdA,C);
end

% dEda = dEdA.*dAda, dAda = (1/s - a/(s^2))
%--------------------------------------------------------------------------
dEda   = dEdA.*(1 - A)/s;


return


function I  = spm_MI(A)
% expected information gain of joint distribution
%--------------------------------------------------------------------------
I    =      A(:)'*spm_log(A(:)) - ...
        sum(A,1) *spm_log(sum(A,1)') - ...
        sum(A,2)'*spm_log(sum(A,2));

return
