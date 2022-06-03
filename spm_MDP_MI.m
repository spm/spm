function [E,dEda,dEdA] = spm_MDP_MI(a,C)
% expected information gain (i.e., mutual information)
% FORMAT [I,dIda,dIdA] = spm_MDP_MI(a,C)
%
% a    - Dirichlet parameters of a joint distribution
% C    - log preferences
%
% E    - expected free energy (information gain minus cost)
% dEda - derivative with respect to Dirichlet parameters (a)
% dEdA - derivative with respect to joint density: A = a/sum(a(:));
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_KL_dir.m 8183 2021-11-04 15:25:19Z guillaume $

% deal with tensors
%--------------------------------------------------------------------------
a      = a(:,:);

% expected information gain
%--------------------------------------------------------------------------
s      = sum(a(:));
A      = a/s;
E      = spm_MI(A);

% expected (negative) cost
%--------------------------------------------------------------------------
if nargin > 1
    E = E + C*sum(A,2);
end

if nargout < 2, return, end

% dEdA = spm_log(A) - 1 - spm_log(sum(A,2))*ones(1,m) - ones(n,1)*spm_log(sum(A,1));
%--------------------------------------------------------------------------
dEdA   = spm_log(A./(sum(A,2)*sum(A,1))) - 1;

% dEda = dEdA/sum(a(:)) - sum(dEdA(:).*A(:))/(sum(a(:)));
%--------------------------------------------------------------------------
dEda   = (dEdA - sum(sum(dEdA.*A)))/s;

% expected (negative) cost: dCda = C/s - sum(C*sum(a,2))/(s^2)
%--------------------------------------------------------------------------
if nargin > 1
    dEdA = bsxfun(@plus,dEdA,C');
    dEda = bsxfun(@plus,dEda,(C - sum(C*sum(A,2)))'/s);
end

return


function I  = spm_MI(A)
% expected information gain of joint distribution
%--------------------------------------------------------------------------
I    =  A(:)'*spm_log(A(:)) - ...
        sum(A,1)*spm_log(sum(A,1)') - ...
        sum(A,2)'*spm_log(sum(A,2));

return
