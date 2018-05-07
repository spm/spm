function [G] = spm_MDP_G(A,x)
% auxiliary function for Bayesian suprise or mutual information
% FORMAT [G] = spm_MDP_G(A,x)
%
% A   - likelihood array (probability of outcomes given causes)
% x   - probability density of causes
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_G.m 7306 2018-05-07 13:42:02Z karl $


% get Bayesian surprise or mutual information
%==========================================================================

% preclude numerical overflow
%--------------------------------------------------------------------------
spm_log = @(x)log(x + exp(-16));

% probability distribution over the hidden causes: i.e., Q(x)
%--------------------------------------------------------------------------
qx    = spm_cross(x);

% accumulate expectation of entropy: i.e., E[lnP(o|x)]
%--------------------------------------------------------------------------
G     = 0;
qo    = 0;
for i = find(qx > exp(-16))'
    
    % probability over outcomes for this combination of causes
    %----------------------------------------------------------------------
    po   = 1;
    for g = 1:numel(A)
        po = spm_cross(po,A{g}(:,i));
    end
    po = po(:);
    qo = qo + qx(i)*po;
    G  = G  + qx(i)*po'*spm_log(po);
    
end

% subtract entropy of expectations: i.e., E[lnQ(o)]
%--------------------------------------------------------------------------
G  = G - qo'*spm_log(qo);


