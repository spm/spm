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
% $Id: spm_MDP_G.m 6812 2016-06-18 11:16:21Z karl $


% get Bayesian surprise or mutual information
%==========================================================================

% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
p0   = exp(-16);
G    = 0;
qo   = 0;
qx   = 1;

% probability distribution over the hidden causes
%--------------------------------------------------------------------------
for f = 1:numel(x)
    qx = spm_cross(qx,x{f});
end

% accumulate expectation of entropy
%--------------------------------------------------------------------------
for i = 1:size(A{1},2)
    for j = 1:size(A{1},3)
        for k = 1:size(A{1},4)
            for l = 1:size(A{1},5)
                
                % probability over outcomes for this combination causes
                %----------------------------------------------------------
                po   = 1;
                for g = 1:numel(A)
                    po = spm_cross(po,A{g}(:,i,j,k,l));
                end
                po = po(:);
                qo = qo + qx(i,j,k,l)*po;
                G  = G  + qx(i,j,k,l)*po'*(log(po + p0));
                
            end
        end
    end
end

% subtract entropy of expectations
%--------------------------------------------------------------------------
G  = G - qo'*(log(qo + p0));
