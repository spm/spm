function [bor,F0,F1] = spm_BMS_bor(L,posterior,priors,C)
% Compute Bayes Omnibus Risk
% FORMAT [bor,F0,F1] = spm_BMS_bor(L,posterior,priors,C)
%
% L         Log model evidence table (models x  subjects)
% posterior .a model counts, .r model-subject probs
% priors    .a model counts
% C         if this field is specified then BOR under family prior 
%           is computed, otherwise BOR under model prior is computed.
%           C(k,f) = 1 if model k belongs to family f (0 otherwise)
%
% REFERENCES:
%
% Rigoux, L, Stephan, KE, Friston, KJ and Daunizeau, J. (2014)
% Bayesian model selection for group studies - Revisited. 
% NeuroImage 84:971-85. doi: 10.1016/j.neuroimage.2013.08.065
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_BMS_bor.m 6444 2015-05-21 11:15:48Z guillaume $


if nargin < 4
    options.families = 0;
    % Evidence of null (equal model freqs)
    F0 = FE_null(L,options); 
else
    options.families = 1;
    options.C = C;
    % Evidence of null (equal model freqs) under family prior
    [tmp,F0] = FE_null(L,options); 
end

% Evidence of alternative
F1 = FE(L,posterior,priors); 

% Implied by Eq 5 (see also p39) in Rigoux et al.
% See also, last equation in Appendix 2
bor = 1/(1+exp(F1-F0)); 


function [F,ELJ,Sqf,Sqm] = FE(L,posterior,priors)
% derives the free energy for the current approximate posterior
% This routine has been copied from the VBA_groupBMC function
% of the VBA toolbox http://code.google.com/p/mbb-vb-toolbox/ 
% and was written by Lionel Rigoux and J. Daunizeau
%
% See equation A.20 in Rigoux et al. (should be F1 on LHS)

[K,n] = size(L);
a0 = sum(posterior.a);
Elogr = psi(posterior.a) - psi(sum(posterior.a));
Sqf = sum(gammaln(posterior.a)) - gammaln(a0) - sum((posterior.a-1).*Elogr);
Sqm = 0;
for i=1:n
    Sqm = Sqm - sum(posterior.r(:,i).*log(posterior.r(:,i)+eps));
end
ELJ = gammaln(sum(priors.a)) - sum(gammaln(priors.a)) + sum((priors.a-1).*Elogr);
for i=1:n
    for k=1:K
        ELJ = ELJ + posterior.r(k,i).*(Elogr(k)+L(k,i));
    end
end
F = ELJ + Sqf + Sqm;


function [F0m,F0f] = FE_null (L,options)
% Free energy of the 'null' (H0: equal frequencies)
%
% F0m       Evidence for null (ie. equal probs) over models 
% F0f       Evidence for null (ie. equal probs) over families
%
% This routine derives from the VBA_groupBMC function
% of the VBA toolbox http://code.google.com/p/mbb-vb-toolbox/ 
% written by Lionel Rigoux and J. Daunizeau
%
% See Equation A.17 in Rigoux et al.

[K,n] = size(L);
if options.families
    f0 = options.C*sum(options.C,1)'.^-1/size(options.C,2);
    F0f = 0;
else
    F0f = [];
end
F0m = 0;
for i=1:n
    tmp = L(:,i) - max(L(:,i));
    g = exp(tmp)./sum(exp(tmp));
    for k=1:K
        F0m = F0m + g(k).*(L(k,i)-log(K)-log(g(k)+eps));
        if options.families
            F0f = F0f + g(k).*(L(k,i)-log(g(k)+eps)+log(f0(k)));
        end
    end
end
