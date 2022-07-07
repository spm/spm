function [F_samp,F_bound] = spm_BMS_F(alpha,lme,alpha0)
% Compute two lower bounds on model evidence p(y|r) for group BMS
% FORMAT [F_samp,F_bound] = spm_BMS_F(alpha,lme,alpha0)
% 
% alpha    - parameters of p(r|y)
% lme      - array of log model evidences 
%              rows:    subjects
%              columns: models (1..Nk)
% alpha0   - priors of p(r)
% 
% F_samp   - sampling estimate of <ln p(y_n|r)>
% F_bound  - lower bound on lower bound of <ln p(y_n|r)>
% 
% Reference:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. Neuroimage 2009 46(4):1004-17
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


alpha0 = sort(alpha0);
if alpha0(1) ~= alpha0(end)
    error('alpha0 should have identical values.')
end
alpha0 = alpha0(1);

a_sum    = sum(alpha);
psi_sum  = psi(a_sum);
psi_diff = psi(alpha) - psi_sum;
gm       = gammaln(alpha);

[s_samp,s_bound] = spm_BMS_F_smpl(alpha,lme,alpha0);

K = length(alpha);
F = 0;
for k = 1:K
    F = F - (alpha(k) - alpha0)*psi_diff(k) + gm(k);
end
F = F - gammaln(a_sum);

F_bound = F + s_bound;
F_samp  = F + s_samp;
