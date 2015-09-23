function [BMC,M] = spm_dcm_peb_rnd_search(DCM,M,field)
% Re-randomisation testing for empirical Bayes and DCM
% FORMAT [BMC,M] = spm_dcm_peb_rnd_search(DCM,M,field)
%
% DCM   - {N x 1} structure DCM array of (M) DCMs from (N) subjects
% -------------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% M.X    - second level design matrix, where X(:,1) = ones(N,1) [default]
% M.pE   - second level prior expectation of parameters
% M.pC   - second level prior covariances of parameters
% M.hE   - second level prior expectation of log precisions
% M.hC   - second level prior covariances of log precisions
% M.Q    - covariance components: {'single','fields','all','none'}
%
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields
% 
% BMC    - Bayesian model comparison structure 
% -------------------------------------------------------------
%     BMC.F    - free energy over joint model space
%     BMC.P    - posterior probability over models
%     BMC.Px   - posterior probability over 1st level models
%     BMC.Pw   - posterior probability over 2nd level models
%     BMC.M    - second level model
%     BMC.K    - model space
%__________________________________________________________________________
%
% This routine calls spm_dcm_peb_rnd to assess the distribution of log Bayes
% factors for different hyperpriors on between subject precision. It is
% assumed that the best hyperpriors maximise the entropy of the null
% distribution of ensuing p-values. This type of prior is then used to
% perform Bayesian model comparison. The optimised priors are in the second
% level model (M.hE, M.hC) in the output arguments.
%
% See also: spm_dcm_peb_rnd.m and spm_dcm_loo.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb_rnd_search.m 6561 2015-09-23 20:41:57Z karl $


% Set up
%==========================================================================
hE    = linspace(-4,2,8);
M.hC  = 1/32;
M.N   = 64;
bins  = 1/40:1/20:1;
for i = 1:numel(hE)
    rng('default');
    M.hE      = hE(i);
    [p,P,f,F] = spm_dcm_peb_rnd(DCM,M,field);
    CP(i)     = p;
    BP(i)     = 1 - f;
    H(:,i)    = hist(P,bins)';
    G(:,i)    = F';
end

% histogram of p-values and entropy (S)
%--------------------------------------------------------------------------
H     = H/M.N;
S     = -sum(H.*log(H + 1e-6));
[s,i] = max(S);

% repeat with maximum entropy hyperprior
%--------------------------------------------------------------------------
M.hE  = hE(i);
spm_dcm_peb_rnd(DCM,M,field);

subplot(3,2,3)
imagesc(hE,bins,H),    hold on
plot(hE,H(end,:),'w'), hold on
plot([hE(1) hE(end)],[1 1]/20,':w'), hold off
xlabel('hyperprior'), ylabel('p-value')
title('Null distribution of p-values','FontSize',16)

subplot(3,2,5)
plot(hE,CP,  'b',[hE(1) hE(end)],[1 1]/20,     ':b'), hold on
plot(hE,BP,'-.b',[hE(1) hE(end)],[1 1]*log(20),':r'), hold on
plot(hE,S,   'r',[hE(i) hE(i)  ],[0 s],        ':r'), hold off
xlabel('hyperprior'), ylabel('p-value and entropy')
title('p-values and entropy','FontSize',16)

% Bayesian model comparison (at first and second levels)
%--------------------------------------------------------------------------
BMC   = spm_dcm_bmc_peb(DCM,M,field);

