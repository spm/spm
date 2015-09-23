function [BMC,M] = spm_dcm_peb_test(DCM,M,field)
% BMC over first and second level models with classical hyperpriors
% FORMAT [BMC,M] = spm_dcm_peb_test(DCM,M,field)
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
% distribution of ensuing p-values. This hyperprior is then used to
% perform Bayesian model comparison. The optimised priors are in the second
% level model (M.hE, M.hC) in the output arguments.
%
% this (efficient) version simply tracks the base factor of an unlikely
% null model to find the prior expectations of between subject precision
% that renders the Bayes factorconsistent with a classical p-value (i.e.,
% resolves Lindley's paradox)
%
% See also: spm_dcm_bmc_peb.m and spm_dcm_loo.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb_test.m 6561 2015-09-23 20:41:57Z karl $


% Set up
%==========================================================================
rng('default');
M.hC  = 1/32;
M.N   = 64;

% find randomisation with the largest classical p-value
%--------------------------------------------------------------------------
[p,P,f,F,X] = spm_dcm_peb_rnd(DCM,M,field);


M0         = M;
M0.noplot  = 1;
[d,i]      = max(F);
M0.X(:,2)  = X(:,i);

% evaluate the Bayes factor over different hyperpriors
%--------------------------------------------------------------------------
hE    = linspace(-4,2,64);
for i = 1:numel(hE)
    M0.hE  = hE(i);
    bmc    = spm_dcm_bmc_peb(DCM,M0,field);
    G(i,:) = bmc.F;
end

j     = find(bmc.K(:,2));
P     = spm_softmax(G');
P     = sum(P(j,:),1);
G     = log(P./(1 - P));

% find hyperprior that is consistent with classical inference 
%--------------------------------------------------------------------------
p     = 1/M.N;
g     = log((1 - p)/p);
[d,i] = find(G > g,1);

% repeat with maximum entropy hyperprior
%--------------------------------------------------------------------------
if isempty(i)
    M.hE = hE(end);
else
    M.hE = hE(i);
end

spm_dcm_peb_rnd(DCM,M,field);

subplot(3,2,3)
plot(hE,G,'b',[hE(1) hE(end)],[g g],':b',[M.hE M.hE],[min(G) max(G)],'-.')
text(-2,g,sprintf('p < %-2.3f',p),'FontSize',10)
xlabel('Hyperprior'), ylabel('Log-Bayes factor')
title('Threshold under null','FontSize',16)


% Bayesian model comparison (at first and second levels)
%--------------------------------------------------------------------------
BMC   = spm_dcm_bmc_peb(DCM,M,field);

