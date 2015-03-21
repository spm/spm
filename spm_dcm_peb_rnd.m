function [p] = spm_dcm_peb_rnd(DCM,M,field)
% Re-randomisation testing for empirical Bayes and DCM
% FORMAT [p] = spm_dcm_peb_rnd(DCM,M,field)
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
% field - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%         'All' will invoke all fields
% 
% p      - classical (re-randomization) p - value
%__________________________________________________________________________
%
% This routine uses the posterior  density over the coefficients of
% between subject effects encoded by a design matrix X. It is assumed
% that the second column of X contains classification or predictor variables.
% The significance of group effects is assessed using re-randomization by
% permuting the element s(of the second) explanatory variable. This provides a
% null distribution for the relative free energy and a posterior
% probability over random permutations of the second level model.
%
% See also: spm_dcm_peb.m and spm_dcm_loo.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb_rnd.m 6385 2015-03-21 12:06:22Z karl $


% Set up
%==========================================================================

% parameter fields
%--------------------------------------------------------------------------
if nargin < 3;
    field  = {'A','B'};
end
if strcmpi(field,'all');
    field = fieldnames(DCM(1,1).M.pE);
end


% re-randomisation
%--------------------------------------------------------------------------
Ns  = size(M.X,1);
N   = 32;
M0  = M;
for i = 1:N
    M0.X(:,2) = M.X(randperm(Ns),2);
    bmc       = spm_dcm_bmc_peb(DCM,M0,field);
    F(i,:)    = bmc.F;
end

j   = 1 + size(bmc.K,1)/2;
F   = F(:,1) - F(:,j);

% Bayesian model comparison
%--------------------------------------------------------------------------
BMC = spm_dcm_bmc_peb(DCM,M,field);
G   = BMC.F(1) - BMC.F(j);

p   = (sum(F > G) + 1)/(N + 1);
r   = sort(F);
r   = r(fix((1 - 0.05)*N));

% show results
%--------------------------------------------------------------------------
spm_figure('GetWin','PEB-BMC');
subplot(3,2,1)
hist(F,32), hold on
plot([G G],[0 N/4],'--r'), hold on
plot([r r],[0 N/4],'--b'), hold on
text(G,N/6,sprintf('p < %-2.3f',p),   'FontSize',10), hold off
text(r,N/6,sprintf('p < %-2.3f',0.05),'FontSize',10), hold off
xlabel('Log Bayes Factor'), ylabel('Frequency')
title('Null distribution','FontSize',16)
axis square

subplot(3,2,2); title('Posterior','FontSize',16)
subplot(3,2,3); delete(gca)
subplot(3,2,5); delete(gca)

set(gcf,'Tag','PEB-RND','name','PEB_RND')
