function [HP] = spm_dcm_peb_rnd_search(DCM,M,field)
% Re-randomisation testing for empirical Bayes and DCM
% FORMAT [HP] = spm_dcm_peb_rnd_search(DCM,M,field)
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
% HP     - maximum entropy hyperprior (precision 1/M.hC)
%__________________________________________________________________________
%
% This routine calls spm_dcm_peb_rnd to assess the distribution of log Bayes
% factors for different hyperpriors on between subject precision. It is
% assumed that the best hyperpriors maximise the entropy of the null
% distribution of ensuing p-values
%
% See also: spm_dcm_peb_rnd.m and spm_dcm_loo.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb_rnd_search.m 6474 2015-06-06 10:41:55Z karl $


% Set up
%==========================================================================
hP    = 4:4:32;
M.hE  = 0;
M.N   = 128;
bins  = 1/40:1/20:1;
for i = 1:numel(hP)
    rng('default');
    M.hC     = 1/hP(i);
    [p,P,f]  = spm_dcm_peb_rnd(DCM,M,field);
    E(i)     = p;
    F(i)     = 1 - f;
    H(:,i)   = hist(P,bins)';
end

% histogram of p-values and entropy (S)
%--------------------------------------------------------------------------
H     = H/M.N;
S     = -sum(H.*log(H + 1e-6));
[s,i] = max(S);

% repeat with maximum entropy hyperprior
%--------------------------------------------------------------------------
HP    = hP(i);
M.hC  = 1/hP(i);
spm_dcm_peb_rnd(DCM,M,field);

subplot(3,2,3)
imagesc(hP,bins,H),    hold on
plot(hP,H(end,:),'w'), hold on
plot([hP(1) hP(end)],[1 1]/20,':w'), hold off
xlabel('hyperprior'), ylabel('p-value')
title('Null distribution of p-values','FontSize',16)

subplot(3,2,5)
plot(hP,E,  'b',[hP(1) hP(end)],[1 1]/20,     ':b'), hold on
plot(hP,F,'-.b',[hP(1) hP(end)],[1 1]*log(20),':r'), hold on
plot(hP,S,  'r',[hP(i) hP(i)  ],[0 s],        ':r'), hold off
xlabel('hyperprior'), ylabel('p-value and entropy')
title('p-value and (null) entropy','FontSize',16)

