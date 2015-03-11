function [qE,qC,Q] = spm_dcm_loo(DCM,X,field)
% Leave-one-out cross validation for empirical Bayes and DCM
% FORMAT [qE,qC,Q] = spm_dcm_loo(DCM,X,field)
%
% DCM   - {N [x M]} structure DCM array of (M) DCMs from (N) subjects
% -------------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% X     - second level design matrix, where X(:,1) = ones(N,1) [default]
% field - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%         'All' will invoke all fields
% 
% qE    - posterior predictive expectation (group effect)
% qC    - posterior predictive covariances (group effect)
% Q     - posterior probability over unique levels of X(:,2)
%__________________________________________________________________________
%
% This routine uses the posterior predictive density over the coefficients 
% of between subject effects encoded by a design matrix X. It is assumed
% that the second column of X contains classification predictor variables.
% A CV scheme is used to estimate the mixture of parameters at the first
% (within subject) level that are conserved over subjects in terms of a
% constant (first column of X) and differences (second column of X). Using
% a leave-one-out scheme, the predictive posterior density of the
% predictive variable is used to assess cross validation accuracy.
%
% For multiple models, this procedure is repeated for each model in the
% columns of the DCM array.
%
% See also: spm_dcm_peb.m and spm_dcm_ppd.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_loo.m 6373 2015-03-11 17:10:54Z karl $


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

% Repeat for each column if TEST is an array
%==========================================================================
if size(DCM,2) > 1
    
    % loop over models in each column
    %----------------------------------------------------------------------
    for i = 1:size(DCM,2)
        [p,q,r] = spm_dcm_loo(DCM(:,i),X,field);
        qE{i}   = p;
        qC{i}   = q;
        Q{i}    = r;
    end
    return
end

% Leave-one-out scheme
%==========================================================================
Ns    = numel(DCM);
for i = 1:Ns

    % get posterior predictive density for each subject
    %----------------------------------------------------------------------
    j         = 1:Ns;
    j(i)      = [];
    [Ep,Cp,P] = spm_dcm_ppd(DCM(i),DCM(j,1),X(i,:),X(j,:),field);
    qE(i)     = Ep;
    qC(i)     = Cp;
    Q(:,i)    = P;

end

% show results
%--------------------------------------------------------------------------
spm_figure('GetWin','LOO cross-validation');clf
subplot(2,2,1), spm_plot_ci(qE,qC), hold on
plot(X(:,2),':'), hold off
xlabel('subject'), ylabel('group effect')
title('Out of sample estimates','FontSize',16)
spm_axis tight, axis square

% classical inference on classification accuracy
%--------------------------------------------------------------------------
[T,df] = spm_ancova(X(:,1:2),[],qE(:),[0;1]);
r      = corr(qE(:),X(:,2));
p      = 1 - spm_Fcdf(T^2,df);
str    = sprintf('corr(df:%-2.0f) = %-0.2f: p = %-0.5f',df(2),r,p);

subplot(2,2,2)
plot(X(:,2),qE,'o','Markersize',8)
xlabel('group effect'), ylabel('estimate')
title(str,'FontSize',16)
spm_axis tight, set(gca,'XLim',[min(X(:,2))-1 max(X(:,2)+1)]), axis square


if size(Q,1) > 2
    subplot(2,1,2), imagesc(Q)
    xlabel('subject'), ylabel('levels of group effect')
    title('Posterior probability','FontSize',16)
    axis square
else
    subplot(2,1,2), bar(Q(2,:))
    xlabel('subject'), ylabel('posterior probability')
    title('Group effect','FontSize',16)
    axis([0 (Ns + 1) 0 1]), axis square
    axis square
end



