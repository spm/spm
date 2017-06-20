function [qE,qC,Q] = spm_dcm_loo(DCM,M,field)
% Leave-one-out cross-validation for empirical Bayes and DCM
% FORMAT [qE,qC,Q] = spm_dcm_loo(DCM,M,field)
%
% DCM   - {N [x M]} structure DCM array of (M) DCMs from (N) subjects
% -------------------------------------------------------------------
%     DCM{i}.M.pE   - prior expectation of parameters
%     DCM{i}.M.pC   - prior covariances of parameters
%     DCM{i}.Ep     - posterior expectations
%     DCM{i}.Cp     - posterior covariance
%
% M.X       - second level design matrix, where X(:,1) = ones(N,1) [default]
% field     - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%             'All' will invoke all fields
% 
% qE        - posterior predictive expectation (group effect)
% qC        - posterior predictive covariances (group effect)
% Q         - posterior probability over unique levels of X(:,2)
% 
% This routine uses the posterior predictive density over the coefficients
% of between-subject effects encoded by a design matrix X. It is assumed
% that the second column of X contains classification or predictor
% variables. A cross-validation scheme is used to estimate the mixture of
% parameters at the first (within-subject) level that are conserved over
% subjects in terms of a constant (first column of X) and differences
% (second column of X). Using a leave-one-out scheme, the predictive
% posterior density of the predictive variable is used to assess
% cross-validation accuracy. For multiple models, this procedure is
% repeated for each model in the columns of the DCM array.
% 
% See also: spm_dcm_peb.m and spm_dcm_ppd.m
%__________________________________________________________________________
% Copyright (C) 2015-2017 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_loo.m 7120 2017-06-20 11:30:30Z spm $


% Set up
%==========================================================================

% parameter fields
%--------------------------------------------------------------------------
if nargin < 3
    try
        field = DCM{1}.field;
    catch
        field = {'A','B'};
    end
end
if strcmpi(field,'all')
    field = fieldnames(DCM{1}.M.pE);
end
if isnumeric(M)
    M = struct('X',M);
end

% Repeat for each column if TEST is an array
%==========================================================================
if size(DCM,2) > 1
    
    % loop over models in each column
    %----------------------------------------------------------------------
    for i = 1:size(DCM,2)
        [p,q,r] = spm_dcm_loo(DCM(:,i),M,field);
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
    [Ep,Cp,P] = spm_dcm_ppd(DCM(i),DCM(j,1),M.X(i,:),M.X(j,:),field);
    qE(i)     = Ep;
    qC(i)     = Cp;
    Q(:,i)    = P;

end

% Show results
%==========================================================================
if spm_get_defaults('cmdline'), return; end

spm_figure('GetWin','LOO cross-validation');clf
subplot(2,2,1), spm_plot_ci(qE,qC), hold on
plot(M.X(:,2),':'), hold off
xlabel('subject','FontSize',12), ylabel('group effect','FontSize',12)
title('Out of sample estimates','FontSize',16)
spm_axis tight, axis square

% classical inference on classification accuracy
%--------------------------------------------------------------------------
[T,df] = spm_ancova(M.X(:,1:2),[],qE(:),[0;1]);
r      = corrcoef(qE(:),M.X(:,2));
r      = full(r(1,2));

if isnan(T)
    p = NaN;
else
    p = 1 - spm_Tcdf(T,df(2));
end
str = sprintf('corr(df:%-2.0f) = %-0.2f: p = %-0.5f',df(2),r,p);

subplot(2,2,2)
plot(M.X(:,2),qE,'o','Markersize',8)
xlabel('group effect','FontSize',12), ylabel('estimate','FontSize',12)
title(str,'FontSize',16)
axis square

if size(Q,1) > 2
    % Continuous prediction
    subplot(2,1,2), imagesc(1:Ns,unique(M.X(:,2)),Q), hold on
    plot(M.X(:,2),'.c','MarkerSize',32), hold off
    xlabel('Subject','FontSize',12);
    ylabel('levels of group effect','FontSize',12)
    title('Predictive posterior (and true values)','FontSize',16)
    axis square xy
else
    % Binary classification
    subplot(4,1,3), bar(Q(1,:))
    xlabel('subject','FontSize',12);
    ylabel('posterior probability','FontSize',12);
    title('Group effect (group 1)','FontSize',16);
    axis([0 (Ns + 1) 0 1]);
    line([0 Ns+0.5],[0.95 0.95],'LineStyle','--','Color','r');
    
    subplot(4,1,4), bar(Q(2,:))
    xlabel('subject','FontSize',12);
    ylabel('posterior probability','FontSize',12);
    title('Group effect (group 2)','FontSize',16);
    axis([0 (Ns + 1) 0 1]);
    line([0 Ns+0.5],[0.95 0.95],'LineStyle','--','Color','r');
end
