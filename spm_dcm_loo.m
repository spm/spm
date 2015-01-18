function [qE,qC] = spm_dcm_loo(DCM,X,field)
% leave your now cross validation for empirical Bayes and DCM
% FORMAT [qE,qC] = spm_dcm_loo(DCM,X,field)
%
% DCM  - {N [x M]} structure DCM array of (M) DCMs from (N) subjects
% ------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% X      - second level design matrix, where X(:,1) = ones(N,1) [default]
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields
% 
% qE    - posterior predictive expectation (group effect)
% qC    - posterior predictive covariances (group effect)
%
%--------------------------------------------------------------------------
% This routine uses the posture predictive density over the coefficient of
% between subject effects encoded by the designer matrix X. It is assumed
% that the second column of X contains classification predictor variables.
% 81 and scheme is used to estimate the mixture of parameters at the first
% (within subject) level that are conserved over subjects in terms of a
% constant (first column of X) anddifferences (second column of X. Using a
% leave one out scheme, the predictive posterior density of the predictive
% variable is used to assess cross validation accuracy.
%
% For multiple models, this procedure is repeated for each model iin the
% columns of the DCM array.
%
% see also: spm_dcm_peb.m and spm_dcm_ppd.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb.m 6305 2015-01-17 12:40:51Z karl $

%  set up
%==========================================================================

% parameter fields
%--------------------------------------------------------------------------
if nargin < 4;
    field  = {'A','B'};
end
if strcmpi(field,'all');
    field = fieldnames(DCM(1,1).M.pE);
end

% repeat for each column if TEST is an array
%==========================================================================
if size(DCM,2) > 1
    
    % loop over models in each column
    %----------------------------------------------------------------------
    for i = 1:size(DCM,2)
        [p q] = spm_dcm_loo(DCM(:,i),X,field);
        qE{i} = p;
        qC{i} = q;
    end
    return
end

% leave one out scheme
%==========================================================================
Ns    = numel(DCM);
for i = 1:Ns

    % get posted predictive density for the subject
    %----------------------------------------------------------------------
    j       = 1:Ns;
    j(i)    = [];
    [Ep Cp] = spm_dcm_ppd(DCM(i),DCM(j,1),X(j,:),field);
    qE(i)   = Ep(2);
    qC(i)   = Cp(2,2);

end


% show results
%--------------------------------------------------------------------------
spm_plot_ci(qE,qC)
xlabel('subject')
ylabel('group effect')
title('Out of sample estimates','FontSize',16)
spm_axis tight






