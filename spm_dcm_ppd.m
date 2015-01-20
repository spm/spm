function [qE,qC] = spm_dcm_ppd(TEST,TRAIN,X,field)
% Posterior predictive density for empirical Bayes and DCM
% FORMAT [qE,qC] = spm_dcm_ppd(TEST,TRAIN,X,field)
%
% TEST   - {1 [x M]} structure DCM array of new subject
% TRAIN  - {N [x M]} structure DCM array of (M) DCMs from (N) subjects
% --------------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% X      - second level design matrix, where X(:,1) = ones(N,1) [default]
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields
% 
% qE     - posterior predictive expectation
% qC     - posterior predictive covariances
%__________________________________________________________________________
%
% This routine inverts a hierarchical DCM using variational Laplace and
% Bayesian model reduction. In essence, it optimises the empirical priors
% over the parameters of a training set of first level DCMs, using 
% between subject constraints specified in the design matrix X. These
% optimised empirical priors are then used to parameterise a model of
% between subject effects for a single (test) subject. Usually, the second
% level of the design matrix specifies group differences and the posterior
% predictive density over this group effect can be used for classification
% or cross validation.
%
% See also: spm_dcm_peb.m and spm_dcm_loo.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_peb.m 6305 2015-01-17 12:40:51Z karl $


% Set up
%==========================================================================

% parameter fields
%--------------------------------------------------------------------------
if nargin < 4;
    field  = {'A','B'};
end
if strcmpi(field,'all');
    field = fieldnames(TEST(1,1).M.pE);
end

% Repeat for each column if TEST is an array
%==========================================================================
if size(TRAIN,2) > 1
    
    % loop over models in each column
    %----------------------------------------------------------------------
    for i = 1:size(TRAIN,2)
        [p,q]  = spm_dcm_ppd(TEST(1,i),TRAIN(:,i),X,field);
        qE{i}  = p;
        qC{i}  = q;
    end
    return
end

% Posterior predictive density
%==========================================================================

% evaluate empirical priors from training set
%--------------------------------------------------------------------------
PEB  = spm_dcm_peb(TRAIN,X,field);

% and estimate their contribution to the test subject
%--------------------------------------------------------------------------
peb  = spm_dcm_peb(TEST,PEB.Ep,field);
qE   = peb.Ep;
qC   = peb.Cp;
