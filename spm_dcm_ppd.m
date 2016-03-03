function [qE,qC,P] = spm_dcm_ppd(TEST,TRAIN,Y,X,field,iX)
% Posterior predictive density for empirical Bayes and DCM
% FORMAT [qE,qC,P] = spm_dcm_ppd(TEST,TRAIN,Y,X,field,i)
%
% TEST   - {1 [x M]} structure DCM array of new subject
% TRAIN  - {N [x M]} structure DCM array of (M) DCMs from (N) subjects
% --------------------------------------------------------------------
%     DCM{i}.M.pE - prior expectation of parameters
%     DCM{i}.M.pC - prior covariances of parameters
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% Y      - known values of design matrix for the test subject
% X      - second level design matrix, where X(:,1) = ones(N,1) [default]
% field  - parameter fields in DCM{i}.Ep to optimise [default: {'A','B'}]
%          'All' will invoke all fields (these constitute random effects)
% iX     - column of design matrix to be predicted   [default: iX=2]
% 
% qE     - posterior predictive expectation
% qC     - posterior predictive covariances
% P      - posterior probability over unique values of X(:,2)
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
% or cross validation. it is assumed that the unknown  explanatory
% variable in the design matrix  pertains to the second column unless
% otherwise specified
%
% See also: spm_dcm_peb.m and spm_dcm_loo.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ppd.m 6737 2016-03-03 12:05:51Z karl $


% Set up
%==========================================================================

% explanatory variable (between subject effect) of interest
%--------------------------------------------------------------------------
if nargin < 6;
    iX    = 2;
end

% parameter fields
%--------------------------------------------------------------------------
if nargin < 5;
    field = {'A','B'};
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
        [p,q,r] = spm_dcm_ppd(TEST(1,i),TRAIN(:,i),X,field);
        qE{i}   = p;
        qC{i}   = q;
        P{i}    = r;
    end
    return
end

% Posterior predictive density
%==========================================================================

% evaluate empirical priors from training set
%--------------------------------------------------------------------------
M.X   = X;
PEB   = spm_dcm_peb(TRAIN,M,field);

% and estimate their contribution to the test subject
%--------------------------------------------------------------------------
nX    = size(X,2);               % number of explanatory variables
bC    = var(X(:,iX))*4;
bC    = sparse(iX,1,bC,nX,1);    % prior covariances (variables)
M.X   = 1;                       % no between subject effects
M.W   = PEB.Ep;                  % emprical prior expectations
M.pC  = PEB.Ce;                  % emprical prior covariance (parameters)
M.bE  = Y; M.bE(:,iX) = 0;       % prior expectation (variables)
M.bC  = diag(bC + 1);            % prior covariances (variables)

for i = 1:4
    peb      = spm_dcm_peb(TEST,M,field);
    M.bE(iX) = peb.Ep(iX);
    M.bC     = diag(bC + exp(-i));
end

qE    = peb.Ep(iX);
qC    = peb.Cp(iX,iX);
pE    = peb.M.pE(iX);
pC    = peb.M.pC(iX,iX);

% Bayesian model reduction over levels of (second) explanatory variables
%--------------------------------------------------------------------------
x     = unique([X(:,iX); Y(:,iX)]);
for j = 1:length(x)
    F(j,1) = spm_log_evidence(qE,qC,pE,pC,x(j),0);
end
P     = spm_softmax(F);
