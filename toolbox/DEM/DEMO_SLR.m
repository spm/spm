function DEMO_SLR
% Demo of sparse logistic regression. This demonstration routines
% illustrates the use of sparse acoustic regression (as implemented in 
% spm_sparse_regression.m) to recover a small number of causes that best
% explain some dependent variable. For example, imagine we thought that a
% small number (e.g., four) of SNPs or copy number variants a generic study
% had sufficiently large effect sizes on some phenotypic measure to be
% worth pursuing. We might treat these (rare variants) as being expressed
% in the context of random effects (e.g., common variants and phenotypic
% measurement noise); however, we have many more potential causes than
% observations. This problem is addressed using (logistic) regression
% under sparsity constraints specified in terms of hyper priors over the
% precision (i.e. inverse variance) of model parameters. This provides the
% Bayesian shrinkage estimators of the regression coefficients that,
% crucially, can then be subject to Bayesian model reduction. Bayesian
% model reduction effectively eliminates redundant parameters that provided
% the optimal balance between accuracy and complexity.
%
% in the example below, we assume that we have 32 subjects with 128 
% independent variables (e.g., following some initial dimension reduction).
% The simulated data is generated with just four of the independent to see 
% whether these can be identified using sparse logistic regression and
% Bayesian model reduction.
%
% If the dependent variables are classes or probabilities a logistic
% transform is automatically applied. However, one can also use this
% routine for continuous (i.e., parametric phenotypes). the graphics
% produced by this demo report the results of sparse logistic regression
% using variational Laplace (i.e., approximate Bayesian inference and hyper
% priors). In addition, it reports the results and summary of the
% subsequent Bayesian model reduction.
%
% see also: spm_sparse_regression.m
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEMO_SLR.m 7454 2018-10-19 19:38:50Z karl $


% Genomic data (G)
%--------------------------------------------------------------------------
Ny   = 32;                                % number of observations
Nx   = 128;                               % number of explanatory variables

% Design matrix (U)
%--------------------------------------------------------------------------
X    = randn(Ny,Nx);
X0   = ones(Ny,1);

% Relevant explanatory variables (j)
%--------------------------------------------------------------------------
j    = 1:4;
P    = zeros(Nx,1);
P(j) = randn(length(j),1);

% Simulate (classification probability) data with a softmax
%--------------------------------------------------------------------------
y    = X*P + randn(Ny,1)/8;
y    = spm_softmax(y);


% Bayesian model reduction
%==========================================================================
RCM  = spm_sparse_regression(y,X,X0);

% show results
%--------------------------------------------------------------------------
spm_figure('Getwin','Model posterior (over families)'); clf

subplot(2,2,1)
spm_plot_ci(RCM.Ep(j),diag(RCM.Cp(j,j)));
xlabel('parameter')
title('Posterior density','FontSize',16)
axis square

subplot(2,2,2)
bar(P(j));
xlabel('parameter')
title('True (nonzero) parameters','FontSize',16)
axis square

subplot(2,2,3)
bar(RCM.Ep)
title('Recovered (sparse) profile','FontSize',16)
xlabel('parameter')
axis square, a = axis;

subplot(2,2,4)
bar(P)
title('True profile','FontSize',16)
xlabel('parameter')
axis square, axis(a)

