function DEMO_SLR
% demo of sparse logistic regression
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEMO_SLR.m 6306 2015-01-18 20:50:38Z karl $


% Genomic data (G)
%--------------------------------------------------------------------------
Ny   = 32;                                  % number of observations
Nx   = 128;                                 % number of explanatory variables

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
spm_plot_ci(RCM.Ep(j),RCM.Vp(j));
xlabel('parameter')
title('Posterior density','FontSize',16)
axis square

subplot(2,2,2)
bar(P(j));
xlabel('parameter')
title('True (nonzero) parameters','FontSize',16)
axis square

subplot(2,2,4)
bar(RCM.Ep)
title('True profile','FontSize',16)
xlabel('parameter')
axis square

subplot(2,2,3)
bar(P)
title('Recovered (sparse) profile','FontSize',16)
xlabel('parameter')
axis square
