function [F,Ep,Cp,pE,pC,Eh] = spm_COVID(Y,pE,pC,hC)
% Variational inversion of COVID model
% FORMAT [F,Ep,Cp,pE,pC,Eh] = spm_COVID(Y,pE,pC,hC)
% Y   - timeseries data
% pE  - prior expectation of parameters
% pC  - prior covariances of parameters
% hC  - prior covariances of precisions
% 
% F   - log evidence (negative variational free energy)
% Ep  - posterior expectation of parameters
% Cp  - posterior covariances of parameters
% pE  - prior expectation of parameters
% pC  - prior covariances of parameters
%
% This routine inverts a generative model of some timeseries data (Y),
% returning a variational (free energy) bound on log model evidence, and
% posterior densities of the model parameters (in terms of posterior
% expectations and covariances). This inversion uses standard variational
% Laplace; i.e., a (natural) gradient ascent on variational free energy
% under the Laplace assumption (i.e.,Gaussian priors and likelihood
% model).
%
% Model inversion entails specifying the generative model in terms of a log
% likelihood function and priors. These priors cover the model parameters
% and precision parameters that determine the likelihood of any given data.
% The precision priors (sometimes referred to as hyper priors) are
% specified in terms of the expectation and covariance of the log precision
% of random fluctuations about the predicted outcome variable. In this
% instance, the outcome variables are campus. This means that a square root
% transform allows a Gaussian approximation to the implicit (Poisson)
% likelihood distribution over observations.
% 
% The log likelihood function is provided as a subroutine in the (Matlab)
% code (spm_COVID_LL) below. However, because of Gaussian assumptions about
% the likelihood, we can use a simpler scheme, using the predicted outcomes
% from spm_COVID_gen, following a square root transform. The square root
% transform is treated as a feature selection or link function; please see
% the subroutine spm_COVID_FS.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% Gaussian priors over model parameters
%--------------------------------------------------------------------------
if nargin < 3
    [pE,pC] = spm_COVID_priors;
end
if nargin < 4
    hC = 1/256;
end

% complete model specification
%--------------------------------------------------------------------------
M.G   = @spm_COVID_gen;           % generative function
M.L   = @spm_COVID_LL;            % log-likelihood function
M.FS  = @spm_COVID_FS;            % feature selection (link function)
M.pE  = pE;                       % prior expectations (parameters)
M.pC  = pC;                       % prior covariances  (parameters)
M.hE  = 0;                        % prior expectation  (log-precision)
M.hC  = hC;                       % prior covariances  (log-precision)
M.T   =   size(Y,1);              % number of samples
U     = 1:size(Y,2);              % number of response variables

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);

% alternative inversion with Variational Laplace (curvature of log likelihood)
%--------------------------------------------------------------------------
% [Ep,Cp,F] = spm_nlsi_Newton(M,U,Y); Eh = 1;

return



% likelihood model
%__________________________________________________________________________

function L = spm_COVID_LL(P,M,U,Y)
% log-likelihood function
% FORMAT L = spm_COVID_LL(P,M,U,Y)
% P    - model parameters
% M    - model structure
% U    - inputs or control variables
% Y    - outputs or response variables
%
% This auxiliary function evaluates the log likelihood of a sequence of
% count data under Poisson assumptions
%__________________________________________________________________________

% generate prediction
%--------------------------------------------------------------------------
[T,N]  = size(Y);
M.T    = T;
y      = M.G(P,M,U);
y      = y(1:T,1:N);

% ensure all counts are greater than zero
%--------------------------------------------------------------------------
i      = logical(Y < 1); Y(i) = 1;
i      = logical(y < 1); y(i) = 1;

% place MDP in trial structure
%--------------------------------------------------------------------------
p      = spm_Npdf(y(:),Y(:),Y(:));
L      = sum(log(p + 1e-6));

% link function for feature selection (square root transform)
%__________________________________________________________________________

function Y = spm_COVID_FS(Y)
% feature selection for COVID model
% FORMAT Y = spm_COVID_FS(Y)
% P    - model parameters
% M    - model structurel
% U    - inputs or control variables
% Y    - outputs or response variables
%
% This auxiliary function takes appropriate gradients and performs a square
% root transform
%__________________________________________________________________________

% ensure all counts are greater than zero
%--------------------------------------------------------------------------
i      = logical(Y < 1); Y(i) = 1;

% square root transform
%--------------------------------------------------------------------------
Y      = sqrt(Y);

return
