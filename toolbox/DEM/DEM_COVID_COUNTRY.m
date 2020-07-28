function DEM_COVID_COUNTRY(country)
% FORMAT DEM_COVID_COUNTRY(country)
% data    - data    to model [default: data = DATA_COVID_JHU]
% country - country to model [default: 'United Kingdom')
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine illustrates Bayesian model comparison using a line search
% over periods of imunity and pooling over countries. In brief,32 countries
% are inverted and 16 with the most informative posterior over the period
% of immunity are retained for Bayesian parameter averaging. The Christian
% predictive densities are then provided in various formats for the average
% country and (16) individual countries.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_COUNTRY.m 7909 2020-07-28 19:15:44Z karl $

% get figure and data
%--------------------------------------------------------------------------
Fsi     = spm_figure('GetWin','SI'); clf;
data    = DATA_COVID_JHU(168);
i       = find(ismember({data.country},country));

% get and set priors
%--------------------------------------------------------------------------
[pE,pC] = spm_COVID_priors;
pE.N    = log(data(i).pop/1e6);
pC.N    = 0;

% data for this country (here, and positive test rates)
%--------------------------------------------------------------------------
set(Fsi,'name',data(i).country)
Y       = [data(i).death, data(i).cases];

% model specification
%==========================================================================
M.G     = @spm_COVID_gen;       % generative function
M.FS    = @(Y)real(sqrt(Y));    % feature selection  (link function)
M.pE    = pE;                   % prior expectations (parameters)
M.pC    = pC;                   % prior covariances  (parameters)
M.hE    = [0 -2];               % prior expectation  (log-precision)
M.hC    = 1/512;                % prior covariances  (log-precision)
M.T     = size(Y,1);            % number of samples
U       = [1 2];                % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp] = spm_nlsi_GN(M,U,Y);

% posterior predictions
%==========================================================================
spm_figure('GetWin',country); clf;
M.T     = 385;
[Z,X]   = spm_COVID_gen(Ep,M,[1 2]);
spm_COVID_plot(Z,X,Y)

spm_figure('GetWin','confidence intervals'); clf;

[S,CS,Y,C] = spm_COVID_ci(Ep,Cp,Y,[1 2],M)




