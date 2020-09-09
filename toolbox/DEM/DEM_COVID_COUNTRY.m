function DCM = DEM_COVID_COUNTRY(country,T)
% FORMAT DCM = DEM_COVID_COUNTRY(country,T)
% country - country to model [default: 'United Kingdom')
% T       - prediction period (days)
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
% $Id: DEM_COVID_COUNTRY.m 7939 2020-09-09 11:02:14Z karl $

% set up and preliminaries
%==========================================================================
if nargin < 1, country = 'United Kingdom'; end
if nargin < 2, T       = 18*32;            end

% get figure and data
%--------------------------------------------------------------------------
Fsi     = spm_figure('GetWin','SI'); clf;
data    = DATA_COVID_JHU(168);
i       = find(ismember({data.country},country));

% get and set priors
%--------------------------------------------------------------------------
[pE,pC] = spm_SARS_priors;
pE.N    = log(data(i).pop/1e6);
pC.N    = 0;

% data for this country (here, and positive test rates)
%--------------------------------------------------------------------------
set(Fsi,'name',data(i).country)
Y       = [data(i).death, data(i).cases];

% model specification
%==========================================================================
M.date  = datestr(datenum(data(i).date,'m/dd/yy'),'dd-mm-yyyy');
M.G     = @spm_SARS_gen;        % generative function
M.FS    = @(Y)real(sqrt(Y));    % feature selection  (link function)
M.pE    = pE;                   % prior expectations (parameters)
M.pC    = pC;                   % prior covariances  (parameters)
M.hE    = 4;                    % prior expectation  (log-precision)
M.hC    = 1/512;                % prior covariances  (log-precision)
M.T     = size(Y,1);            % number of samples
U       = [1 2];                % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
Ep      = spm_nlsi_GN(M,U,Y);

% repeat with reduced likelihood precision
%--------------------------------------------------------------------------
M.P     = Ep;
M.hE    = 2;
[Ep,Cp] = spm_nlsi_GN(M,U,Y);

DCM.M   = M;
DCM.Ep  = Ep;
DCM.Cp  = Cp;
DCM.Y   = Y;

% posterior predictions
%==========================================================================
spm_figure('GetWin',country); clf;
M.T     = T;
[Z,X]   = spm_SARS_gen(Ep,M,[1 2]);
spm_SARS_plot(Z,X,Y)

spm_figure('GetWin','confidence intervals'); clf;
%--------------------------------------------------------------------------
spm_SARS_ci(Ep,Cp,Y,1,M);

spm_figure('GetWin','reproduction ratio'); clf
%--------------------------------------------------------------------------
spm_SARS_ci(Ep,Cp,[],4,M);

% repeat with rapid loss of immunity
%==========================================================================

% spm_figure('GetWin',[country ': 6 months']); clf;
% %------------------------------------------------------------------------
% Ep.Tim   = log(6);
% [Z,X]    = spm_SARS_gen(Ep,M,[1 2]);
% spm_SARS_plot(Z,X,Y)
% 
% spm_figure('GetWin','confidence intervals'); hold on
% %------------------------------------------------------------------------
% spm_SARS_ci(Ep,Cp,[],1,M);
% Ep.Tim = DCM.Ep.Tim;

% repeat with efficient FTTIS
%==========================================================================
spm_figure('GetWin','confidence intervals'); hold on
Ep.ttt = log(1/4);
M.TTT  = datenum(date) - datenum(M.date,'dd-mm-yyyy');
spm_SARS_ci(Ep,Cp,[],1,M);











