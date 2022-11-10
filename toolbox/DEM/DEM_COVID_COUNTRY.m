function DCM = DEM_COVID_COUNTRY(country)
% FORMAT DCM = DEM_COVID_COUNTRY(country)
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
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up and preliminaries
%==========================================================================
if nargin < 1, country = 'United Kingdom'; end

% get figure and data
%--------------------------------------------------------------------------
Fsi     = spm_figure('GetWin','SI'); clf;
data    = DATA_COVID_JHU(256);
c       = find(ismember({data.country},country));

% get and set priors
%--------------------------------------------------------------------------
[pE,pC] = spm_SARS_priors;
pE.N    = log(data(c).pop/1e6);
pC.N    = 0;

% augment priors with fluctuations
%--------------------------------------------------------------------------
onset   = datenum(data(c).date,'m/dd/yy');
M.date  = datestr(onset - 8,'dd-mm-yyyy');
i       = floor((datenum(date) - onset)/32);
j       = floor((datenum(date) - onset)/48);
k       = floor((datenum(date) - onset)/64);
pE.tra  = zeros(1,k);          % increases in transmission strength
pC.tra  = ones(1,k)/8;         % prior variance

% pE.pcr  = zeros(1,j);          % testing
% pC.pcr  = ones(1,j)/8;         % prior variance
% 
% pE.mob  = zeros(1,i);          % mobility
% pC.mob  = ones(1,i)/8;         % prior variance

% data for this country (here, and positive test rates)
%--------------------------------------------------------------------------
set(Fsi,'name',data(c).country)
Y       = [data(c).death, data(c).cases];
h       = log(sum(Y));
h       = mean(h) - h;
h(1)    = h(1) + 2;

% model specification
%==========================================================================
M.Nmax  = 256;                  % maximum number of iterations
M.G     = @spm_SARS_gen;        % generative function
M.FS    = @(Y)real(sqrt(Y));    % feature selection  (link function)
M.pE    = pE;                   % prior expectations (parameters)
M.pC    = pC;                   % prior covariances  (parameters)
M.hE    = h;                    % prior expectation  (log-precision)
M.hC    = 1/512;                % prior covariances  (log-precision)
M.T     = size(Y,1);            % number of samples
U       = [1 2];                % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp,Eh] = spm_nlsi_GN(M,U,Y);
DCM.M      = M;
DCM.Ep     = Ep;
DCM.Eh     = Eh;
DCM.Cp     = Cp;
DCM.Y      = Y;

% posterior predictions
%==========================================================================
spm_figure('GetWin',country); clf;
M.T     = datenum(date) + 256 - datenum(M.date,'dd-mm-yyyy');
[Z,X]   = spm_SARS_gen(DCM.Ep,M,[1 2 15]);
spm_SARS_plot(Z,X,DCM.Y,[1 2 15])

spm_figure('GetWin','confidence intervals'); clf;
% death rates
%--------------------------------------------------------------------------
subplot(3,1,1)
spm_SARS_ci(Ep,Cp,Y,1,M);hold on

% certified deaths
%--------------------------------------------------------------------------
spm_SARS_ci(Ep,Cp,[],15,M);

% notification rates
%--------------------------------------------------------------------------
subplot(3,1,2)
spm_SARS_ci(Ep,Cp,[],2,M);

% reproduction ratio
%--------------------------------------------------------------------------
subplot(3,1,3)
spm_SARS_ci(Ep,Cp,[],4,M);
plot(get(gca,'XLim'),[1,1],'-.r')
plot(datenum(date)*[1,1],get(gca,'YLim'),'-.b')

return

% repeat with rapid loss of immunity
%==========================================================================

% spm_figure('GetWin',[country ': 6 months']); clf;
% %------------------------------------------------------------------------
% Ep.Tim   = log(6);
% [Z,X]    = spm_SARS_gen(Ep,M,[1 2]);
% spm_SARS_plot(Z,X,Y)
% 
% spm_figure('GetWin','death rates'); hold on
% %------------------------------------------------------------------------
% spm_SARS_ci(Ep,Cp,[],1,M);
% Ep.Tim = DCM.Ep.Tim;

% repeat with efficient FTTIS
%==========================================================================
subplot(3,1,1); hold on
M.FTT = 1/4;
for c = [0 4]*7
    M.TTT = datenum(date) - datenum(M.date,'dd-mm-yyyy') + c;
    spm_SARS_ci(Ep,Cp,Y,1,M);
end

return

% repeat for Circuit break
%==========================================================================
spm_figure('GetWin','Circuit break'); clf
%--------------------------------------------------------------------------
CBT   = datenum(date) - datenum(DCM.M.date,'dd-mm-yyyy') + 4;
CBD   = 14;                         % duration of circuit breaker

M     = DCM.M;
M.T   = datenum('01-1-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
[Z,X] = spm_SARS_gen(DCM.Ep,M,[1 2 3]);
spm_SARS_plot(Z,X,DCM.Y)
for j = 1:6
    subplot(3,2,j), hold on
    set(gca,'ColorOrderIndex',1);
end

M.CBT = CBT;
M.CBD = CBD;
[Z,X] = spm_SARS_gen(DCM.Ep,M,[1 2 3]);
spm_SARS_plot(Z,X,DCM.Y)

% fatalities in confidence intervals
%--------------------------------------------------------------------------
spm_figure('GetWin','Circuit break - deaths'); clf
M     = DCM.M;
M.T   = datenum('01-01-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');
spm_SARS_ci(DCM.Ep,DCM.Cp,DCM.Y,1,M); hold on

M.CBT = CBT;
M.CBD = CBD;
spm_SARS_ci(DCM.Ep,DCM.Cp,[],1,M);


M.CBT = CBT;
M.CBD = CBD;
M.FTT = 1/4;
M.TTT = CBT + 28;
spm_SARS_ci(DCM.Ep,DCM.Cp,[],1,M);

M.CBT = CBT;
M.CBD = CBD;
M.FTT = .8;
M.TTT = CBT + 28;
spm_SARS_ci(DCM.Ep,DCM.Cp,[],1,M);

t   = datenum('01-1-2021','dd-mm-yyyy');
t0  = t - 7*24;
set(gca,'XLim',[t0,t])
set(gca,'YLim',[0,200])
set(gca,'XTick',[t0:14:t])
datetick('x','mmm-dd','keeplimits','keepticks')




