function DCM = DEM_COVID_WID(country)
% FORMAT DCM = DEM_COVID_WID(country)
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
% $Id: DEM_COVID_COUNTRY.m 8129 2021-08-02 18:08:36Z karl $

% set up and preliminaries
%==========================================================================

% get figure and data
%--------------------------------------------------------------------------
Fsi   = spm_figure('GetWin','SI'); clf;
D     = DATA_WID_data;

% triage data
%--------------------------------------------------------------------------
i     = zeros(1,numel(D));
for r = 1:numel(D)
    try
        d = D(r).vaccine / D(r).population;
        if max(d) > 1/1000 && isfinite(D(r).hospital_beds_per_thousand)
            i(r) = 1;
        end
    end 
end
D     = D(find(i));

% find(ismember({D.country},'United Kingdom'))

for r = 1:numel(D)
    
    fprintf('%d out of %d\n',r,numel(D));
    try, clear Y; end
    
    % create data structure
    %----------------------------------------------------------------------
    Y(1).type = 'new cases';
    Y(1).unit = 'number/day';
    Y(1).U    = 2;
    Y(1).date = datenum(D(r).date);
    Y(1).Y    = D(r).cases;
    Y(1).h    = 2;

    Y(2).type = 'Daily deaths';
    Y(2).unit = 'number/day';
    Y(2).U    = 1;
    Y(2).date = datenum(D(r).date);
    Y(2).Y    = D(r).death;
    Y(2).h    = 0;
    
    Y(3).type = 'Positivity';
    Y(3).unit = 'percent';
    Y(3).U    = 23;
    Y(3).date = datenum(D(r).date);
    Y(3).Y    = D(r).positive * 100;
    Y(3).h    = 0;
    
    Y(4).type = 'new tests'; % daily PCR tests performed
    Y(4).unit = 'number/day';
    Y(4).U    = 6;
    Y(4).date = datenum(D(r).date);
    Y(4).Y    = D(r).tests;
    Y(4).h    = 0;
    
    Y(5).type = 'number vaccinated'; % percent vaccinated (England)
    Y(5).unit = 'percent';
    Y(5).U    = 22;
    Y(5).date = datenum(D(r).date);
    Y(5).Y    = 100 * D(r).vaccine / D(r).population;
    Y(5).h    = 2;
    
    Y(6).type = 'R-ratio'; % the production ratio
    Y(6).unit = 'ratio';
    Y(6).U    = 4;
    Y(6).date = datenum(D(r).date);
    Y(6).Y    = D(r).ratio;
    Y(6).h    = 2;
    
    
    % remove NANs, smooth and sort by date
    %----------------------------------------------------------------------
    [Y,S] = spm_COVID_Y(Y,'01-02-2020',16);
    
    % data structure with vectorised data and covariance components
    %----------------------------------------------------------------------
    xY.y  = spm_vec(Y.Y);
    xY.Q  = spm_Ce([Y.n]);
    hE    = spm_vec(Y.h);
   
    % get and set priors
    %----------------------------------------------------------------------
    [pE,pC] = spm_SARS_priors;
    pE.N    = log(D(r).population/1e6);
    pC.N    = 0;
    
    % augment priors with fluctuations
    %----------------------------------------------------------------------
    k       = floor((datenum(date) - datenum('01-02-2020','dd-mm-yyyy'))/64);
    pE.tra  = zeros(1,k);          % increases in transmission strength
    pC.tra  = ones(1,k)/8;         % prior variance
    
    % country specific priors
    %----------------------------------------------------------------------
    rol       = log( 0.01 * max(D(r).vaccine / D(r).population) );
    pE.rol(1) = rol;
    pE.fol(1) = -16;
    pC.fol    = spm_zeros(pC.fol);
    pE.hos    = pE.hos + log(D(r).hospital_beds_per_thousand / 2.54);
    pE.vac    = log(512);
    pE.vef    = log(.2);
    
    % model specification
    %======================================================================
    M.Nmax = 32;                   % maximum number of iterations
    M.date = '01-02-2020';         % start date
    M.G    = @spm_SARS_gen;        % generative function
    M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
    M.pE   = pE;                   % prior expectations (parameters)
    M.pC   = pC;                   % prior covariances  (parameters)
    M.hE   = hE;                   % prior expectation  (log-precision)
    M.hC   = 1/512;                % prior covariances  (log-precision)
    M.T    = Y;                    % data structure
    
    U      = [Y.U];                % outputs to model
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);
    
    % save in DCM structure
    %----------------------------------------------------------------------
    DCM(r).M  = M;
    DCM(r).Ep = Ep;
    DCM(r).Eh = Eh;
    DCM(r).Cp = Cp;
    DCM(r).F  = F;
    DCM(r).Y  = S;
    DCM(r).U  = U;
    DCM(r).D  = D(r);

    
    % posterior predictions
    %======================================================================
    spm_figure('GetWin',D(r).country); clf;
    M.T     = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
    [Z,X]   = spm_SARS_gen(DCM(r).Ep,M,DCM(r).U);
    spm_SARS_plot(Z,X,DCM(r).Y,U)
    
  
    
end

return

spm_figure('GetWin','confidence intervals'); clf;

for r = 1:numel(D)
    
    
    
    % death rates
    %----------------------------------------------------------------------
    subplot(3,1,1)
    spm_SARS_ci(Ep,Cp,Y,1,M);hold on
    
    % certified deaths
    %----------------------------------------------------------------------
    spm_SARS_ci(Ep,Cp,[],15,M);
    
    % notification rates
    %----------------------------------------------------------------------
    subplot(3,1,2)
    spm_SARS_ci(Ep,Cp,[],2,M);
    
    % reproduction ratio
    %----------------------------------------------------------------------
    subplot(3,1,3)
    spm_SARS_ci(Ep,Cp,[],4,M);
    plot(get(gca,'XLim'),[1,1],'-.r')
    plot(datenum(date)*[1,1],get(gca,'YLim'),'-.b')
    
end
    
return

