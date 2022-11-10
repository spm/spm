function DEM_SARS
% FORMAT DEM_SARS
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


% Get data (see DATA_COVID): an array with a structure for each country
%==========================================================================
data  = DATA_COVID_JHU(16);
for i = 1:numel(data)
    p    = data(i).death;
    p    = p/sum(p);
    t    = (1:numel(p));
    t    = (t - mean(t)).^2;
    L(i) = t*p;
end

% retain ten countries with a transient epidemic
%--------------------------------------------------------------------------
N     = 10;
[d,j] = sort(L,'ascend');
data  = data(j(1:N));

% Inversion (i.e., fitting) of empirical data
%==========================================================================
Fsi     = spm_figure('GetWin','SI'); clf;

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC] = spm_SARS_priors;

% Bayesian inversion (placing posteriors in a cell array of structures)
%--------------------------------------------------------------------------
DCM   = cell(N,1);
for i = 1:numel(data)
    
    % initialisation (to previous posterior)
    %----------------------------------------------------------------------
    if i == 1
        M.P = pE;
    else
        M.P = Ep;
    end
    pE.N  = log(data(i).pop/1e6);
    pC.N  = 0;
    
    % data for this country (here, and positive test rates)
    %----------------------------------------------------------------------
    set(Fsi,'name',data(i).country)
    Y = [data(i).death, data(i).cases];
    
    % model specification
    %======================================================================
    M.G    = @spm_SARS_gen;       % generative function
    M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
    M.pE   = pE;                   % prior expectations (parameters)
    M.pC   = pC;                   % prior covariances  (parameters)
    M.hE   = 2;                    % prior expectation  (log-precision)
    M.hC   = 1/512;                % prior covariances  (log-precision)
    M.T    = size(Y,1);            % number of samples
    U      = [1 2];                % outputs to model
    
    % model inversion with Variational Laplace (Gauss Newton)
    %----------------------------------------------------------------------
    [Ep,Cp,Eh,Fp] = spm_nlsi_GN(M,U,Y);
        
    % assemble prior and posterior estimates (and log evidence)
    %----------------------------------------------------------------------
    GCM.M  = M;
    GCM.Ep = Ep;
    GCM.Cp = Cp;
    GCM.Eh = Eh;
    GCM.F  = Fp;
    GCM.Y  = Y;
    DCM{i} = GCM;
    
    % and plot latent or hidden states
    %----------------------------------------------------------------------
    spm_figure('GetWin',['Latent cases: ' data(i).country]); clf;
    %----------------------------------------------------------------------
    M.T   = 32*12;
    [Z,X] = spm_SARS_gen(DCM{i}.Ep,M,[1 2]);
    spm_SARS_plot(Z,X,DCM{i}.Y)
    
end

% save
%--------------------------------------------------------------------------
clear Fsi ans
save SARS
load SARS

% illustrate accuracy
%==========================================================================
spm_figure('GetWin','data fits'); clf

% death rate
%--------------------------------------------------------------------------
M.T   = 300;
for i = 1:N
    
    [Y,X] = spm_SARS_gen(DCM{i}.Ep,M,[1 2]);
    
    % death rate
    %----------------------------------------------------------------------
    subplot(2,1,1)
    d  = find(X{2}(:,2) > 1/1000,1,'first');
    t  = (1:size(Y,1)) - d;
    plot(t,cumsum(Y(:,1))), hold on
    t  = (1:size(DCM{i}.Y,1)) - d;
    plot(t,cumsum(DCM{i}.Y(:,1)),'k.')
    set(gca,'XLim',[0 (M.T + 32)])
    text(M.T - d,sum(Y(:,1)),data(i).country,'Fontweight','bold','FontSize',8)
    
    % new cases
    %----------------------------------------------------------------------
    subplot(2,1,2)
    t  = (1:size(Y,1)) - d;
    plot(t,cumsum(Y(:,2))), hold on
    t  = (1:size(DCM{i}.Y,1)) - d;
    plot(t,cumsum(DCM{i}.Y(:,2)),'k.')
    set(gca,'XLim',[0 (M.T + 32)])
    text(M.T - d,sum(Y(:,2)),data(i).country,'Fontweight','bold','FontSize',8)
    drawnow
    
end

subplot(2,1,1), xlabel('Time (days)'),ylabel('number')
title('Cumulative deaths','FontSize',16), box off
subplot(2,1,2), xlabel('Time (days)'),ylabel('number')
title('Cumulative cases','FontSize',16), box off


% Bayesian parameter averaging
%==========================================================================
BPA = spm_dcm_bpa(DCM,'nocd');
disp(spm_vecfun(BPA.Ep,@exp))

% Bayesian model averaging: empirical priors
%==========================================================================
[BMR,BMC] = spm_dcm_bmr_all(BPA,'All');



