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
cd('C:\Users\karl\Dropbox\Coronavirus')
close all, clear all

D      = DATA_WID_data;
[v,vi] = sort([D.vrate],'ascend');
[n,ni] = sort([D.population],'descend');
n      = 8;
ri     = unique([ni(1:n) vi(1:n)]);

%% free parameters of local model (fixed effects)
%==========================================================================
free  = {'n','o','m','exp','sde','qua','inn','lim','ons','rat','pcr','mob','pro','tra','rel'};

% (empirical) priors from the United Kingdom
%--------------------------------------------------------------------------
PCM   = load('DCM_UK.mat','DCM');
PCM   = PCM.DCM;
pE    = PCM.Ep;
pC    = spm_zeros(PCM.M.pC);
for i = 1:numel(free)
    pC.(free{i}) = PCM.M.pC.(free{i});
end

% remove age-specific vaccination rollout
%--------------------------------------------------------------------------
pE.rol = log([0.0001 (365 + 0)   32;
    0.001  (365 + 0)   32;
    0.02   (365 + 0)   32;
    0.03   (365 + 0)   32]);
pE.fol = log([exp(-16)   (365 + 512) 32;
    exp(-16) (365 + 128) 32;
    exp(-16) (365 + 64 ) 32;
    exp(-16) (365 + 32 ) 32]);
pC.rol = spm_zeros(pE.rol) + exp(-4);
pC.fol = spm_zeros(pE.fol);

PCM.Ep.rol = pE.rol;


r     = find(ismember({D.country},'Argentina'));
UK    = find(ismember({D.country},'United Kingdom'));
for r = ri
    
    fprintf('%d out of %d\n',r,numel(D));
    disp(D(r).country)
    try, clear Y; end
    
    % create data structure
    %----------------------------------------------------------------------
    Y(1).type = 'New cases';
    Y(1).unit = 'number/day';
    Y(1).U    = 2;
    Y(1).date = datenum(D(r).date);
    Y(1).Y    = D(r).cases;
    Y(1).h    = 0;
    
    Y(2).type = 'Daily deaths';
    Y(2).unit = 'number/day';
    Y(2).U    = 1;
    Y(2).date = datenum(D(r).date);
    Y(2).Y    = D(r).death;
    Y(2).h    = 2;
    
    Y(3).type = 'Positivity';
    Y(3).unit = 'percent';
    Y(3).U    = 23;
    Y(3).date = datenum(D(r).date);
    Y(3).Y    = D(r).positive * 100;
    Y(3).h    = 0;
    
    Y(4).type = 'New tests'; % daily PCR tests performed
    Y(4).unit = 'number/day';
    Y(4).U    = 6;
    Y(4).date = datenum(D(r).date);
    Y(4).Y    = D(r).tests;
    Y(4).h    = 0;
    
    Y(5).type = 'Percent vaccinated'; % percent vaccinated (England)
    Y(5).unit = 'percent';
    Y(5).U    = 22;
    Y(5).date = datenum(D(r).date);
    Y(5).Y    = 100 * D(r).vaccine / D(r).population;
    Y(5).h    = 0;
    
    Y(6).type = 'R-ratio'; % the production ratio
    Y(6).unit = 'ratio';
    Y(6).U    = 4;
    Y(6).date = datenum(D(r).date);
    Y(6).Y    = D(r).ratio;
    Y(6).h    = 0;
    
    % remove NANs, smooth and sort by date
    %----------------------------------------------------------------------
    [Y,S] = spm_COVID_Y(Y,'01-02-2020',16);
    
    % data structure with vectorised data and covariance components
    %----------------------------------------------------------------------
    xY.y  = spm_vec(Y.Y);
    xY.Q  = spm_Ce([Y.n]);
    hE    = spm_vec(Y.h);
    
    % country specific priors
    %======================================================================
    
    % population
    %----------------------------------------------------------------------
    N(4)    = 1e6* D(r).aged_70_older/D(r).population;
    s       = (1 - N(4) - 1/2)/(70 - D(r).median_age);
    N(2)    = (70 - 15)*s;
    N(3)    = N(2)*35/(35 + 20);
    N(2)    = N(2) - N(3);
    N(1)    = 1 - N(4) - N(3) - N(2);
    pE.N    = log(N(:)*D(r).population/1e6);
    
    % initial seeding
    %----------------------------------------------------------------------
    onset   = find(S(1:256,1) > max(S(1:256,1))/8,1) - 52;
    pE.n    = PCM.Ep.n - onset/8;
    pE.ons  = log(exp(PCM.Ep.ons) + onset);
    
    
    % testing rates
    %----------------------------------------------------------------------
    rel     = max(D(r).cases/D(r).population)/max(D(UK).cases/D(UK).population);
    pE.lim  = PCM.Ep.lim + log(rel);
    
    % vaccination rates
    %----------------------------------------------------------------------
    rel         = D(r).vrate;
    pE.rol(:,1) = PCM.Ep.rol(:,1) + log(rel);

    % hospital capacity
    %----------------------------------------------------------------------
    rel         = D(r).hospital_beds_per_thousand/D(UK).hospital_beds_per_thousand;
    if rel < 1
        hos     = erf(exp(PCM.Ep.hos))*rel;
        pE.hos  = log(erfinv(hos));
    elseif isnan(rel)
        pE.hos  = PCM.Ep.hos + log(1/512);
    end
    
    % initialise based on previous inversions
    %======================================================================
    try, load DCM_OWID DCM, end
    clear M
    M.P   = pE;
    try
        % use previous inversion
        %------------------------------------------------------------------
        field = fieldnames(DCM(r).Ep);
        for i = 1:numel(field)
            if all(size(pE.(field{i})) == size(DCM(r).Ep.(field{i})))
                M.P.(field{i}) = DCM(r).Ep.(field{i});
            end
        end
        disp('initialising with previous posteriors')
    end
    
    % model specification
    %======================================================================
    M.Nmax = 16;                   % maximum number of iterations
    M.date = '01-02-2020';         % start date
    M.G    = @spm_SARS_gen;        % generative function
    M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
    M.pE   = pE;                   % prior expectations (parameters)
    M.pC   = pC;                   % prior covariances  (parameters)
    M.hE   = hE + 2;               % prior expectation  (log-precision)
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
    DCM(r).S  = S;
    DCM(r).U  = U;
    DCM(r).D  = D(r);
    
    save('DCM_OWID.mat','DCM')
    
    % posterior predictions
    %======================================================================
    spm_figure('GetWin',D(r).country); clf;
    M.T     = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
    v       = [1 2 23];
    u       = [find(U == 1) find(U == 2) find(U == 23)];   
    [Z,X]   = spm_SARS_gen(DCM(r).Ep,M,v);
    spm_SARS_plot(Z,X,DCM(r).S(:,u),v)
    
    % data fits
    %----------------------------------------------------------------------
    spm_figure('GetWin',[D(r).country ' data fits']); clf;
    for i = 1:numel(U)
        subplot(4,2,i)
        spm_SARS_ci(Ep,Cp,S(:,i),U(i),M);
        title(Y(i).type,'FontSize',14), ylabel(Y(i).unit) 
    end
    
end

return



%% load
%==========================================================================
spm_figure('GetWin','Accuracy of fits'); clf;
load DCM_OWID

for r = 1:numel(DCM)
    i(r) = numel(DCM(r).F);
end
r     = find(i);

DCM   = DCM(r);
D     = DATA_WID_data;
D     = D(r);


%% check accuracy
%==========================================================================
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    Y   = DCM(r).M.T;
    
    M.T  = size(S,1);
    Z    = spm_SARS_gen(DCM(r).Ep,M,DCM(r).U);
    i    = isfinite(S(:,2));
    c    = corrcoef(S(i,2),Z(i,2));
    A(r) = c(1,2)^2;
    
end

subplot(2,2,1), bar(A)
subplot(2,2,2), histogram(A,8)

% data fits
%----------------------------------------------------------------------
[d,r] = max(A)
spm_figure('GetWin',[D(r).country ' best data fit']); clf;
[d,r] = max(A)
for i = 1:numel(U)
    subplot(4,2,i)
    spm_SARS_ci(DCM(r).Ep,DCM(r).Cp,DCM(r).S(:,i),DCM(r).U(i),M);
    title(DCM(r).M.T(i).type,'FontSize',14), ylabel(Y(i).unit)
end

% data fits
%----------------------------------------------------------------------
[d,r] = min(A)
spm_figure('GetWin',[D(r).country ' worst data fit']); clf;
for i = 1:numel(U)
    subplot(4,2,i)
    spm_SARS_ci(DCM(r).Ep,DCM(r).Cp,DCM(r).S(:,i),DCM(r).U(i),M);
    title(DCM(r).M.T(i).type,'FontSize',14), ylabel(Y(i).unit)
end



%% equable vaccination rates
%==========================================================================
spm_figure('GetWin','vaccination rates'); clf;

d      = datenum(date) - datenum(DCM(1).M.date,'dd-mm-yyyy') % duration of epidemi
pE.fol = log([0.0001   (d + 512) 32;
    0.001 (d + 128) 32;
    0.01 (d + 64 ) 32;
    0.01  (d + 32 ) 32]);

I     = find(ismember({D.country},'United Kingdom'));
Z     = 0;
E     = 0;
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');

    z   = spm_SARS_gen(Ep,M,22);
    
    if min(z) < 0, keyboard, end
    Z   = Z + z*D(r).population/100;
 
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    e   = spm_SARS_gen(Ep,M,22);
    E   = E + e*D(r).population/100;
    
    
    % graphics
    %----------------------------------------------------------------------
    subplot(3,1,1), plot(1:M.T,Z,1:M.T,E)
    xlabel('days'), ylabel('cumulative first doses')
    title('Cumulative first doses (global)','FontSize',14)
    legend({'current','equable'})
    
    subplot(3,1,2), plot(1:M.T,z,'k'), hold on
    title('Current rollout (country-specific)','FontSize',14)
    xlabel('days'), ylabel('cumulative first doses')
    i  = find(ismember(U,22));
    plot(S(:,i),'.k')
    
    subplot(3,1,3), plot(1:M.T,e,'k'), hold on
    title('Equable rollout (country-specific)','FontSize',14)
    xlabel('days'), ylabel('cumulative first doses')
    drawnow

end

%% scenario modelling - fatalities
%==========================================================================
spm_figure('GetWin','new deaths'); clf;
Z     = 0;
E     = 0;
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
        
    % certified deaths and 28 day death rates
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',1)
    
    % certified and 28 day deaths
    %----------------------------------------------------------------------
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],   15,M); hold on
              spm_SARS_ci(Ep,Cp,S(:,2),1,M);
    Z = Z + y;
    title('Certified and 28 day deaths','FontSize',14)
    
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    % certified deaths
    %----------------------------------------------------------------------
    subplot(3,1,2), set(gca,'ColorOrderIndex',1)
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],15,M);
    title('Certified deaths with equable rollout','FontSize',14)

    E  = E + y;
    
    % lives saved
    %----------------------------------------------------------------------
    subplot(3,1,3)
    plot(1:M.T,cumsum(Z),1:M.T,cumsum(E))
    xlabel('days'), ylabel('cumulative deaths')
    title('Total deaths (global)','FontSize',14)
    legend({'current','equable'})
    
end


%% scenario modelling - cases
%==========================================================================
spm_figure('GetWin','new cases'); clf;
Z     = 0;
E     = 0;
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
        
    % new cases
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',1)
    
    % certified and 28 day deaths
    %----------------------------------------------------------------------
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],2,M); hold on
    Z = Z + y;
    
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    % new cases
    %----------------------------------------------------------------------
    subplot(3,1,2), set(gca,'ColorOrderIndex',1)
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],2,M);
    title('with equable rollout','FontSize',14)

    E  = E + y;
    
    % cases avoided
    %----------------------------------------------------------------------
    subplot(3,1,3)
    plot(1:M.T,cumsum(Z),1:M.T,cumsum(E))
    xlabel('days'), ylabel('cumulative cases')
    title('cumulative cases (global)','FontSize',14)
    legend({'current','equable'})
    
end


%% scenario modelling - economy
%==========================================================================
spm_figure('GetWin','new deaths'); clf;
Z     = 0;
E     = 0;
for r = 1:numel(DCM)
    
    % unpack model and posterior expectations
    %----------------------------------------------------------------------
    M   = DCM(r).M;                                 % model (priors)
    Ep  = DCM(r).Ep;                                % posterior expectation
    Cp  = DCM(r).Cp;                                % posterior covariances
    S   = DCM(r).S;                                 % smooth timeseries
    U   = DCM(r).U;                                 % indices of outputs
    
    M.T = datenum(date) + 512 - datenum('01-02-2020','dd-mm-yyyy');
        
    % certified deaths and 28 day death rates
    %----------------------------------------------------------------------
    subplot(3,1,1), set(gca,'ColorOrderIndex',1)
    
    % mobility (transport)
    %----------------------------------------------------------------------
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],14,M); hold on
    Z = Z + y;
    title('Economic activity','FontSize',14)
    
    % with equable rollout
    %----------------------------------------------------------------------
    Ep.fol = pE.fol;
    
    % mobility
    %----------------------------------------------------------------------
    subplot(3,1,2), set(gca,'ColorOrderIndex',1)
    [~,~,y] = spm_SARS_ci(Ep,Cp,[],15,M);
    title('with equable rollout','FontSize',14)
    E  = E + y;
    
    % lives saved
    %----------------------------------------------------------------------
    subplot(3,1,3)
    plot(1:M.T,cumsum(Z),1:M.T,cumsum(E))
    xlabel('days'), ylabel('Cumulative economic cost')
    title('Total economic cost','FontSize',14)
    legend({'current','equable'})
    
end

return

