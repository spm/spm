function DEM_COVID_I
% FORMAT DEM_COVID_I
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
% $Id: DEM_COVID_I.m 7870 2020-06-09 15:02:12Z karl $

% download data if required
%__________________________________________________________________________
if false
    url = 'https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/';
    urlwrite([url,'time_series_covid19_confirmed_global.csv'],'time_series_covid19_confirmed_global.csv');
    urlwrite([url,'time_series_covid19_deaths_global.csv'],   'time_series_covid19_deaths_global.csv');
    urlwrite([url,'time_series_covid19_recovered_global.csv'],'time_series_covid19_recovered_global.csv');
    
    url = 'https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/'
    urlwrite([url,'covid-19-tests-uk.csv'],'covid-19-tests-uk.csv');
end
%__________________________________________________________________________


% Get data (see DATA_COVID): an array with a structure for each country
%==========================================================================
data    = DATA_COVID_JHU(32);

% Inversion (i.e., fitting) of empirical data
%==========================================================================
Fsi     = spm_figure('GetWin','SI'); clf;

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC] = spm_COVID_priors;
pC.Tim  = 1/512;                       % tight shrinkage priors on immunity

% Bayesian inversion (placing posteriors in a cell array of structures)
%--------------------------------------------------------------------------
Tim   = 1:32;
GCM   = cell(size(data(:)));
for i = 1:numel(data)
    for j = 1:numel(Tim)
        
        % data for this country (here, and positive test rates)
        %------------------------------------------------------------------
        set(Fsi,'name',data(i).country)
        Y = [data(i).death, data(i).cases];
        
        % prior expectation of the period of immunity
        %------------------------------------------------------------------
        pE.Tim = log(Tim(j));
        
        % model specification
        %==================================================================
        M.G    = @spm_COVID_gen;       % generative function
        M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
        M.pE   = pE;                   % prior expectations (parameters)
        M.pC   = pC;                   % prior covariances  (parameters)
        M.hE   = 0;                    % prior expectation  (log-precision)
        M.hC   = 1/64;                 % prior covariances  (log-precision)
        M.T    = size(Y,1);            % number of samples
        U      = [1 2];                % outputs to model
        
        % initialisation (two previous posterior for this country)
        %------------------------------------------------------------------
        if j == 1
            M.P = pE;
        else
            M.P = Ep;
        end
        
        % model inversion with Variational Laplace (Gauss Newton)
        %------------------------------------------------------------------
        [Ep,Cp,Eh,Fp] = spm_nlsi_GN(M,U,Y);
        
        
        % assemble prior and posterior estimates (and log evidence)
        %------------------------------------------------------------------
        DCM.M  = M;
        DCM.Ep = Ep;
        DCM.Cp = Cp;
        DCM.Eh = Eh;
        DCM.F  = Fp;
        DCM.Y  = Y;
        
        % save this country (in a cell array)
        %----------------------------------------------------------------------
        GCM{i,j} = DCM;
        F(i,j)   = Fp
    end
    
end

% save
%--------------------------------------------------------------------------
clear Fsi ans
save COVID_I


% characterise second wave
%==========================================================================
load COVID_I

% remove China
%--------------------------------------------------------------------------
j     = ~ismember({data.country},'China');
F     = F(j,:);
GCM   = GCM(j,:);
data  = data(j);

% retain minimum entropy countries
%==========================================================================
for i = 1:size(F,1)
    p    = spm_softmax(spm_vec(F(i,:)));
    E(i) = Tim*p;
    H(i) = p'*log(p);
end

[G,J] = sort(H,'descend');
J     = J(1:16);
F     = F(J,:);
GCM   = GCM(J,:);
data  = data(J);

% posterior over a period of immunity
%==========================================================================
spm_figure('GetWin','period of immunity'); clf;

p   = spm_softmax(spm_vec(sum(F,1)));
c   = cumsum(p);
c0  = Tim(find(c < 0.05,1,'last'));
c1  = Tim(find(c > 0.95,1,'first'));
e   = Tim*p;

subplot(2,2,1)
semilogy(Tim,log(p),[1,1]*e,[min(log(p)), max(log(p))],'b-.')
xlabel('period of immunity (months)'),ylabel('log probability')
title('Log probability','FontSize',16),axis square, box off

subplot(2,2,2)
try, plot(Tim,c,[1,1]*c0,[0,1],'b-.',[1,1]*c1,[0,1],'b-.'), end
xlabel('period of immunity (months)'),ylabel('probability')
title('Cumulative probability','FontSize',16),axis square, box off


% Bayesian parameter average
%==========================================================================
DCM   = spm_dcm_bpa(GCM,'nocd');

% or use UK China
%--------------------------------------------------------------------------
j     = ismember({data.country},'United Kingdom');
DCM   = GCM(j,:);

% predictions for three periods of immunity
%--------------------------------------------------------------------------
spm_figure('GetWin','predictions'); clf;

I     = find(c > 0.5,1,'first');
for i = [2,I,numel(p)]
    M      = DCM{i}.M;
    M.T    = 365*1.5;
    M.date = '25-Jan-2020';
    spm_COVID_ci(DCM{i}.Ep,DCM{i}.Cp,[],1,M);
end
title(sprintf('Death rates (%.0f,%.0f, and %.0f months)',Tim(1),Tim(I),Tim(end)))
datetick('x','mmmyy')


% and plot latent or hidden states
%--------------------------------------------------------------------------
spm_figure('GetWin','latent causes'); clf;
%--------------------------------------------------------------------------
% (latent causes of observed consequences). The upper panels reproduce the
% expected trajectories of the previous figure. Here, the expected death rate is shown in
% blue, new cases in red, predicted recovery rate in orange and CCU
% occupancy in puce. The black dots correspond to empirical data. The lower
% four panels show the evolution of latent (ensemble) dynamics, in terms of
% the expected probability of being in various states.
%--------------------------------------------------------------------------
[Z,X] = spm_COVID_gen(DCM{I}.Ep,M,1:2);
spm_COVID_plot(Z,X)

% and plot confidence intervals around reproduction rate
%--------------------------------------------------------------------------
spm_figure('GetWin','reproduction rate'); clf;
spm_COVID_ci(DCM{I}.Ep,DCM{I}.Cp,[],4,M);
subplot(2,1,1), hold on
plot([datenum('01-Feb-2020'), datenum('01-Aug-2021')],[1,1],'r-.')


% illustrate accuracy
%==========================================================================
spm_figure('GetWin','data fits'); clf

% death rate
%--------------------------------------------------------------------------
M.T   = 180;
for i = 1:16
    
    [Y,X] = spm_COVID_gen(GCM{i,I}.Ep,M,1:2);
    % spm_COVID_plot(Y,X,GCM{i,I}.Y)
 
    % death rate
    %----------------------------------------------------------------------
    subplot(2,1,1)
    d  = find(X{2}(:,2) > 1/1000,1,'first');
    t  = (1:size(Y,1)) - d;
    plot(t,cumsum(Y(:,1))), hold on
    t  = (1:size(GCM{i,I}.Y,1)) - d;
    plot(t,cumsum(GCM{i,I}.Y(:,1)),'k.')
    set(gca,'XLim',[0 200])
    text(M.T - d,sum(Y(:,1)),data(i).country,'Fontweight','bold','FontSize',8)
    
    % new cases
    %----------------------------------------------------------------------
    subplot(2,1,2)
    t  = (1:size(Y,1)) - d;
    plot(t,cumsum(Y(:,2))), hold on
    t  = (1:size(GCM{i,I}.Y,1)) - d;
    plot(t,cumsum(GCM{i,I}.Y(:,2)),'k.')
    set(gca,'XLim',[0 200])
    text(M.T - d,sum(Y(:,2)),data(i).country,'Fontweight','bold','FontSize',8)
    drawnow
    
end

subplot(2,1,1), xlabel('Time (days)'),ylabel('number')
title('Cumulative deaths','FontSize',16), box off
subplot(2,1,2), xlabel('Time (days)'),ylabel('number')
title('Cumulative cases','FontSize',16), box off

% table of first and second peaks
%==========================================================================
spm_figure('GetWin','variability'); clf
M.T   = 365*2;
for i = 1:16
    
    % posterior predictions of first and second peaks
    %----------------------------------------------------------------------
    [Y,X] = spm_COVID_gen(GCM{i,I}.Ep,M,1);
    spm_COVID_plot(Y,X,GCM{i,I}.Y),drawnow
    subplot(3,2,1), hold on
    set(gca,'YLim',[0 1000])
    
    % table of first and second peaks
    %----------------------------------------------------------------------
    y        = GCM{i,I}.Y(:,1);
    m        = find(diff(Y(1:end - 1)) > 0 & diff(Y(2:end)) < 0);
    m        = round([m; numel(Y)]);
    m(1)     = find(y == max(y),1);
    dat(i,:) = datenum(data(i).date) + m(1:2);
    tab{i,1} = data(i).country;
    tab{i,2} = datestr(dat(i,1));
    tab{i,3} = datestr(dat(i,2));
    
    tab{i,4} = round(max(X{2}(:,4))*100);
    tab{i,5} = round(y(m(1)));
    tab{i,6} = round(Y(m(2)));
    tab{i,7} = round(exp(GCM{i,I}.Ep.res)*100);
    tab{i,8} = round(exp(GCM{i,I}.Ep.N));
    tab{i,9} = round(data(i).pop/1e6);
    
end

% reorder table (mortality at first peak)
%--------------------------------------------------------------------------
[d,i] = sort([tab{:,5}],'descend');
tab   = tab(i,:);
dat   = dat(i,:);

vstr  = {'country','first','second'};
Tab   = cell2table(tab);
table(Tab(:,1),Tab(:,2),Tab(:,3),'VariableNames',vstr)

vstr  = {'country','immunity','initial','secondary','resistant','Effective','Total'};
table(Tab(:,1),Tab(:,4),Tab(:,5),Tab(:,6),Tab(:,7),Tab(:,8),Tab(:,9),'VariableNames',vstr)

% effective population as a percentage of total population
%--------------------------------------------------------------------------
% round(spm_vec(100*[tab{:,8}]./[tab{:,9}]))



