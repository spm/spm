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

% retain eight countries with a transient epidemic
%--------------------------------------------------------------------------
[d,j] = sort(L,'ascend');
data  = data(j(1:8));

% Inversion (i.e., fitting) of empirical data
%==========================================================================
Fsi     = spm_figure('GetWin','SI'); clf;

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC] = spm_SARS_priors;

% Bayesian inversion (placing posteriors in a cell array of structures)
%--------------------------------------------------------------------------
Tim   = 32;
GCM   = cell(size(data(:)));
for i = 1:numel(data)
    
    pE.N  = log(data(i).pop/1e6);
    pC.N  = 0;
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
        M.G    = @spm_SARS_gen;        % generative function
        M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
        M.pE   = pE;                   % prior expectations (parameters)
        M.pC   = pC;                   % prior covariances  (parameters)
        M.hE   = 2;                    % prior expectation  (log-precision)
        M.hC   = 1/512;                % prior covariances  (log-precision)
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
        %------------------------------------------------------------------
        GCM{i,j} = DCM;
        F(i,j)   = Fp
    end
    
end

% save
%--------------------------------------------------------------------------
clear Fsi ans
save LANCET


% characterise second wave
%==========================================================================
load LANCET

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


% use UK
%--------------------------------------------------------------------------
j     = ismember({data.country},'United Kingdom');
DCM   = GCM(j,:);

% predictions for three periods of immunity
%--------------------------------------------------------------------------
spm_figure('GetWin','predictions'); clf;

I      = find(c > 0.5,1,'first');
M      = DCM{1}.M;
M.T    = 365*1.5;
M.date = '25-Jan-2020';
for  i = [I,numel(DCM)]
    spm_SARS_ci(DCM{i}.Ep,DCM{i}.Cp,DCM{i}.Y,1,M);
end
title(sprintf('Death rates (%.0f, and %.0f months)',Tim(I),Tim(end)))
datetick('x','mmm-yy')

% and plot confidence intervals around reproduction rate
%--------------------------------------------------------------------------
spm_figure('GetWin','reproduction ratio'); clf;
spm_SARS_ci(DCM{I}.Ep,DCM{I}.Cp,[],4,M);
subplot(2,1,1), hold on
plot([datenum('01-Feb-2020'), datenum('01-Aug-2021')],[1,1],'r-.')


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
[Z,X] = spm_SARS_gen(DCM{I}.Ep,M,[1 2]);
spm_SARS_plot(Z,X,DCM{I}.Y)


% illustrate accuracy
%==========================================================================
spm_figure('GetWin','data fits'); clf

% death rate
%--------------------------------------------------------------------------
N     = 8;
M.T   = 365;
for i = 1:N
    
    [Y,X] = spm_SARS_gen(GCM{i,I}.Ep,M,[1,2]);
    % spm_SARS_plot(Y,X,GCM{i,I}.Y)
 
    % death rate
    %----------------------------------------------------------------------
    subplot(2,1,1)
    d  = find(X{2}(:,2) > 1/1000,1,'first');
    t  = (1:size(Y,1)) - d;
    plot(t,cumsum(Y(:,1))), hold on
    t  = (1:size(GCM{i,I}.Y,1)) - d;
    plot(t,cumsum(GCM{i,I}.Y(:,1)),'k.')
    set(gca,'XLim',[0 (M.T + 32)])
    text(M.T - d,sum(Y(:,1)),data(i).country,'Fontweight','bold','FontSize',8)
    
    % new cases
    %----------------------------------------------------------------------
    subplot(2,1,2)
    t  = (1:size(Y,1)) - d;
    plot(t,cumsum(Y(:,2))), hold on
    t  = (1:size(GCM{i,I}.Y,1)) - d;
    plot(t,cumsum(GCM{i,I}.Y(:,2)),'k.')
    set(gca,'XLim',[0 (M.T + 32)])
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
global CHOLD
CHOLD = 0;
M.T   = 365*1.5;
for i = 1:N
    
    % posterior predictions of first and second peaks
    %----------------------------------------------------------------------
    [Y,X]    = spm_SARS_gen(GCM{i,I}.Ep,M,1);
    spm_SARS_plot(Y,X),drawnow
    subplot(3,2,1), hold on
    set(gca,'YLim',[0 1000])

    % table of first and second peaks
    %----------------------------------------------------------------------
    y         = GCM{i,I}.Y(:,1);
    m         = find(diff(Y(1:end - 1)) > 0 & diff(Y(2:end)) < 0);
    m         = round([m; numel(Y)]);
    m(1)      = find(y == max(y),1);
    dat(i,:)  = datenum(data(i).date) + m(1:2);
    tab{i,1}  = data(i).country;                 % country
    tab{i,2}  = datestr(dat(i,1));               % date first
    tab{i,3}  = datestr(dat(i,2));               % date second
    
    
    n         = 1 - X{1}(end,4);
    tab{i,4}  = ceil(max(X{2}(:,4))*100);       % population immunity
    tab{i,5}  = ceil(y(m(1)));                  % death rate (first)
    tab{i,6}  = ceil(Y(m(2)));                  % death rate (second)
    tab{i,7}  = ceil(exp(GCM{i,I}.Ep.res)*100); % non-infectious proportion
    tab{i,8}  = ceil(exp(GCM{i,I}.Ep.r)*100);   % non-susceptible proportion
    tab{i,9}  = ceil(exp(GCM{i,I}.Ep.N)*n);     % effective population
    tab{i,10} = ceil(data(i).pop/1e6);          % total population
    tab{i,11} = ceil(100*n);                    % effective proportion

end
subplot(3,2,1), legend({data.country}), legend 'boxoff'
CHOLD = 1;

% reorder table (mortality at first peak)
%--------------------------------------------------------------------------
vstr  = {'country','first','second'};
Tab   = cell2table(tab);
table(Tab(:,1),Tab(:,2),Tab(:,3),'VariableNames',vstr)

vstr  = {'country','herd','initial','secondary','noncon','nonsus','Effective','Total'};
table(Tab(:,1),Tab(:,4),Tab(:,5),Tab(:,6),Tab(:,7),Tab(:,8),Tab(:,9),Tab(:,10),'VariableNames',vstr)

% effective population as a percentage of total population
%--------------------------------------------------------------------------
% round(spm_vec(100*[tab{:,8}]./[tab{:,9}]))


% reproduction of the Lancet paper results
%==========================================================================
spm_figure('GetWin','LANCET'); clf
M.T   = 365;
for i = 1:N
    
    subplot(3,1,1), hold on
    t     = (1:M.T) + datenum(data(i).date,'m/dd/yy');
    T     = t(1:GCM{i,I}.M.T);
    pop   = data(i).pop/1e6;
    [Y,X] = spm_SARS_gen(GCM{i,I}.Ep,M,1);
    plot(t,cumsum(Y(:,1))/pop,T,cumsum(GCM{i,I}.Y(:,1))/pop,'.b')
    title('Cumulative deaths per million','FontSize',12)
    xlabel('date'),ylabel('per million')
    datetick('x','mmm')
    axis square, box off
    
    subplot(3,1,2), hold on
    T     = find(X{1}(:,2) < X{1}(8,2)/2,1,'first');
    N1    = sum(Y(1:T,1))/pop;
    N2    = sum(Y(T + (1:6*7),1))/pop;
    plot(log(N1),log(N2),'o')
    title('Cumulative deaths and lockdown','FontSize',12)
    xlabel('log deaths per million before')
    ylabel('log deaths per million after')
    axis square, box off
    
    subplot(3,1,3), hold on
    plot(X{2}(:,4)*100,cumsum(Y(:,1))/pop)
    title('Deaths and seroprevalence','FontSize',12)
    xlabel('seroprevalence (%)')
    ylabel('Cumulative deaths per million')
    axis square, box off

    
end
legend({data.country})



% model comparison
%==========================================================================
spm_figure('GetWin','LANCET BMR'); clf
G     = zeros(5,N);
for i = 1:N
    
    % populations proportions
    %----------------------------------------------------------------------
    Pop(1,i) = exp(GCM{i,I}.Ep.N);
    Pop(2,i) = exp(GCM{i,I}.Ep.m)*Pop(1,i);
    
    r(i)     = exp(GCM{i,I}.Ep.r);
    res(i)   = exp(GCM{i,I}.Ep.res);
    
    % full and reduced priors
    %----------------------------------------------------------------------
    pE     = GCM{i,I}.M.pE;
    pC     = GCM{i,I}.M.pC;
    qE     = GCM{i,I}.Ep;
    qC     = GCM{i,I}.Cp;
   
    % remove non-exposed population
    %----------------------------------------------------------------------
    rE     = pE;
    rE.m   = log(1/64);
    rC     = pC;
    F      = spm_log_evidence(qE,qC,pE,pC,rE,rC);
    G(1,i) = F;
    
    % remove susceptibility heterogeneity
    %----------------------------------------------------------------------
    rE     = pE;
    rE.r   = log(1/64);
    rC     = pC;
    F      = spm_log_evidence(qE,qC,pE,pC,rE,rC);
    G(2,i) = F;
    
    % remove transmission heterogeneity
    %----------------------------------------------------------------------
    rE     = pE;
    rE.res = log(1/64);
    rC     = pC;
    F      = spm_log_evidence(qE,qC,pE,pC,rE,rC);
    G(3,i) = F;
    
    % suppress lockdown
    %----------------------------------------------------------------------
    rE     = pE;
    rE.sde = log(1/2);
    rC     = pC;
    F      = spm_log_evidence(qE,qC,pE,pC,rE,rC);
    G(4,i) = F;
    
    % suppress herd immunity
    %----------------------------------------------------------------------
    rE     = pE;
    rE.Tim = log(1/32);
    rC     = pC;
    rC.Tim = 1/256;
    F      = spm_log_evidence(qE,qC,pE,pC,rE,rC);
    G(5,i) = F;
    
end

subplot(2,1,1), bar(G)
legend({data.country})
title('Bayesian model comparison','FontSize',12)
str = {'exposure','susceptibility','transmission','no lockdown','no immunity'};
set(gca,'XtickLabel',str)
xlabel('epidemiological model')
ylabel('log evidence')
box off, legend 'boxoff'

subplot(2,3,4), bar(Pop')
title('Population size','FontSize',12)
ylabel('millions')
xlabel('country')
axis square, box off, legend('total','effective'), legend 'boxoff'

subplot(2,3,5), bar(r*100)
title('non-susceptible','FontSize',12)
ylabel('proportion (%)')
xlabel('country')
axis square, box off

subplot(2,3,6), bar(res*100)
title('Non-contagious','FontSize',12)
ylabel('proportion (%)')
xlabel('country')
axis square, box off

