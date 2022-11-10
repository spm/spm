function [DCM] = DEM_COVID_T
% FORMAT [DCM] = DEM_COVID_T
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This demo routine focuses on surveillance and mitigation strategies in
% the UK. It first estimate the parameters of a dynamic causal model for
% the epidemic in the United Kingdom. Crucially, in this inversion the data
% are supplemented with the total number of cases (in addition to positive
% cases and daily deaths). It rests upon an augmented DCM that includes a
% state of self isolation. Moving from a state of isolation depends upon a
% negative test. Tracing and tracking in this model is reflected as a small
% percentage of being tested if infected and asymptomatic. Otherwise, the
% baseline testing is in play. We will consider the effects of changing
% baseline testing and tracing and tracking at various phases of the
% outbreak.
%
% Finally, this routine performs a brief comparative analysis with Germany
% to see if the differences with the UK can be explained in terms of
% surveillance or clinical management.
% 
% NB: annotated notes appended to this routine illustrate a number of
% analyses simulations relevant to various containment, suppression and
% mitigation strategies. For example, the effect of early lockdowns, the
% effect of maintaining a tracking and tracing policy at the inception of
% the pandemic. In addition, there are notes showing how to incorporate
% serological data during inversion of the dynamic causal model.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Get data for the United Kingdom (including total tests R)
%==========================================================================
[Y,R] = DATA_COVID_UK('United Kingdom');
i     = min(size(Y,1),size(R,1));
Y     = [Y(1:i,:),R(1:i,:)];      % supplement cases with tests
R     = R/R(end);                 % empirical test rate

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% priors for this analysis (use size of population and estimate resistance)
%--------------------------------------------------------------------------
pE.N  = log(66);                 % population of UK (M)
pC.N  = 1/256;

% variational Laplace (estimating log evidence (F) and posteriors)
%==========================================================================

% complete model specification
%--------------------------------------------------------------------------
M.G   = @spm_COVID_gen;           % generative function
M.FS  = @(Y)sqrt(Y);              % feature selection  (link function)
M.pE  = pE;                       % prior expectations (parameters)
M.pC  = pC;                       % prior covariances  (parameters)
M.hE  = 0;                        % prior expectation  (log-precision)
M.hC  = 1/256;                    % prior covariances  (log-precision)
M.T   = size(Y,1);                % number of samples
M.R   = R;                        % empirical test rate
U     = [1 2 6];                  % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp] = spm_nlsi_GN(M,U,Y);

% assemble prior and posterior estimates (and log evidence)
%--------------------------------------------------------------------------
M.date = '25-01-2020';            % date of first time point
DCM.M  = M;
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.Y  = Y;

% save
%--------------------------------------------------------------------------
clear ans, save COVID_UK

% Predictions
%==========================================================================
M.T    = 256;                                    % six-month period
Y      = DCM.Y;                                  % empirical data
Ep     = DCM.Ep;                                 % posterior expectations
Cp     = DCM.Cp;                                 % posterior covariances

% show projections in terms of confidence intervals and superimpose data
%--------------------------------------------------------------------------
spm_figure('GetWin','UK'); clf;
%--------------------------------------------------------------------------
% (predicted outcomes). predictions with confidence intervals for deaths,
% positive cases and total counts
%--------------------------------------------------------------------------
spm_COVID_ci(Ep,Cp,Y,[1 2 6],M);

% add seasonal flu rates
%--------------------------------------------------------------------------
FLU = [1692,28330];            % death rate for seasonal flu (per season)
subplot(2,2,2), hold on
x   = get(gca,'XLim');
plot(x,[FLU(1) FLU(1)],'-.r',x,[FLU(2) FLU(2)],'-.r')
spm_axis tight

% supplemented with effective reproduction ratio (R) and herd immunity
%--------------------------------------------------------------------------
spm_figure('GetWin','R'); clf;
spm_COVID_ci(Ep,Cp,[],[4 5],M);


% and plot latent or hidden states
%--------------------------------------------------------------------------
spm_figure('GetWin','Predictions: UK'); clf;
%--------------------------------------------------------------------------
% (latent causes of observed consequences). The upper panels reproduce the
% expected trajectories of the previous figures, for an example country
% (here the United Kingdom). Here, the expected death rate is shown in
% blue, new cases in red, predicted recovery rate in orange. The black dots
% correspond to empirical data. The lower four panels show the evolution of
% latent (ensemble) dynamics, in terms of the expected probability of being
% in various states.
%--------------------------------------------------------------------------
M.T   = 365*1.5;                       % assume vaccination in 18 months
[Z,X] = spm_COVID_gen(Ep,M, [1 2]);
spm_COVID_plot(Z,X,DCM.Y,[],[1 2]);

% illustration of social distancing basal probability of leaving home
%--------------------------------------------------------------------------
spm_COVID_dashboard(DCM);


% Sensitivity analysis: which factors determine cumulative deaths?
%==========================================================================
spm_figure('GetWin','Sensitivity: UK'); clf
%--------------------------------------------------------------------------
% (sensitivity analysis). These panels show the change in outcome measures
% (i.e., death rate) as a function of time. The bar charts are the
% derivatives of final outcomes with respect to each testing parameter: 

% names{18} = 'track and trace'; %**
% names{19} = 'testing latency'; %**
% names{20} = 'test delay'; %**
% names{21} = 'test selectivity'; %**
% names{22} = 'sustained testing'; %**
% names{23} = 'baseline testing'; %**

% sensitivity analysis in terms of partial derivatives
%--------------------------------------------------------------------------
Np  = spm_length(Ep);
p   = [18 19 20 21 22 23];
V   = sparse(p,p,1,Np,Np); V = V(:,p);
dY  = spm_diff(@(P,M,U)spm_COVID_gen(P,M,U),Ep,M,1,1,{V});
t   = (1:size(dY,1))/7;
T   = 200;

subplot(2,1,1)
plot(t,dY,[T,T]/7,[min(dY(:)),max(dY(:))],'--')
xlabel('time (weeks)')
ylabel('First-order sensitivity'), box off
title('Effect on death rates','FontSize',16),
spm_axis tight, box off
legend(str.names{p}), legend boxoff

% cumulative effects over timea
%--------------------------------------------------------------------------
subplot(2,2,3), bar(sum(dY(1:T,:)))
xlabel('testing parameter')
ylabel('cumulative deaths')
title('Short-term effect','FontSize',16),
axis square, box off

subplot(2,2,4), bar(sum(dY(1:end,:)))
xlabel('testing parameter')
ylabel('cumulative deaths')
title('Long-term effect','FontSize',16),
axis square, box off
a   = axis; subplot(2,2,3), axis(a);


% Illustrate the effect of social distancing
%==========================================================================
spm_figure('GetWin','Testing: UK'); clf;
%--------------------------------------------------------------------------
% (the effects of surveillance). Here, trajectories are reproduced under
% different levels of tracking and tracing (in 32 steps). this parameter
% determines a probability of being tested if you are infected that
% asymptomatic. Crucially, it can be introduced before or after the first
% wave. Here, we look at the effect on the second wave

% testing parameters
%--------------------------------------------------------------------------
% P.ttt = 1/10000;              % test, track and trace
% P.ont = 2;                    % testing latency (months)
% P.del = 2;                    % test delay (days)
% P.tes = 1;                    % test selectivity (for infection)
% P.sus = 1/10000;              % sustained testing
% P.bas = 8/10000;              % baseline testing

% increase track and trace probability over 32 levels
%--------------------------------------------------------------------------
M.TTT = 20*7;                              % late start (20 weeks)
par   = linspace(1e-6,1,32);               % track, trace and test
S     = par;
R     = par;
for i = 1:numel(par)
    
    % social distancing threshold
    %----------------------------------------------------------------------
    P     = Ep;                            % expansion point
    P.ttt = log(par(i));                   % parameter
    [Y,X] = spm_COVID_gen(P,M,[1 2 9]);    % outcomes and causes
    S(i)  = sum(Y(:,1));                   % cumulative deaths
    R(i)  = max(Y(:,3))*par(i);            % infected and asymptomatic

    % plot results and hold graph
    %----------------------------------------------------------------------
    if rem(i,2)
        spm_COVID_plot(Y(:,1:2),X)
        for j = 1:6
            subplot(3,2,j), hold on
            set(gca,'ColorOrderIndex',1);
        end
    end
end

% cumulative deaths as a function of surveillance
%--------------------------------------------------------------------------
spm_figure('GetWin','Testing: UK II'); clf;

subplot(2,2,1), plot(par,S)
title('Cumulative deaths','FontSize',16),
xlabel('efficacy of track and trace')
ylabel('total deaths')
axis square, box off

subplot(2,2,2)
semilogy(par,R,par,0*par + 100000,':')
title('Peak identification','FontSize',16),
xlabel('efficacy of track and trace')
ylabel('identifiable cases')
axis square, box off
legend ({'peak identification','capacity'}), legend boxoff


% repeat with tracking and testing at the onset of the outbreak
%==========================================================================
M.TTT = 0;                       % onset
for i = 1:numel(par)
    
    % social distancing threshold
    %----------------------------------------------------------------------
    P     = Ep;
    P.ttt = log(par(i));
    Y     = spm_COVID_gen(P,M,[1,9]);
    S(i)  = sum(Y(:,1));
    R(i)  = max(Y(:,2))*par(i);

end

% cumulative deaths as a function of social distancing
%--------------------------------------------------------------------------
spm_figure('GetWin','Testing: UK II');
subplot(2,2,1), hold on, plot(par,S,'-.')
subplot(2,2,2), hold on, plot(par,R,'-.')
legend({'peak testing','capacity','peak testing (early)'}), legend boxoff



% Illustrate the effect at (quasi) endemic equilibrium
%==========================================================================
spm_figure('GetWin','Endemic equilibria'); clf;
%--------------------------------------------------------------------------
% This analysis focuses on the endpoint of the pandemic; namely, the
% endemic equilibrium, or at least that caused by an efficacious
% vaccination programme at 18 months. It considers increasing the
% sensitivity to infection, the selectivity for infected cases and track
% and trace, respectively.

% explore different testing strategies
%--------------------------------------------------------------------------
subplot(2,1,1), hold off
M.TTT = 20*7;                                   % late start
for i = 1:4
    
    % testing parameters
    %----------------------------------------------------------------------
    P    = Ep;                                  % expansion point
    if i == 2, P.bas = Ep.bas + 1;  end         % response
    if i == 3, P.tes = Ep.tes + 1;  end         % selective
    if i == 4, P.ttt = log(1/4);    end         % tracking
    
    % integrate and plot
    %----------------------------------------------------------------------
    y{i} = spm_COVID_gen(P,M,U);
    subplot(2,2,1), plot(y{i}(1:end,1)), hold on
    subplot(2,2,2), loglog(y{i}(32:end,1),y{i}(32:end,2)), hold on
    
end

% illustrate trajectories in terms of timeseries and in phase space
%--------------------------------------------------------------------------
subplot(2,2,1), title('Predicted deaths','FontSize',16)
hold on, plot(M.T,y{1}(end,1),'.','MarkerSize',32)
hold on, plot(M.T,y{4}(end,1),'.','MarkerSize',32)
xlabel('time (days)')
ylabel('daily new cases')
axis square, box off, spm_axis tight, hold off
legend({'current','enhanced testing','selective testing','tracing'}), legend boxoff

subplot(2,2,2), title('Endemic equilibria','FontSize',16)
hold on, plot(y{1}(end,1),y{1}(end,2),'.','MarkerSize',32)
hold on, plot(y{4}(end,1),y{4}(end,2),'.','MarkerSize',32)
xlabel('daily deaths')
ylabel('daily new cases')
axis square, box off, hold off

% cumulative deaths as a function of surveillance
%--------------------------------------------------------------------------
subplot(2,3,4), hold off
bar([y{1}(end,1), y{2}(end,1), y{3}(end,1), y{4}(end,1)])
title('Daily deaths','FontSize',16),
xlabel('testing strategy')
ylabel('death rate')
axis square, box off

subplot(2,3,5), hold off
bar([y{1}(end,3); y{2}(end,3); y{3}(end,3); y{4}(end,3)]/1e4)
title('Daily tests','FontSize',16),
xlabel('testing strategy')
ylabel('daily test (x 10,000)')
axis square, box off

subplot(2,3,6), hold off
S    = sum([y{1}(:,1), y{2}(:,1), y{3}(:,1), y{4}(:,1)]);
bar(max(S) - S)
title('Lives saved','FontSize',16),
xlabel('testing strategy')
ylabel('lives saved')
axis square, box off


% comparative analysis with Germany
%==========================================================================
% The final analysis considers the differences between the UK and Germany
% in terms of parameters that underwrite surveillance as opposed to
% clinical management. Using the same priors, the model is inverted for
% German data and the respective parameters are then compared with the UK.

% get data from Germany
%--------------------------------------------------------------------------
GCM.Y      = DATA_COVID_UK('Germany');

% set priors (in terms of population)
%--------------------------------------------------------------------------
GCM.M      = DCM.M;
GCM.M.pE.n = 4;
GCM.M.pE.N = log(83);

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
GCM.M.T = size(GCM.Y ,1);                % number of samples
GCM.M   = rmfield(GCM.M,'R');

U       = [1,2];
[Ep,Cp] = spm_nlsi_GN(GCM.M,U,GCM.Y);
GCM.Ep  = Ep;
GCM.Cp  = Cp;

% comparison of latent variables
%==========================================================================
spm_figure('GetWin','UK and Germany'); clf;
%--------------------------------------------------------------------------

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
GCM.M.T = 365;                            % one year
[Y,X]   = spm_COVID_gen(GCM.Ep,GCM.M,U);
spm_COVID_plot(Y,X,GCM.Y,[],U);
for j = 1:6
    subplot(3,2,j), hold on
    set(gca,'ColorOrderIndex',1);
end
DCM.M.T = 365;                            % one year
[Y,X] = spm_COVID_gen(DCM.Ep,DCM.M,U);
spm_COVID_plot(Y,X,DCM.Y,[],U);

% plot results
%--------------------------------------------------------------------------
spm_figure('GetWin','UK and Germany - parameters'); clf;
%--------------------------------------------------------------------------

% show 12 parameters with the greatest relative difference
%--------------------------------------------------------------------------
EP    = spm_vec(DCM.Ep);
GP    = spm_vec(GCM.Ep);
[d,q] = sort(abs(EP - GP),'descend');
q     = q(1:12);

param = fieldnames(Ep);
param = param(q);
EP    = EP(q);
GP    = GP(q);

subplot(2,2,1)
bar([EP GP])
set(gca,'XTick',1:numel(param),'Xticklabel',param,'FontSize',10)
set(gca,'XTickLabelRotation',90)
ylabel('log-parameters'), axis square, box off
title('Log parameters','FontSize',16)

subplot(2,2,2)
bar(exp([EP GP]))
set(gca,'XTick',1:numel(param),'Xticklabel',param,'FontSize',10)
set(gca,'XTickLabelRotation',90)
ylabel('scale-parameters'), axis square, box off
title('Scale parameters','FontSize',16)
legend({'UK','Germany'}), legend boxoff

subplot(2,2,3)
bar(GP - EP),                  hold on
plot([0 numel(q)],[0,0],'-.'), hold off
set(gca,'XTick',1:numel(param),'Xticklabel',param,'FontSize',10)
set(gca,'XTickLabelRotation',90)
ylabel('relative log-differences '), axis square, box off
title('Log differences','FontSize',16)

subplot(2,2,4)
bar(exp(GP - EP)),             hold on
plot([0 numel(q)],[1,1],'-.'), hold off
set(gca,'XTick',1:numel(param),'Xticklabel',param,'FontSize',10)
set(gca,'XTickLabelRotation',90)
ylabel('relative differences '), axis square, box off
title('Scale differences','FontSize',16)
set(gca,'YLim',[0 8])


% epilogue
%==========================================================================
% recent serologic testing in Germany can now be compared with the
% predictions of the above model inversion:
%--------------------------------------------------------------------------
spm_figure('GetWin','Heisberg'); clf;
%--------------------------------------------------------------------------
% A German antibody survey was the first out of the gate several weeks ago.
% At a press conference on 9 April, virologist Hendrik Streeck from the
% University of Bonn announced preliminary results from a town of about
% 12,500 in Heinsberg, a region in Germany that had been hit hard by
% COVID-19. He told reporters his team had found antibodies to the virus in
% 14% of the 500 people tested. By comparing that number with the recorded
% deaths in the town, the study suggested the virus kills only 0.37% of the
% people infected. (The rate for seasonal influenza is about 0.1%.) The
% team concluded in a two-page summary that "15% of the population can no
% longer be infected with SARS-CoV-2, and the process of reaching herd
% immunity is already underway." They recommended that politicians start to
% lift some of the regions' restrictions.
%--------------------------------------------------------------------------
T     = datenum('25-Jan-2020') - datenum('6-Apr-2020') + size(DCM.Y,1);
[Y,X] = spm_COVID_gen(GCM.Ep,GCM.M,5);
X     = X{2};

% get infection status predictions for this state
%--------------------------------------------------------------------------
j     = size(GCM.Y,1);
X     = X(1:j,2:4);

subplot(2,1,1)
plot(X*100)
ylabel('prevalence (%)'), title('Infection status','FontSize',16)
xlabel('time (days)'), set(gca,'YLim',[0 32])
axis square, box off

% superimpose empirical estimates
%--------------------------------------------------------------------------
t   = j - T; hold on
plot([t,t - 7],[1,1]*15.1,'Linewidth',8)
legend({'infected','infectious','immune','Streeck et al'})





% Study by the Office for National Statistics (ONS)
%==========================================================================
% this demonstration shows how to include serological data as point data
% during the time series. These kinds of data are included in a timeseries
% by setting the precision of the remaining time points to 0.
%
% The study, by the Office for National Statistics (ONS), tested 10,705
% people in more than 5,000 households and estimated 0.27% of the
% population in England were currently positive for Covid-19. The analysis
% suggests about 148,000 people across the entire population would have
% tested positive on any day between 27 April and 10 May 2020.
%
% As of 9 May 2020, the Office for National Statistics (ONS) had received
% the results of swab tests collected from 7,087 individual participants in
% the coronavirus (COVID-19) Infection Survey in England between 26 April
% and 8 May 2020.
%
% It is estimated that 0.24% of the population in England tested positive
% for COVID-19 (95% confidence interval: 0.14% to 0.40%).
% Almost one in five people in London - 17 per cent - have already had the
% coronavirus, according to surveillance testing, meaning that around
% 1.53million people have been infected with the virus and recovered
%
% After making adjustments for the accuracy of the assay and the age and
% gender distribution of the population, the overall adjusted prevalence in
% London increased from 1.5% in week 13 to 12.3% in weeks 15-16 and 17.5%
% in week 18.
%
% See: https://www.bmj.com/content/369/bmj.m1808
% https://www.gov.uk/government/publications/national-covid-19-surveillance-reports/sero-surveillance-of-covid-19
%--------------------------------------------------------------------------
load COVID_UK
%--------------------------------------------------------------------------
Y    = struct;
M    = DCM.M;
Y.y  = DCM.Y;
n    = size(Y.y,1);
i    = (1:n);
t    = i + datenum('25-Jan-2020');

% serological data
%--------------------------------------------------------------------------
ST   = [13 15 18]*7 + datenum('01-Jan-2020') - datenum('25-Jan-2020');
SY   = [1.5 12.3 17.5];

% add PCR and serological data
%--------------------------------------------------------------------------
PT   = (datenum('20-Apr-2020'):datenum('10-May-2020')) - datenum('25-Jan-2020');
Y.y  = [Y.y sparse(ST,1,SY,n,1) sparse(PT,1,0.24,n,1)];
U    = [1 2 6 5 8];

% precision components (with high precision for serology)
%--------------------------------------------------------------------------
Q    = cell(5,1);
Q{1} = speye(n,n);                         % death rate
Q{2} = speye(n,n);                         % positive cases
Q{3} = speye(n,n);                         % total tests
Q{4} = sparse(ST,ST,exp(8),n,n);           % antibodies
Q{5} = sparse(PT,PT,exp(4),n,n);           % prevalence
Q    = spm_cat(spm_diag(Q));

% add PCR and serological data
%--------------------------------------------------------------------------
Y.Q     = Q;
[Ep,Cp] = spm_nlsi_GN(M,U,Y); 


% graphics
%==========================================================================
spm_figure('GetWin','ONS: predictions'); clf;
%--------------------------------------------------------------------------
M.T     = 365;
[Z,X]   = spm_COVID_gen(Ep,M,[1 2]);
spm_COVID_plot(Z,X,DCM.Y(:,1:2));

% get infection status predictions for this state
%--------------------------------------------------------------------------
spm_figure('GetWin','ONS: serology'); clf;
%--------------------------------------------------------------------------
X     = X{2};
X     = X(1:M.T,2:4);
t     = (1:M.T) + datenum('25-Jan-2020');
subplot(2,1,1), hold off, plot(t,X*100), datetick('x')
ylabel('prevalence (%)'), title('Infection status','FontSize',16)
xlabel('time (date)'), set(gca,'YLim',[0 32])
axis square, box off

% superimpose empirical estimates
%--------------------------------------------------------------------------
hold on
T = PT + datenum('25-Jan-2020');
plot(T,ones(size(T))*0.24,'-b','Linewidth',2)

T = ST + datenum('25-Jan-2020');
plot(T ,SY,'or')
legend({'infected','infectious','immune','0.24%','17%'})
legend('boxoff')


% Independent SAGE: predictions for FTTIS
%==========================================================================
% these notes illustrate simulations of enhanced FTTI S from the date
% specified below. These plots are fatality rates and confidence intervals
% under 50% tracking and tracing efficacy.
%--------------------------------------------------------------------------
load COVID_UK, spm_figure('GetWin','Indie_SAGE: predictions'); clf
%--------------------------------------------------------------------------
M       = DCM.M;
Ep      = DCM.Ep;
Cp      = DCM.Cp;

subplot(2,1,1), hold on
M.T     = 380;
M.TTT   = datenum('01-Jul-2020') - datenum('25-Jan-2020');
Ep.ttt  = log(1/10000);
spm_COVID_ci(Ep,Cp,DCM.Y(:,1),1,M);
Y0      = spm_COVID_gen(Ep,M,1);

Ep.ttt  = log(1/2);
spm_COVID_ci(Ep,Cp,DCM.Y(:,1),1,M);
Y1      = spm_COVID_gen(Ep,M,1);

sprintf('lives saved %.0f',sum(Y0) - sum(Y1))


% Channel 4: predictions
%==========================================================================
% these notes illustrate simulations of early lockdown policies modelled in
% terms of decreasing the threshold (on the prevalence of infection) that
% elicits a social distancing response. This early intervention is modelled
% by reducing the social distancing thresholdduring the initial phase of
% the outbreak (as parameterised by M.T)
%--------------------------------------------------------------------------
load COVID_UK, spm_figure('GetWin','Channel 4: predictions'); clf
%--------------------------------------------------------------------------
M       = DCM.M;
Ep      = DCM.Ep;
Cp      = DCM.Cp;
M.T     = 365;
spm_COVID_ci(Ep,Cp,DCM.Y(:,1),1,M);

% and plot latent or hidden states
%--------------------------------------------------------------------------
spm_figure('GetWin','Channel 4: latent causes'); clf;
%--------------------------------------------------------------------------
[Z,X] = spm_COVID_gen(Ep,M,1);
spm_COVID_plot(Z,X,DCM.Y,[],1);

% hold plots for subsequent overlay and illustrate onset of lockdown
%--------------------------------------------------------------------------
for j = 1:6
    subplot(3,2,j), hold on
    set(gca,'ColorOrderIndex',1);
end
u     = 2/100;                              % threshold for lockdown
t     = (1:M.T) + datenum('25-Jan-2020');   % time
i     = find(X{1}(:,2) < u,1,'first');      % onset of lockdown
j     = sum(X{1}(:,2) < u);                 % duration of lockdown
subplot(3,2,3), plot([i,i]/7,[0 32],'r:')
text(i/7,32,datestr(t(i)))

% increase social distancing for the first 64 days
%--------------------------------------------------------------------------
M.TT  = 64;
[Z,X] = spm_COVID_gen(Ep,M,1);
spm_COVID_plot(Z,X,DCM.Y,[],1);

i     = find(X{1}(:,2) < u,1,'first');      % onset of lockdown
j     = sum(X{1}(:,2) < u);                 % duration of lockdown
subplot(3,2,3), plot([i,i]/7,[0 36],'r:')
text(i/7,36,datestr(t(i)))

% return to confidence interval graphics and overlay predictions
%--------------------------------------------------------------------------
spm_figure('GetWin','Channel 4: predictions')
subplot(2,1,1), hold on
spm_COVID_ci(Ep,Cp,[],1,M);



% tracking and testing
%==========================================================================
% repeat a similar kind of analysis in terms of the efficiency of FTTI
%--------------------------------------------------------------------------
spm_figure('GetWin','Channel 4: tracking and testing'); clf;
%--------------------------------------------------------------------------
subplot(2,1,1), hold on
M.T     = 365;
M.TT    = 0;
M.TTT   = datenum('01-Feb-2020') - datenum('25-Jan-2020');
Ep.ttt  = log(1/10000);
[Z,X]   = spm_COVID_gen(Ep,M,1);
spm_COVID_plot(Z,X,DCM.Y,[],1);

% hold the plots and repeat with an effective FTTI (of 50%)
%--------------------------------------------------------------------------
for j = 1:6
    subplot(3,2,j), hold on
    set(gca,'ColorOrderIndex',1);
end
Ep.ttt = log(1/2);
[Z,X]  = spm_COVID_gen(Ep,M,1);
spm_COVID_plot(Z,X,DCM.Y,[],1);


% tracking and testing - over a range of FTTI latencies
%--------------------------------------------------------------------------
spm_figure('GetWin','Channel 4: tracking and testing I'); clf;
%--------------------------------------------------------------------------
M      = DCM.M;
Ep     = DCM.Ep;
Cp     = DCM.Cp;
Ep.ttt = log(1/4);
M.TT   = 0;
M.T    = 180;

TTT   = linspace(0,128,8);
for i = 1:numel(TTT)
    M.TTT = TTT(i);
    [Z,X] = spm_COVID_gen(Ep,M,1);
    spm_COVID_plot(Z,X,DCM.Y,[],1);
    for j = 1:6
        subplot(3,2,j), hold on
        set(gca,'ColorOrderIndex',1);
    end
end

% tracking and testing - over a range of FTTI efficacies
%--------------------------------------------------------------------------
spm_figure('GetWin','Channel 4: tracking and testing II'); clf;
%--------------------------------------------------------------------------
M      = DCM.M;
Ep     = DCM.Ep;
Cp     = DCM.Cp;
M.TTT  = 0;
M.TT   = 0;
M.T    = 356*2;

TTT   = linspace(0,1,8);
for i = 1:numel(TTT)
    Ep.ttt = log(TTT(i));
    [Z,X]  = spm_COVID_gen(Ep,M,1);
    spm_COVID_plot(Z,X,DCM.Y,[],1);
    for j = 1:6
        subplot(3,2,j), hold on
        set(gca,'ColorOrderIndex',1);
    end
end


% risk to children returning to school
%==========================================================================
% these notes illustrate how to evaluate the risk of morbidity and
% mortality for children returning to a classroom of 15. This is based upon
% the prevalence of infection and the number of contacts. Risk here is
% assessed in relation to the probability of being involved in a fatal road
% traffic accident.
%
% The overall death rate from covid-19 has been estimated at 0.66%, rising
% sharply to 7.8% in people aged over 80 and declining to 0.0016% in
% children aged 9 and under.
%
% It reported that 0.04% of 10-19 year olds would probably require hospital
% care - as would 1.0% of people in their 20s, 3.4% of people aged 30-39,
% 4.3% aged 40-49, 8.2% aged 50-59, 11.8% in their 60s, 16.6% in their 70s,
% and 18.4% of those over 80.
%
% Verity R, Okell LC, Dorigatti I, et al. Estimates of the severity of
% coronavirus disease 2019: a model-based analysis. Lancet Infect Dis 2020
% Mar 30.
%--------------------------------------------------------------------------
load COVID_UK, spm_figure('GetWin','Risk'); clf;
%--------------------------------------------------------------------------

% 2018: fatalities, 1782 injuries, 160378
%--------------------------------------------------------------------------
rta    = 1782/66/365; disp(rta)
M      = DCM.M;
Ep     = DCM.Ep;
Cp     = DCM.Cp;
M.date = '25-Jan-2020';
M.T    = datenum('01-Sep-2020') - datenum(M.date);
M.TTT  = datenum('01-Jun-2020') - datenum(M.date);
%%% Ep.ttt = log(1/4);   % uncomment to implement tracking and tracing

% get confidence intervals for levels of infection
%--------------------------------------------------------------------------
[S,CS,Y,C] = spm_COVID_ci(Ep,Cp,[],7,M);
subplot(2,1,1)
set(gca,'XLim',[datenum('01-Apr-2020') datenum('01-Aug-2020')])
set(gca,'YLim',[0 24])
hold on

% plot an arbitrary threshold of 1% 
%--------------------------------------------------------------------------
d = (1:M.T) + datenum(M.date);
plot([datenum('01-Apr-2020') datenum('01-Aug-2020')],[1 1],'r-.')
i = find(Y > 1,1,'last');
plot([1,1]*d(i),[0 8],'r')
text(d(i),10,datestr(d(i)))
xlabel('date'), ylabel('percent risk')

% create a table of various risks (please see below)
%--------------------------------------------------------------------------
[Y,X]   = spm_COVID_gen(Ep,M,5);
P       = spm_vecfun(Ep,@exp);
datstr  = {'1-Jun-2020','15-Jun-2020','1-Sep-2020'};
dstr    = {'Jun1','Jun15','Sep1'};
R       = [15 P.Rin];

for n = 1:numel(R)
    for d = 1:numel(datstr)
        
        % prevalence of contagion
        %------------------------------------------------------------------
        T    = datenum(datstr{d}) - datenum('25-Jan-2020');
        p    = X{2}(T,3);
        
        % probability of contagious child in a class of 15
        %------------------------------------------------------------------
        Prev = 1 - (1 - p)^R(n);
        
        % probability of contracting virus
        %------------------------------------------------------------------
        Pinf = 1 - (1 - P.trn*p)^R(n);
        
        % probability of dying from virus
        %------------------------------------------------------------------
        Pfat   = Pinf*0.0016/100;
        
        Tab(1,n,d) = Prev*100;
        Tab(2,n,d) = Pinf*100;
        Tab(3,n,d) = Pfat*1e6;
    end
end

RowNames = {'contagious (%)',' contracting (%)','fatality (per million)'};
table(Tab(:,:,1),Tab(:,:,2),Tab(:,:,3),'RowNames',RowNames,'VariableNames',dstr)


return

% exemplar estimates of reproduction ratio based upon conventional
% modelling (these can be overlaid on the instantaneous estimates produced
% above)
% source: https://mrc-ide.github.io/covid19estimates/#/download
%--------------------------------------------------------------------------
% date: 13-02-2020; vs 25/1/2020 (19 days)

R   = [3.935988661
3.935988661
3.935988661
3.935988661
3.935988661
3.935988661
3.935965106
3.935957086
3.935947327
3.935935166
3.935919767
3.93590016
3.935875162
3.935843279
3.935802608
3.935750713
3.935684477
3.935599911
3.935491909
3.935353932
3.935177606
3.934952202
3.934663973
3.934295293
3.933823566
3.933219811
3.932446858
3.93145702
3.766451904
3.76496748
3.763071699
3.760660245
3.513081893
3.509724729
3.505488422
3.500177012
3.493538739
3.311507008
3.302189455
3.290611429
0.705844158
0.705043758
0.704098029
0.703073459
0.702027324
0.700995102
0.699995328
0.699036107
0.698119907
0.697246465
0.69641435
0.695621729
0.694866698
0.694147417
0.693462143
0.692809232
0.692187131
0.691594372
0.691029557
0.690491359
0.689978513
0.689489814
0.689024115
0.68858032
0.688157388
0.687754324
0.687370184
0.687004064
0.686655107
0.686322496
0.68600545
0.685703229
0.685415127
0.68514047
0.684878621
0.684628968
0.684390931
0.68416396
0.683947527
0.683741132
0.683544299
0.683356575
0.683177527
0.683006746
0.682843839];

% overlay Bayesian regression estimates
%--------------------------------------------------------------------------
spm_figure('GetWin','R');
subplot(2,2,1), hold on
t = (19 + (1:numel(R)))/7;
plot(t,R), spm_axis tight, set(gca,'YLim',[0 4])









