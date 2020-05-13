function [DCM] = DEM_COVID_T
% FORMAT [DCM] = DEM_COVID_T
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% this demo routine focuses on surveillance and mitigation strategies in
% the UK. It first estimate the parameters of a dynamic causal model for
% the epidemic in the United Kingdom. Crucially, in this inversion the data
% are supplemented with the total number of cases (in addition to positive
% cases and daily deaths). It rests upon an augmented DCM that includes a
% state of self isolation. Moving from a state of isolation depends upon a
% negative test. Tracing and tracking in this model by a small percentage
% of being tested if infected that asymptomatic. Otherwise, the baseline
% testing is in play. We will consider the effects of changing baseline
% testing and tracing and tracking at various phases of the outbreak.
%
% Finally, this routine performs in brief comparative analysis with Germany
% to see if the differences between the UK can be explained in terms of
% surveillance  or clinical management.

%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: DEM_COVID_T.m 7849 2020-05-13 19:48:29Z karl $

% Get data for the United Kingdom
%==========================================================================
[Y,R] = DATA_COVID_UK('United Kingdom');
i     = min(size(Y,1),size(R,1));
Y     = [Y(1:i,:),R(1:i,:)];       % supplement cases with tests
R     = R/R(end);                  % empirical test rate

% assemble (Gaussian) priors over model parameters
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% priors for this analysis (use size of population and estimate resistance)
%--------------------------------------------------------------------------
pE.Tim = log(16);                 % period of immunity 
pE.N   = log(66);                 % population of UK (M)
pE.r   = log(1/3);                % non-susceptible proportion

pC.N   = 0;
pC.r   = 1/256;

% variational Laplace (estimating log evidence (F) and posteriors)
%==========================================================================

% complete model specification
%--------------------------------------------------------------------------
M.G   = @spm_COVID_gen;           % generative function
M.FS  = @(Y)sqrt(Y);              % feature selection  (link function)
M.pE  = pE;                       % prior expectations (parameters)
M.pC  = pC;                       % prior covariances  (parameters)
M.hE  = [0 0 -4];                 % prior expectation  (log-precision)
M.hC  = 1/64;                     % prior covariances  (log-precision)
M.T   = size(Y,1);                % number of samples
M.R   = R;                        % empirical test rate
U     = [1 2 6];                  % outputs to model

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
[Ep,Cp] = spm_nlsi_GN(M,U,Y);

% assemble prior and posterior estimates (and log evidence)
%--------------------------------------------------------------------------
DCM.M  = M;
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.Y  = Y;

% save
%--------------------------------------------------------------------------
clear ans, save COVID_UK

% Predictions
%==========================================================================
M.T    = 180;                                    % six-month period
Y      = DCM.Y;                                  % empirical data
Ep     = DCM.Ep;                                 % posterior expectations
Cp     = DCM.Cp;                                 % posterior covariances

% show projections in terms of confidence intervals and superimpose data
%--------------------------------------------------------------------------
spm_figure('GetWin','UK'); clf;
%--------------------------------------------------------------------------
% (predicted outcomes). This figure reports predicted new deaths and cases
% for an exemplar country; here, the United Kingdom. The panels on the left
% shows the predicted outcomes as a function of weeks. The blue line
% corresponds to the expected trajectory, while the shaded areas are 90%
% Bayesian credible intervals. The black dots represent empirical data,
% upon which the parameter estimates are based. The lower right panel shows
% the parameter estimates for the country in question. As in previous
% figures, the prior expectations are shown as in bars over the posterior
% expectations (and credible intervals). The upper right panel illustrates
% the equivalent expectations in terms of cumulative deaths.
%--------------------------------------------------------------------------
spm_COVID_ci(Ep,Cp,Y,[1 2 6],M)

% add seasonal flu rates
%--------------------------------------------------------------------------
FLU = [1692,28330];            % death rate for seasonal flu (per season)
subplot(2,2,2), hold on
x   = get(gca,'XLim');
plot(x,[FLU(1) FLU(1)],'-.r',x,[FLU(2) FLU(2)],'-.r')
spm_axis tight

% supplemented with effective reproduction ratio and herd immunity
%--------------------------------------------------------------------------
spm_figure('GetWin','R'); clf;
spm_COVID_ci(Ep,Cp,[],[4 5],M)


% and plot latent or hidden states
%--------------------------------------------------------------------------
spm_figure('GetWin','Predictions: UK'); clf;
%--------------------------------------------------------------------------
% (latent causes of observed consequences). The upper panels reproduce the
% expected trajectories of the previous figure, for an example country
% (here the United Kingdom). Here, the expected death rate is shown in
% blue, new cases in red, predicted recovery rate in orange. The black dots
% correspond to empirical data. The lower four panels show the evolution of
% latent (ensemble) dynamics, in terms of the expected probability of being
% in various states.
%--------------------------------------------------------------------------
M.T   = 365*1.5;                       % assume vaccination in 18 months
[Z,X] = spm_COVID_gen(Ep,M, [1 2]);
spm_COVID_plot(Z,X,DCM.Y,[],[1 2]);

% Sensitivity analysis: which factors determine cumulative deaths?
%==========================================================================
spm_figure('GetWin','Sensitivity: UK'); clf
%--------------------------------------------------------------------------
% (sensitivity analysis). These panels show the change in outcome measures
% (i.e., death rate) as a function of time. The bar charts are the
% derivatives of final outcomes with respect to each testing parameter: 

% names{18} = 'trace and test'; %**
% names{19} = 'response testing'; %**
% names{20} = 'test delay'; %**
% names{21} = 'test selectivity'; %**
% names{22} = 'sustained testing'; %**
% names{25} = 'baseline testing'; %**

% sensitivity analysis in terms of partial derivatives
%--------------------------------------------------------------------------
Np  = spm_length(Ep);
p   = [18 19 20 21 25];
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
% P.ttt = 1/10000;               % track, trace and test
% P.sen = 1/10000;               % response testing
% P.del = 2;                     % test delay (days)
% P.tes = 2;                     % test selectivity (for infection)
% P.exp = 1/10000;               % sustained testing
% 
% % immunity
% %--------------------------------------------------------------------------
% P.Tim = 32;                    % period of immunity (months)
% P.r   = 1e-6;                  % proportion resistant cases
% P.bas = 8/10000;               % baseline testing

% increase tracking test probability over 32 levels
%--------------------------------------------------------------------------
M.TTT = 20*7;                    % late start (20 weeks)
par   = linspace(1e-6,1,32);     % track, trace and test
S     = par;
R     = par;
for i = 1:numel(par)
    
    % social distancing threshold
    %----------------------------------------------------------------------
    P     = Ep;                  % expansion point
    P.ttt = log(par(i));         % parameter
    [Y,X] = spm_COVID_gen(P,M,U);
    S(i)  =  sum(Y(:,1));        % cumulative deaths
    R(i)  =  max(Y(:,3));        % peak test rates
    A(i)  = mean(Y(:,3));        % average test rates
    
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
xlabel('track and trace')
ylabel('total deaths')
axis square, box off

subplot(2,2,2)
plot(par,R,par,0*par + 250e3,':')
title('Peak testing rates','FontSize',16),
xlabel('track and trace')
ylabel('tests per day')
axis square, box off
legend ({'peak testing','capacity'}), legend boxoff


% repeat with tracking and testing at the onset of the outbreak
%==========================================================================
M.TTT = 0;                       % onset
for i = 1:numel(par)
    
    % social distancing threshold
    %----------------------------------------------------------------------
    P     = Ep;                  % expansion point
    P.ttt = log(par(i));         % parameter
    Y     = spm_COVID_gen(P,M,U);
    S(i)  =  sum(Y(:,1));        % cumulative deaths
    R(i)  =  max(Y(:,3));        % peak test rates
    A(i)  = mean(Y(:,3));        % average test rates

end

% cumulative deaths as a function of social distancing
%--------------------------------------------------------------------------
spm_figure('GetWin','Testing: UK II');
subplot(2,2,1), hold on, plot(par,S,'-.')
subplot(2,2,2), hold on, plot(par,R,'-.')
legend({'peak testing','capacity','peak (first wave)'}), legend boxoff



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
ylabel('cumulative deaths')
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


% epilogue
%==========================================================================
% recent serologic testing in Germany can now be compared with the
% predictions of the above model inversion:

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
% team concluded in a two-page summary that “15% of the population can no
% longer be infected with SARS-CoV-2, and the process of reaching herd
% immunity is already underway.” They recommended that politicians start to
% lift some of the regions’ restrictions.
%--------------------------------------------------------------------------
T     = 33;                    % days since Ab testing 6/3/2020 to date
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


return





% https://mrc-ide.github.io/covid19estimates/#/download

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









