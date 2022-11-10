function Tab = DEM_vaccination
% FORMAT Tab = DEM_vaccination
%
% Demonstration of COVID-19 modelling using variational Laplace
%__________________________________________________________________________
%
% This routine evaluates outcomes under some intervention over a specified
% set of dates. The outcomes are then tabulated and displayed in the MATLAB
% window. specify the duration and (parametric) nature of the intervention
% by editing the code below; namely, the non-pharmacological intervention
% structure NPI.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% sent plotting colours for overlay (with without interventions)
%--------------------------------------------------------------------------
global CHOLD
CHOLD = 1;

%% get empirical priors from fitting national data
%==========================================================================
if false
    DCM = DEM_COVID_UK;
    save('DCM_UK.mat','DCM')
else
    DCM = load('DCM_UK.mat','DCM');
    DCM = DCM.DCM;
end

%% parametric intervention (vaccination and contact tracing)
%==========================================================================
period = {DCM.M.date,'01-01-2022'};         % duration of epidemic
for i  = 1:6
    NPI(i).period = period;
    NPI(i).param  = {'vac','ttt'};  %%%% {'vac','ttt'}; 
    NPI(i).Q(1)   = 1e-8 + (i - 1)*1/5;
    NPI(i).Q(2)   = 0.24;
    NPI(i).dates  = {'08-12-2020',period{2}};
    str{i}        = sprintf('%g%s',round(100 * NPI(i).Q(1)),'%');
end

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','testing and cases'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = U([1 2 3 9]);

% quantify the effect of efficient intervention (NPI)
%==========================================================================
for i = 1:numel(NPI)
    
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X]  = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z(:,1:3),X,S,u(1:3))
    
    % record number of deaths
    %----------------------------------------------------------------------
    D(i) = sum(Z(:,2));
    C(i) = mean(Z(:,4));
    H(i) = 100 * (1 - X{2}(end,1));
    
end

% number of people vaccinated
%--------------------------------------------------------------------------
spm_figure('GetWin','vaccinations'); clf;
%--------------------------------------------------------------------------
Z     = spm_SARS_gen(Ep,M,23);
spm_SARS_plot(Z,[],[],23)
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    Z     = spm_SARS_gen(Ep,M,23,NPI(i));
    spm_SARS_plot(Z,[],[],23)
end
title('Vaccinations','FontSize',14)
ylabel('cumulative number')

% number of people dying
%--------------------------------------------------------------------------
subplot(3,2,3), bar(D), xlabel('vaccination rate'),ylabel('total deaths')
set(gca,'XTickLabels',str), title({'Cumulative','deaths'},'FontSize',14), box off

subplot(3,2,4), bar(C), xlabel('vaccination rate'),ylabel('percent')
set(gca,'XTickLabels',str), title({'Workplace','activity (%)'},'FontSize',14), box off

subplot(3,2,5), bar(H), xlabel('vaccination rate'),ylabel('percent')
set(gca,'XTickLabels',str), title({'Effective','herd immunity'},'FontSize',14), box off


% illustrate the effect of vaccination in terms of infection states
%--------------------------------------------------------------------------
i      = 4;                                 % for the i-th vaccination rate

[Z,X]  = spm_SARS_gen(Ep,M,1,NPI(1));
x(1,:) = X{2}(end,[1 4 5]);
[Z,X]  = spm_SARS_gen(Ep,M,1,NPI(i));
x(2,:) = X{2}(end,[1 4 5]);
x(3,:) = x(2,:) - x(1,:);
x      = 100 * x

subplot(3,2,6), bar(x),ylabel('percent'), xstr  = {str{1} str{i} 'diff.'};
set(gca,'XTickLabels',xstr), title({'Seroprevalence'},'FontSize',14), box off
legend({'susceptible','seropositive','seronegative'}), legend boxoff

return

%% quantify the effect of efficient intervention (NPI)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u     = [2 9 8];
nu    = numel(u);
Ylim  = [1200, 150, 4];
for j = 1:nu
    for i = 1:numel(NPI)
        subplot(nu,1,j), spm_SARS_ci(Ep,Cp,S(:,u(j)),U(u(j)),M,NPI(i)); hold on
        set(gca,'YLim',[0 Ylim(j)])
    end
end

return



%% scenario modelling for different kinds of lockdown
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-06-2021'};          % duration of epidemic
for i  = 1:3
    NPI(i).period = period;
    NPI(i).param  = {'sde'};
    NPI(i).Q      = exp(DCM.Ep.sde - (2 - i)/3 );
    NPI(i).dates  = {'01-01-2021',period{2}};
end

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','lockdowns'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 13;

% the effect of intervention (NPI) on latent states
%==========================================================================
for i = 1:numel(NPI)
    
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X]  = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
    
end

% illustrate outcomes
%--------------------------------------------------------------------------
spm_figure('GetWin','hospital cases'); clf;
%--------------------------------------------------------------------------
for i = 1:numel(NPI)
    
    subplot(3,1,1), hold on
    u      = 13;
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    
    subplot(3,1,2), hold on
    u      = 1;
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    
    subplot(3,1,3), hold on
    u      = 27;
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    
end

return


%% scenario modelling for different Lockdown timings
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention(relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'28-02-2021'};         % duration of epidemic
for i  = 1:2
    NPI(i).period = period;
    NPI(i).param  = {'sde'};
    NPI(i).Q      = exp(DCM.Ep.sde - (i - 1)/3 );
    NPI(i).dates  = {'22-12-2020','28-02-2021',};
end

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','lockdowns'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 13;

% the effect of intervention (NPI) on latent states
%==========================================================================
for i = 1:numel(NPI)
    
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X]  = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
    
end

%% illustrate outcomes
%--------------------------------------------------------------------------
spm_figure('GetWin','hospital cases'); clf;
%--------------------------------------------------------------------------
for i = 1:numel(NPI)
    
    subplot(3,1,1), hold on
    u      = 13;
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    
    subplot(3,1,2), hold on
    u      = 1;
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    
    subplot(3,1,3), hold on
    u      = 27;
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    
end


%% scenario modelling for lockdown and vaccination in response to the 
% Prime Minister's statement on April 13, 2021
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-10-2021'};         % duration of epidemic

NPI(1).period = period;
NPI(1).param  = {'qua'};
NPI(1).Q      = 1/32;
NPI(1).dates  = {'01-02-2021','01-10-2021',};

NPI(2).period = period;
NPI(2).param  = {'rol'};
NPI(2).Q      = exp(-16);
NPI(2).dates  = {'01-02-2021','01-10-2021',};

NPI(3).period = period;
NPI(3).param  = {'rol','qua'};
NPI(3).Q      = [exp(-16) 1/32];
NPI(3).dates  = {'01-02-2021','01-10-2021',};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 1;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)


% the effect of intervention (NPI) on latent states
%==========================================================================
for i = 1:numel(NPI)
    
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X]  = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
    
end

spm_figure('GetWin','outcomes'); clf;subplot(2,1,1), hold on
%--------------------------------------------------------------------------
spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M);
for i = 1:numel(NPI)
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
end

return


%% scenario modelling for an early relaxation of lockdown on March 8 2021
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-08-2021'};         % duration of epidemic

% scenario one: relaxation of (third) lockdown on March 8, 2021
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'qua'};
NPI(1).Q      = 8;
NPI(1).dates  = {'08-03-2021','01-08-2021'};

% scenario two: maintenance of (second November 5) lockdown on December 2
% 2020 until January 4,2021
%--------------------------------------------------------------------------
% NPI(1).period = period;
% NPI(1).param  = {'sde'};
% NPI(1).Q      = exp(DCM.Ep.sde + 1/2);
% NPI(1).dates  = {'02-12-2020','04-01-2021'};

% scenario three: early second lockdown on October 26, 2020
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'sde'};
NPI(1).Q      = exp(DCM.Ep.sde + 1);
NPI(1).dates  = {'26-10-2020','27-11-2020'};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 1;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

% the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
for j = 1:2, subplot(3,2,j), hold on, end
for j = 5:12,subplot(6,2,j), hold on, end
[Z,X] = spm_SARS_gen(Ep,M,u,NPI);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
subplot(2,1,1), hold on, u = 1;
[m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M);
[m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI);
subplot(2,1,2), hold on, u = 14;
spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M);
spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI);

m1  = m1(end);
m2  = m2(end);
c1  = c1(end);
c2  = c2(end);
d   = m1 - m2;

sprintf('1st %.0f (%.0f to %.0f)',m1,m1 - 1.64*sqrt(c1),m1 + 1.64*sqrt(c1))
sprintf('2nd %.0f (%.0f to %.0f)',m2,m2 - 1.64*sqrt(c2),m2 + 1.64*sqrt(c2))
sprintf('Dif %.0f (%.0f to %.0f)',d, d  - 1.64*sqrt(c1),d  + 1.64*sqrt(c1))

return



%% scenario modelling for increased transmissibility (80%) on 1 May 2021
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-01-2022'};         % duration of epidemic

% scenario one: relaxation of (third) lockdown on March 8, 2021
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'trn','trm'};
NPI(1).Q      = [exp(DCM.Ep.trn + log(1.8)) exp(DCM.Ep.trm + log(1.8))];
NPI(1).dates  = {'01-05-2021',period{2}};


% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 1;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u(1))),u)

% the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
for j = 1:2, subplot(3,2,j), hold on, end
for j = 5:12,subplot(6,2,j), hold on, end
[Z,X] = spm_SARS_gen(Ep,M,u,NPI);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
subplot(2,1,1), hold on, u = 1;
[m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M);
[m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI);
subplot(2,1,2), hold on, u = 14;
spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M);
spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI);

m1  = m1(end);
m2  = m2(end);
c1  = c1(end);
c2  = c2(end);
d   = m1 - m2;

sprintf('1st %.0f (%.0f to %.0f)',m1,m1 - 1.64*sqrt(c1),m1 + 1.64*sqrt(c1))
sprintf('2nd %.0f (%.0f to %.0f)',m2,m2 - 1.64*sqrt(c2),m2 + 1.64*sqrt(c2))
sprintf('Dif %.0f (%.0f to %.0f)',d, d  - 1.64*sqrt(c1),d  + 1.64*sqrt(c1))

return




%% scenario modelling for relaxing restrictions on 21st of June 2021
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-01-2022'};          % duration of epidemic

% scenario: relaxation of restrictions on 21st of June
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'qua'};
NPI(1).Q      = 1;
NPI(1).dates  = {'21-06-2021','01-10-2021',};

NPI(2).period = period;
NPI(2).param  = {'qua'};
NPI(2).Q      = 1;
NPI(2).dates  = {'21-07-2021','01-10-2021',};

NPI(3).period = period;
NPI(3).param  = {'qua'};
NPI(3).Q      = 1;
NPI(3).dates  = {'21-08-2021','01-10-2021',};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 22;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u(1))),u)

%% the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X] = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
end


%% confidence intervals: fatalities (u = 1) and retail (u = 16)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
subplot(2,1,1), hold on, u = 16;
[m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M);
subplot(2,1,2), hold on, u = 27;
spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M);

for i = 1:numel(NPI)
    subplot(2,1,1), hold on, u = 16;
    [m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    subplot(2,1,2), hold on, u = 27;
    spm_SARS_ci(Ep,Cp,S(:,find(U == u)),u,M,NPI(i));
    
    m1  = m1(end);
    m2  = m2(end);
    c1  = c1(end);
    c2  = c2(end);
    d   = m1 - m2;
    
    sprintf('1st %.0f (%.0f to %.0f)',m1,m1 - 1.64*sqrt(c1),m1 + 1.64*sqrt(c1))
    sprintf('2nd %.0f (%.0f to %.0f)',m2,m2 - 1.64*sqrt(c2),m2 + 1.64*sqrt(c2))
    sprintf('Dif %.0f (%.0f to %.0f)',d, d  - 1.64*sqrt(c1),d  + 1.64*sqrt(c1))
end
return



%% scenario modelling for relaxing restrictions on 19 July 2021
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-04-2022'};          % duration of epidemic

% scenario: relaxation of restrictions on 21st of June
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'qua'};
NPI(1).Q      = 1;
NPI(1).dates  = {'19-07-2021','01-10-2021',};

NPI(2).period = period;
NPI(2).param  = {'qua'};
NPI(2).Q      = 1;
NPI(2).dates  = {'19-08-2021','01-10-2021',};

NPI(3).period = period;
NPI(3).param  = {'qua'};
NPI(3).Q      = 1;
NPI(3).dates  = {'19-09-2021','01-10-2021',};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 2;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

% the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
u     = 2;
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X] = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
end


%% confidence intervals: fatalities (u = 1), cases (u = 2), retail (u = 14)
% admissions (u = 16) and Long COVID (u = 28)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u  = [1 2 14 16 28];
nu = numel(u);

for j = 1:numel(NPI)
    for i = 1:nu

        subplot(nu,1,i), hold on
        [m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M);
        subplot(nu,1,i), hold on
        [m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M,NPI(j));
        
        m1  = m1(end);
        m2  = m2(end);
        c1  = c1(end);
        c2  = c2(end);
        d   = m2 - m1;
        
        sprintf('1st %.0f (%.0f to %.0f)',m1,m1 - 1.64*sqrt(c1),m1 + 1.64*sqrt(c1))
        sprintf('2nd %.0f (%.0f to %.0f)',m2,m2 - 1.64*sqrt(c2),m2 + 1.64*sqrt(c2))
        sprintf('Dif %.0f (%.2f percent)',d, 100*d/m1)
    end
end

return



%% scenario modelling the impact of vaccinating children
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-04-2023'};          % duration of epidemic

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 2;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

% suppress vaccination of children
%--------------------------------------------------------------------------
Ep.rol(1) = -16;
Ep.fol(1) = -16;

for j = 1:2, subplot(3,2,j), hold on, end
for j = 5:12,subplot(6,2,j), hold on, end

[Z,X] = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)


%% confidence intervals: fatalities (u = 1), cases (u = 2), retail (u = 14)
% admissions (u = 16) and Long COVID (u = 28)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u  = [1 2 14 16 28];
nu = numel(u);
for i = 1:nu
    
    Ep.rol(1) = DCM.Ep.rol(1);
    Ep.fol(1) = DCM.Ep.fol(1);
    subplot(nu,1,i), hold on
    [m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M);
    
    % suppress vaccination of children
    %--------------------------------------------------------------------------
    Ep.rol(1) = -16;
    Ep.fol(1) = -16;
    subplot(nu,1,i), hold on
    [m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M);
    
    m1  = m1(end);
    m2  = m2(end);
    c1  = c1(end);
    c2  = c2(end);
    d   = m1 - m2;
    
    sprintf('1st %.0f (%.0f to %.0f)',m1,m1 - 1.64*sqrt(c1),m1 + 1.64*sqrt(c1))
    sprintf('2nd %.0f (%.0f to %.0f)',m2,m2 - 1.64*sqrt(c2),m2 + 1.64*sqrt(c2))
    sprintf('Dif %.0f (%.0f to %.0f)',d, d  - 1.64*sqrt(c1),d  + 1.64*sqrt(c1))
    
end

return

%% scenario modelling for relaxing restrictions on 1 november 2021
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-08-2022'};          % duration of epidemic

% scenario: relaxation of restrictions on 21st of June
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'sde'};
NPI(1).Q      = exp(DCM.Ep.sde + 2);
NPI(1).dates  = {'01-11-2021','01-12-2021',};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 2;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

% the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
u     = 2;
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X] = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
end


%% confidence intervals: fatalities (u = 1), cases (u = 2), retail (u = 14)
% admissions (u = 16) and Long COVID (u = 28)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u  = [1 2 14 16 28 32];
nu = numel(u);

for j = 1:numel(NPI)
    for i = 1:nu

        subplot(nu,1,i), hold on
        [m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M);
        subplot(nu,1,i), hold on
        [m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M,NPI(j));
        
        plot([1 1]*datenum(NPI(j).dates{1},'dd-mm-yyyy'),get(gca,'YLim'),':k')
        plot([1 1]*datenum(NPI(j).dates{2},'dd-mm-yyyy'),get(gca,'YLim'),':k')

        k      = datenum(NPI(j).dates{1},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
        e1     = m1(end) - m1(k);
        e2     = m2(end) - m2(k);
        d(i,j) = e2 - e1;
         
        xlabel(sprintf('Date [difference: %.0f (%.2f%s)]',-d(i,j), -100*d(i,j)/e1, '%'))
    end
end

% hospital admissions
%--------------------------------------------------------------------------
clc
dif = -d(4,1);
sprintf('hospital admissions: %.0f x £50,000 = £%.2fM',dif, dif*50000/1e6)

% Gross domestic product
%--------------------------------------------------------------------------
dif = -d(6,1)/(M.T - k);   % percent GDP loss per day for this quarter
sprintf('Quarterly GDP: %.2f%s x £500B = £%.2fB',dif, '%', 500*dif/100)


return

%% scenario modelling for increased transmissibility on 28 november 2021
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (transmission)
%==========================================================================
period = {DCM.M.date,'01-08-2022'};          % duration of epidemic

% scenario: relaxation of restrictions on 21st of June
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'trn'};
NPI(1).Q      = exp(DCM.Ep.trn + log(1.4));
NPI(1).dates  = {'28-11-2021','01-08-2022'};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 2;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

% the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
u     = 2;
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X] = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
end


% confidence intervals: fatalities (u = 1), cases (u = 2), retail (u = 14)
% admissions (u = 16) and Long COVID (u = 28)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u  = [1 2 14 16 28 32];
nu = numel(u);

for j = 1:numel(NPI)
    for i = 1:nu

        subplot(nu,1,i), hold on
        [m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M);
        subplot(nu,1,i), hold on
        [m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M,NPI(j));
        
        plot([1 1]*datenum(NPI(j).dates{1},'dd-mm-yyyy'),get(gca,'YLim'),':k')
        plot([1 1]*datenum(NPI(j).dates{2},'dd-mm-yyyy'),get(gca,'YLim'),':k')

        k      = datenum(NPI(j).dates{1},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
        e1     = m1(end) - m1(k);
        e2     = m2(end) - m2(k);
        d(i,j) = e2 - e1;
         
        xlabel(sprintf('Date [difference: %.0f (%.2f%s)]',-d(i,j), -100*d(i,j)/e1, '%'))
    end
end




%% scenario modelling for increased transmissibility on 15 Jan 2022
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation of social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-09-2022'};          % duration of epidemic

% scenario: relaxation of restrictions on 21st of June
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'qua'};
NPI(1).Q      = 8;
NPI(1).dates  = {'26-01-2022',period{2}};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 14;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

% and the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X] = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
end


%% confidence intervals: fatalities (u = 1), cases (u = 2), retail (u = 14)
% admissions (u = 16) and Long COVID (u = 28)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u  = [1 2 16 28 13 32];
nu = numel(u);

for j = 1:numel(NPI)
    for i = 1:nu

        subplot(nu/2,2,i), hold on
        [m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M);
        axis square
        subplot(nu/2,2,i), hold on
        [m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M,NPI(j));
        axis square
        
        plot([1 1]*datenum(NPI(j).dates{1},'dd-mm-yyyy'),get(gca,'YLim'),':k')
        plot([1 1]*datenum(NPI(j).dates{2},'dd-mm-yyyy'),get(gca,'YLim'),':k')

        k      = datenum(NPI(j).dates{1},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
        e1     = m1(end) - m1(k);
        e2     = m2(end) - m2(k);
        d(i,j) = e2 - e1;
         
        xlabel(sprintf('Date [difference: %.0f (%.2f%s)]',d(i,j), 100*d(i,j)/e1, '%'))
    end
end

%% hospital admissions
%--------------------------------------------------------------------------
clc
dif = d(3,1);
fprintf('Cost of hospital admissions: %.0f x £50,000 = £%.2fB\n\n',dif, dif*50000/1e9)

% £-QALY = £60,000 (green Book), 1 death = 7.6 QALY
%--------------------------------------------------------------------------
dif = d(1,1);
fprintf('Cost (£-QALY = £60,000) of deaths (1 death = 7.6 QALY): %.0f x £456,000 = £%.2fB\n\n',dif, dif*456000/1e9)


% Gross domestic product
%--------------------------------------------------------------------------
days = M.T - k;
dif = d(6,1)/days;   % percent GDP loss per day
gdp = 560/(365/4);
fprintf('Loss to GDP: %.2f%s per day x %.0f days @ £560B per quarter = £%.2fB\n\n',dif, '%', days, days*gdp*dif/100)





%% scenario modelling for rescinding self-isolation on 16 March 2022
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% parametric intervention (relaxation of social distancing threshold)
%==========================================================================
period = {DCM.M.date,'01-03-2023'};          % duration of epidemic

% scenario: relaxation of restrictions on 21st of June
%--------------------------------------------------------------------------
NPI(1).period = period;
NPI(1).param  = {'iso'};
NPI(1).Q      = 2;
NPI(1).dates  = {'16-03-2022',period{2}};

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
M.T    = datenum(period{2},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
u      = 14;
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S(:,find(U == u)),u)

% and the effect of intervention (NPI) on latent states
%--------------------------------------------------------------------------
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X] = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S(:,find(U == u)),u)
end


%% confidence intervals: fatalities (u = 1), cases (u = 2), retail (u = 14)
% admissions (u = 16) and Long COVID (u = 28)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u  = [1 2 16 28 13 32];
nu = numel(u);

for j = 1:numel(NPI)
    for i = 1:nu

        subplot(nu/2,2,i), hold on
        [m1,c1] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M);
        axis square
        subplot(nu/2,2,i), hold on
        [m2,c2] = spm_SARS_ci(Ep,Cp,S(:,find(U == u(i),1)),u(i),M,NPI(j));
        axis square
        
        plot([1 1]*datenum(NPI(j).dates{1},'dd-mm-yyyy'),get(gca,'YLim'),':k')
        plot([1 1]*datenum(NPI(j).dates{2},'dd-mm-yyyy'),get(gca,'YLim'),':k')

        k      = datenum(NPI(j).dates{1},'dd-mm-yyyy') - datenum(period{1},'dd-mm-yyyy');
        e1     = m1(end) - m1(k);
        e2     = m2(end) - m2(k);
        d(i,j) = e2 - e1;
         
        xlabel(sprintf('Date [difference: %.0f (%.2f%s)]',d(i,j), 100*d(i,j)/e1, '%'))
    end
end

%% hospital admissions
%--------------------------------------------------------------------------
clc
dif = d(3,1);
fprintf('Cost of hospital admissions: %.0f x £50,000 = £%.2fB\n\n',dif, dif*50000/1e9)

% £-QALY = £60,000 (green Book), 1 death = 7.6 QALY
%--------------------------------------------------------------------------
dif = d(1,1);
fprintf('Cost (£-QALY = £60,000) of deaths (1 death = 7.6 QALY): %.0f x £456,000 = £%.2fB\n\n',dif, dif*456000/1e9)


% Gross domestic product
%--------------------------------------------------------------------------
days = M.T - k;
dif = d(6,1)/days;   % percent GDP loss per day
gdp = 560/(365/4);
fprintf('Loss to GDP: %.2f%s per day x %.0f days @ £560B per quarter = £%.2fB\n\n',dif, '%', days, days*gdp*dif/100)



