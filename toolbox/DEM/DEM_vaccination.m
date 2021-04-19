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
% $Id: DEM_Dispatches.m 8015 2020-11-24 10:47:41Z karl $


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


%% scenario modelling for lockdown and vaccinationin response to the 
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

% scenario three: early second lockdown on October 21, 2020
%--------------------------------------------------------------------------
% NPI(1).period = period;
% NPI(1).param  = {'sde'};
% NPI(1).Q      = exp(DCM.Ep.sde + 1);
% NPI(1).dates  = {'21-10-2020','05-11-2020'};

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


