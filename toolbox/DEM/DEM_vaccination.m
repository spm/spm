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

% get empirical priors from fitting national data
%==========================================================================
if false
    DCM = DEM_COVID_UK;
    save('DCM_UK.mat','DCM')
else
    DCM = load('DCM_UK.mat','DCM');
    DCM = DCM.DCM;
end

% parametric intervention(vaccination and contact tracing)
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
    spm_SARS_plot(Z(:,1:3),X,S,[],u(1:3))
    
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
spm_SARS_plot(Z,[],[],[],23)
for i = 1:numel(NPI)
    for j = 1:2, subplot(3,2,j), hold on, end
    Z     = spm_SARS_gen(Ep,M,23,NPI(i));
    spm_SARS_plot(Z,[],[],[],23)
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

% quantify the effect of efficient intervention (NPI)
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






