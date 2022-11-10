function Tab = DEM_Stephen
% FORMAT Tab = DEM_Stephen
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

% get empirical priors from fitting national data
%==========================================================================
if false
    DCM = DEM_COVID_UK;
    save('DCM_UK.mat','DCM')
else
    DCM = load('DCM_UK.mat','DCM');
    DCM = DCM.DCM;
end

% parametric intervention
%==========================================================================
period     = {DCM.M.date,'03-04-2021'};         % duration of epidemic

% NPI.period = period;                          % duration of epidemic
% NPI.param  = 'ttt';                           % parameter to change
% NPI.Q      = 0.24;                            % new value
% NPI.dates  = {'01-11-2020','01-12-2020'};     % dates of implementation

NPI.period = period;
NPI.param  = {'sde','vac'};
NPI.Q      = [exp(DCM.Ep.sde)*1.5 0.5];
NPI.dates  = {'24-12-2020','01-01-2021'};


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
u      = U([1 2 3]);
[Z,X]  = spm_SARS_gen(Ep,M,u);
spm_SARS_plot(Z,X,S,u)

% quantify the effect of efficient intervention (NPI)
%==========================================================================
for i = 1:numel(NPI)
    
    for j = 1:2, subplot(3,2,j), hold on, end
    for j = 5:12,subplot(6,2,j), hold on, end
    
    [Z,X]  = spm_SARS_gen(Ep,M,u,NPI(i));
    spm_SARS_plot(Z,X,S,u)
    
end


% quantify the effect of efficient intervention (NPI)
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes'); clf;
%--------------------------------------------------------------------------
u     = [2 9 8];
nu    = numel(u);
for j = 1:nu
    subplot(nu,1,j), spm_SARS_ci(Ep,Cp,S(:,u(j)),U(u(j)),M);  hold on
    for i = 1:numel(NPI)
        subplot(nu,1,j), spm_SARS_ci(Ep,Cp,S(:,u(j)),U(u(j)),M,NPI(i));
    end
end

% specify outputs of interest ( with indices u)
%==========================================================================
% outputs
%--------------------------------------------------------------------------
% Y(:,1)  - daily deaths
% Y(:,2)  - daily tests
% Y(:,3)  - CCU occupancy
% Y(:,4)  - reproduction ratio (R)
% Y(:,5)  - seropositive immunity (%)
% Y(:,6)  - PCR testing rate
% Y(:,7)  - contagion risk (%)
% Y(:,8)  - prevalence (%)
% Y(:,9)  - new contacts per day
% Y(:,10) - daily incidence
% Y(:,11) - number infected  
% Y(:,12) - number symptomatic
% Y(:,13) - mobility (%)
% Y(:,14) - work (%)
% Y(:,15) - certified deaths/day
% Y(:,16) - hospitalisation
% Y(:,17) - deaths in hospital
% Y(:,18) - deaths in isolation
% Y(:,19) - deaths > 60 yrs
% Y(:,20) - deaths < 60 yrs
%--------------------------------------------------------------------------
NPI       = NPI(end);
[~,~,str] = spm_SARS_priors;
u  = [15 17 18 19 3 6 10 12 13 14 16 4];

% solve dynamic causal model with (Y0) and without (Yi) intervention
%--------------------------------------------------------------------------
Y0 = spm_SARS_gen(Ep,M,u);
Yi = spm_SARS_gen(Ep,M,u,NPI);

% indices of intervention period over which to evaluate outcomes
%--------------------------------------------------------------------------
i  = datenum(NPI.dates{1},'dd-mm-yyyy'):datenum(period{2},'dd-mm-yyyy');
i  = i - datenum(NPI.period{1},'dd-mm-yyyy');
j  = datenum(period{2},'dd-mm-yyyy') - datenum(NPI.period{1},'dd-mm-yyyy');
E0 = mean(Y0(i,:));
Ei = mean(Yi(i,:));
S0 = sum(Y0(i,:));
Si = sum(Yi(i,:));

% construct a table of outcomes and differences
%--------------------------------------------------------------------------
clear Tab
for i = 1:numel(u)
    Tab{i,1}  = str.outcome{u(i)};   
    Tab{i,2}  = Y0(j,i);
    Tab{i,3}  = E0(i);
    Tab{i,4}  = S0(i);
    Tab{i,5}  = Si(i) - S0(i);
    Tab{i,6}  = 100*(Si(i) - S0(i))/S0(i);

end

VariableNames{1}  = 'outcome';
VariableNames{2}  = period{2};
VariableNames{3}  = 'daily average';
VariableNames{4}  = 'cost (sum)';
VariableNames{5}  = 'excess (sum)';
VariableNames{6}  = 'excess (%)';


Tab = cell2table(Tab);
Tab.Properties.Description  = 'parameter estimates';
Tab.Properties.VariableNames = VariableNames;

% display
%--------------------------------------------------------------------------
fprintf('\n    Intervention from %s to %s\n',NPI.dates{1},NPI.dates{2})
fprintf('\n    Averages and sums from %s to %s\n',NPI.dates{1},period{2})
try
    fprintf('\n   (%s = %g versus %s = %g)\n\n',NPI.param,NPI.Q,NPI.param,exp(Ep.(NPI.param)))
end
disp(Tab)

writetable(Tab,'Table.xlsx')

return






