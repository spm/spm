function DCM = DEM_FIN
% FORMAT DCM = DEM_FIN
%
% Demonstration of COVID-19 modelling using variational Laplace (4 groups)
%__________________________________________________________________________
% This demonstration routine illustrates the complex system modelling of
% financial markets using a generic (Helmholtz Hodge) generative model of
% stochastic chaos under the prior are there exists a pullback attractor.
% This endows the dynamics with a compact functional form; especially when
% using second-order approximations (cf., A generalisation of the Laplace
% approximation).
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% set up and get data
%==========================================================================
close all
clear all
clc
spm_figure('GetWin','SI'); clf;
cd('C:\Users\karl\Dropbox\Fintech')

% set final date for retrospective analyses. Full analysis: '01-Dec-2023'
%==========================================================================
global DATE
date = '31/12/2007';
DATE = date;
%--------------------------------------------------------------------------



%% import data
%==========================================================================
data   = importdata('data_yfinance_daily.csv');

prices = {'SPY', 'VEA', 'AGG', 'VNQ', 'DBC', 'BIL'};
ip     = ismember(data.textdata(1,2:end),prices);
d      = 1;

% create data structure
%--------------------------------------------------------------------------
Y(1).type = 'indicator';
Y(1).unit = 'number/day';
Y(1).date = datenum(data.textdata((1 + d):end,1),'dd/mm/yyyy');
Y(1).Y    = data.data(d:end,~ip);
Y(1).h    = 0;

% create data structure
%--------------------------------------------------------------------------
Y(2).type = 'prices';
Y(2).unit = 'number/day';
Y(2).date = datenum(data.textdata((1 + d):end,1),'dd/mm/yyyy');
Y(2).Y    = data.data(d:end,ip);
Y(2).h    = 0;

subplot(2,1,1)
plot(Y(1).date,Y(1).Y)
datetick('x','mmm-yy')

subplot(2,1,2)
plot(Y(2).date,Y(2).Y)
datetick('x','mmm-yy')

Ny    = size(Y(1).Y,1); 
X     = spm_orth([ones(Ny,1) (1:Ny)']);
for j = 1:size(Y(1).Y,2)

    % test for nature of data
    %----------------------------------------------------------------------
    if all(Y(1).Y(:,j) >= 0)
        y(:,j) = spm_log(Y(1).Y(:,j));
    else
        y(:,j) = Y(1).Y(:,j);
    end

    % detrend (in log sapce)
    %----------------------------------------------------------------------
    B(:,j) = X\y(:,j);
    y(:,j) = y(:,j) - X*B(:,j);

end

% modes
%==========================================================================
yy      = spm_en(y);
[u,s,v] = spm_svd(yy);
s       = diag(s);

d  = y(:,26) > 2;

subplot(2,2,1)
plot3(u(:,1),u(:,2),u(:,3)); hold on, plot3(u(d,1),u(d,2),u(d,3),'.r'), hold off
subplot(2,2,2)
plot3(u(:,2),u(:,3),u(:,4)); hold on, plot3(u(d,2),u(d,3),u(d,4),'.r'), hold off



return


% remove NANs, smooth and sort by date
%==========================================================================
M.date  = '01-Jan-2020';
M.end   = date;
[Y,S]   = spm_COVID_Y(Y,M.date,8,M.end);

% get and set priors
%==========================================================================
[pE,pC] = spm_SARS_priors(nN);
pE.N    = log(N(:));
pC.N    = spm_zeros(pE.N);

% age-specific
%--------------------------------------------------------------------------
j       = nN;
pE.mo   = ones(j,1)*[ 0 1];    % coefficients for transport
pC.mo   = ones(j,2)/8;         % prior variance
pE.wo   = ones(j,1)*[ 0 1];    % coefficients for workplace
pC.wo   = ones(j,2)/8;         % prior variance
pE.gd   = ones(j,1)*[-1 1];    % coefficients gross domestic product
pC.gd   = ones(j,2)/8;         % prior variance

% augment priors with fluctuations: number of basis functions
%--------------------------------------------------------------------------
i       = ceil((datenum(M.end) - datenum(M.date))/64);

pE.mob  = zeros(1,i);          % mobility
pC.mob  = ones(1,i)/1024;      % prior variance
pE.pcr  = zeros(1,i);          % testing
pC.pcr  = ones(1,i)/8;         % prior variance

% reporting lags
%--------------------------------------------------------------------------
lag([Y.U]) = [Y.lag];
pE.lag  = spm_zeros(lag);      % reporting delays
pC.lag  = lag; 

% data structure with vectorised data and covariance components
%--------------------------------------------------------------------------
xY.y  = spm_vec(Y.Y);
nY    = numel(Y);
for i = 1:nY
    Q       = spm_zeros({Y.Q});
    Q{i}    = Y(i).Q;
    xY.Q{i} = spm_cat(spm_diag(Q));   
end
hE    = spm_vec(Y.h);

% model specification
%==========================================================================
M.Nmax = 16;                   % maximum number of iterations
M.G    = @spm_SARS_gen;        % generative function
M.FS   = @(Y)real(sqrt(Y));    % feature selection  (link function)
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = hE;                   % prior expectation  (log-precision)
M.hC   = 1/128;                % prior covariances  (log-precision)
M.T    = Y;                    % data structure

U      = [Y.U];                % outputs to model
A      = [Y.age];              % age bands

% initialisation
%--------------------------------------------------------------------------
if true
    
    % load previous posteriors
    %----------------------------------------------------------------------
    load DCM_UK_tmp
    
    % initialise if previous posteriors are the right size
    %----------------------------------------------------------------------
    M.P   = pE;
    field = fieldnames(M.P);
    for i = 1:numel(field)
        try
            if all(size(pE.(field{i})) == size(DCM.Ep.(field{i})))
                M.P.(field{i}) = DCM.Ep.(field{i});
            end
        end
    end
    clear DCM
end

% model inversion with Variational Laplace (Gauss Newton)
%==========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);

% save in DCM structure
%--------------------------------------------------------------------------
DCM.M  = M;
DCM.Ep = Ep;
DCM.Eh = Eh;
DCM.Cp = Cp;
DCM.F  = F;
DCM.Y  = S;
DCM.U  = U;
DCM.A  = A;

% return if just DCM is required
%--------------------------------------------------------------------------
if nargout, return, end

% otherwise save
%--------------------------------------------------------------------------
save('DCM_UK_tmp.mat','DCM')
cd('C:\Users\karl\Dropbox\Coronavirus')
save('DCM_UK_tmp.mat','DCM')

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs
A   = DCM.A;                                 % age cohort


% posterior predictions
%==========================================================================
spm_figure('GetWin','United Kingdom'); clf;
%--------------------------------------------------------------------------
M.T       = 128 + datenum(M.end) - datenum(M.date);
u         = [find(U == 1,1) find(U == 2,1) find(U == 3,1)];
[H,X,~,R] = spm_SARS_gen(Ep,M,[1 2 3]);
spm_SARS_plot(H,X,S(:,u),[1 2 3])

spm_figure('GetWin','outcomes (1)'); clf;
%--------------------------------------------------------------------------
Y     = DCM.M.T;
j     = 0;
k     = 1;
for i = 1:numel(Y)
    
    j = j + 1;
    subplot(4,2,j)
    spm_SARS_ci(Ep,Cp,S(:,i),U(i),M,[],A(i)); hold on
    plot(datenum(M.end)*[1,1],get(gca,'YLim'),'-.b')
    title(Y(i).type,'FontSize',14), ylabel(Y(i).unit)

    
    % add R = 1 and current date
    %----------------------------------------------------------------------
    if Y(i).U == 4
        plot(get(gca,'XLim'),[1,1],'-.r')
        plot(datenum(M.end)*[1,1],get(gca,'YLim'),'-.b')
        set(gca,'YLim',[0 5]), ylabel('ratio')
    end
    
    % hold plot
    %----------------------------------------------------------------------
    if Y(i).hold
        j = j - 1; hold on
    end

    % new figure
    %----------------------------------------------------------------------
    if j == 8
        k = k + 1;
        j = 0;
        spm_figure('GetWin',sprintf('outcomes (%i)',k)); clf;
    end
    
end

% time varying parameters
%==========================================================================

% infection fatality ratios (%)
%--------------------------------------------------------------------------
j     = j + 1;
subplot(4,2,j)
for i = 1:numel(R)
    spm_SARS_ci(Ep,Cp,[],21,M,[],i); hold on
end
ylabel('percent'), title('Infection fatality ratio','FontSize',14)

% transmission risk
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j), hold off
plot([R{end}.Ptrn]), hold on
plot(erf([R{end}.Ptra]*exp(Ep.trn))), spm_axis tight
title('Transmission risk','FontSize',14)
xlabel('days'),ylabel('probability')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), box off

% contact rate
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j), hold on
for i = 1:numel(R)
    plot([R{i}.Pout])
end
spm_axis tight
title('Contact rate','FontSize',14)
xlabel('days'),ylabel('probability')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), box off

% case fatality ratio
%--------------------------------------------------------------------------
j    = j + 1;
subplot(4,2,j), hold on
for i = 1:numel(R)
    plot(100 * [R{i}.Pfat])
end
spm_axis tight
title('Fatality risk | ARDS','FontSize',14)
xlabel('days'),ylabel('percent')
hold on, plot([1,1]*size(DCM.Y,1),[0,1/2],':'), box off
legend({'<15yrs','15-35yrs','35-70yrs','>70yrs'})


%% long-term forecasts (12 months from the current data)
%==========================================================================
spm_figure('GetWin','outcomes (4)'); clf

Ep  = DCM.Ep;
Cp  = DCM.Cp;
M   = DCM.M;
M.T = 30*12 + datenum(M.end) - datenum(M.date);
t   = (1:M.T) + datenum(M.date);

% infection fatality ratios (%)
%--------------------------------------------------------------------------
subplot(3,1,1)
spm_SARS_ci(Ep,Cp,[],19,M); hold on
spm_SARS_ci(Ep,Cp,[],20,M); hold on
ylabel('cases per 100,000'), title('Incidence per 100,000','FontSize',14)
plot(datenum(M.end)*[1,1],get(gca,'YLim'),':')
legend({'CI per day','actual cases per day','CI per week','confirmed cases per week'})

subplot(3,1,2)
spm_SARS_ci(Ep,Cp,[],34,M); hold on
spm_SARS_ci(Ep,Cp,S(:,1),2, M); hold on
title('Actual and confirmed daily cases','FontSize',14)
plot(datenum(M.end)*[1,1],get(gca,'YLim'),':')
legend({'CI per day','actual cases','CI per day','confirmed cases'})

subplot(3,1,3)
spm_SARS_ci(Ep,Cp,[],28,M); hold on
plot(datenum(M.end)*[1,1],get(gca,'YLim'),':')
legend({'CI per day','Incidence'})


%% switch windows
%--------------------------------------------------------------------------
spm_figure('GetWin','long-term (1)'); clf

% fatalities
%--------------------------------------------------------------------------
subplot(2,1,1)

i   = find(DCM.U == 1,1);  D = DCM.Y(:,i); spm_SARS_ci(Ep,Cp,D,1,M);  hold on
i   = find(DCM.U == 15,1); D = DCM.Y(:,i); spm_SARS_ci(Ep,Cp,D,15,M); hold on

plot(datenum(M.end)*[1,1],get(gca,'YLim'),':k')
ylabel('number per day'), title('Daily deaths','FontSize',14)
legend({'CI 28-day','PCR test within 28 days','ONS','CI certified','certified deaths'})
legend boxoff
drawnow

% lockdown and mobility
%--------------------------------------------------------------------------
subplot(2,1,2)
i       = find(DCM.U == 14,1);
try
    [~,~,q] = spm_SARS_ci(Ep,Cp,DCM.Y(:,i),14,M,[],DCM.A(i)); hold on
catch
    [~,~,q] = spm_SARS_ci(Ep,Cp,[],14,M); hold on
end

% thresholds
%--------------------------------------------------------------------------
% q  = spm_SARS_gen(Ep,M,14); plot(t,q); hold on
%--------------------------------------------------------------------------
u1   = datenum('10-May-2020','dd-mmm-yyyy') - t(1) + 1;
u2   = datenum('10-Aug-2020','dd-mmm-yyyy') - t(1) + 1;
u3   = datenum('10-Sep-2020','dd-mmm-yyyy') - t(1) + 1;
U    = sort([0 q(u1) q(u2) q(u3)]); U(end) = 90;
dstr = datestr(t,'dd-mmm');

% loop over levels
%==========================================================================
for i = 1:numel(U)
    
    % intervals for this level
    %----------------------------------------------------------------------
    if i == 1
        j  = find(q <= U(i + 1));
    elseif i == numel(U)
        j  = find(q >= U(i));
    else
        j  = find(q >= U(i) & q <= U(i + 1));
    end
    
    % Timeline
    %----------------------------------------------------------------------
    for k = 1:numel(j)
        try
            fill(t(j(k) + [0 1 1 0]),[0 0 1 1]*32,'r', ...
                'FaceAlpha',(numel(U) - i)/16,'Edgecolor','none')
        end
    end
    
    % label level crossings
    %----------------------------------------------------------------------
    if i <numel(U)
        j = find((q(1:end - 1) <= U(i + 1)) & (q(2:end) > U(i + 1)));
    else
        j = [];
    end
    for k = 1:numel(j)
        text(t(j(k)),i*8,dstr(j(k),:),'Color','k','FontSize',10)
    end
    
    % plot levels
    %----------------------------------------------------------------------
    plot([t(1),t(end)],U(i)*[1,1],':r')
    
end

% UEFA EURO 2020/Dates
%--------------------------------------------------------------------------
% d1 = datenum('11-Jun-2021','dd-mmm-yyyy');
% d2 = datenum('11-Jul-2021','dd-mmm-yyyy');
% plot([d1,d2],[120,120],'k','Linewidth',8)
% text(d1 - 86,120,'EURO 2020','FontSize',10)


ylabel('percent'),  title('Mobility and lockdown','FontSize',14)
legend({'credible interval','mobility (%)','Google workplace'}), legend boxoff
drawnow


%% prevalence and reproduction ratio
%--------------------------------------------------------------------------
spm_figure('GetWin','long-term (2)'); clf

M.T = 30*12   + datenum(M.end) - datenum(M.date);
t   = (1:M.T) + datenum(M.date);

subplot(2,1,1)
i   = find(DCM.U == 4,1);
spm_SARS_ci(Ep,Cp,[],11,M); hold on
[~,~,q,c] = spm_SARS_ci(Ep,Cp,DCM.Y(:,i),4,M); hold on

j   = find(t == datenum(M.end));
q   = q(j);
d   = sqrt(c{1}(j,j))*1.64;
str = sprintf('Prevalence and reproduction ratio (%s): R = %.2f (CI %.2f to %.2f)',datestr(date,'dd-mmm-yy'),q,q - d,q + d);

% attack rate, herd immunity and herd immunity threshold
%--------------------------------------------------------------------------
E         = 1 - mean(erf(exp(Ep.vef)))*mean(exp(Ep.ves));
[H,~,~,R] = spm_SARS_gen(Ep,M,[4 29 26 35]);
i         = 8:32;                           % pre-pandemic period
TRN       = [R{1}.Ptrn];                    % transmission risk
TIC       = [R{1}.Tic];                     % incubation period (days)
TIN       = [R{1}.Tin];                     % exposure time (days)
ST        = H(:,4);                         % serial interval (days)
R0        = max(H(i,1));                    % basic reproduction ratio
RT        = R0*TRN(:)/min(TRN(i));          % effective reproduction ratio
HIT       = 100 * (1 - 1./RT)/E;            % herd immunity threshold

% Add R0
%--------------------------------------------------------------------------
alpha   = datenum('20-Sep-2020','dd-mmm-yyyy');
delta   = datenum('20-Mar-2021','dd-mmm-yyyy');
omicron = datenum('15-Nov-2021','dd-mmm-yyyy');

plot(t,RT)
text(alpha,4,'alpha','FontSize',10,'HorizontalAlignment','center')
text(delta,4,'delta','FontSize',10,'HorizontalAlignment','center')
text(omicron,4,'omicron','FontSize',10,'HorizontalAlignment','center')

% add R = 1 and current dateline
%--------------------------------------------------------------------------
plot(get(gca,'XLim'),[1,1],':k')
plot(datenum(M.end)*[1,1],[0 max(RT)],':k')
set(gca,'YLim',[0 10]), ylabel('ratio or percent')
title(str,'FontSize',14)

legend({'90% CI','Prevalence (%)',' ','Effective R-number',...
    'UKHSA estimate','Basic R_{0}'},'Interpreter','tex','location','northwest')
legend boxoff
drawnow

% attack rate, herd immunity and herd immunity threshold
%--------------------------------------------------------------------------
subplot(2,1,2)
spm_SARS_ci(Ep,Cp,[],25,M); hold on
spm_SARS_ci(Ep,Cp,[],26,M); hold on

% effective immunity threshold at 80% contact rates
%--------------------------------------------------------------------------
plot(t,HIT), hold on
spm_SARS_ci(Ep,Cp,[],29,M); hold on

plot(get(gca,'XLim'),[100 100],':k')
plot(datenum(M.end)*[1,1],[0 100],':k')
ylabel('percent'),  title('Attack rate and immunity','FontSize',14)
legend({' ','Attack rate',' ',...
       'Effective population immunity',...
       'Effective immunity threshold',' '....
       'Effective antibodies'},'location','west')
legend boxoff


%% infection fatality ratio (%)
% probability of becoming symptomatic: (1 - Pinf.^(Q{n}.Tin + Q{n}.Tcn))
%----------------------------------------------------------------------
for n = 1:numel(R)
    
    Psev     = [R{n}.Psev]';
    Pfat     = [R{n}.Pfat]';
    IFR(n,1) = 100 * Psev(j).*Pfat(j);

    Psev = Psev*exp(Ep.lnk(n));         % vaccine efficiency: pathogenicity
    Pfat = Pfat*exp(Ep.lnf(1));         % vaccine efficiency: fatality
    IFR(n,2) = 100 * Psev(j).*Pfat(j);
end

fprintf('\n\n(Symptomatic) IFR: likelihood of dying from COVID-19 (%s) \n','%');
fprintf('            Age group  0-14     15-34    35-69   70+  \n');
fprintf('   recent vaccination: %.3f    %.3f    %.3f   %.3f \n',IFR(:,2));
fprintf('no vaccine protection: %.3f    %.3f    %.3f   %.3f \n\n',IFR(:,1));


% report vaccine efficiency
%--------------------------------------------------------------------------
q   = Ep.vef(end);
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.vef(end))*1.64;
qE  = 100*(1 - erf(exp(q)));
qL  = 100*(1 - erf(exp(q + d)));
qU  = 100*(1 - erf(exp(q - d)));
fprintf('preventing infection: %.1f%s (CI %.1f to %.1f)\n',qE,'%',qL,qU)
q   = Ep.ves(end);
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.ves(end))*1.64;
qE  = 100*(1 - exp(q));
qL  = 100*(1 - exp(q + d));
qU  = 100*(1 - exp(q - d));
fprintf('preventing transmission following infection %.1f%s (CI %.1f to %.1f)\n',qE,'%',qL,qU)
q   = Ep.lnk(2);
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.lnk(2))*1.64;
qE  = 100*(1 - exp(q));
qL  = 100*(1 - exp(q + d));
qU  = 100*(1 - exp(q - d));
fprintf('preventing serious illness when symptomatic (age 15-34) %.1f%s (CI %.1f to %.1f)\n',qE,'%',qL,qU)
q   = Ep.lnk(3);
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.lnk(3))*1.64;
qE  = 100*(1 - exp(q));
qL  = 100*(1 - exp(q + d));
qU  = 100*(1 - exp(q - d));
fprintf('preventing serious illness when symptomatic (age 35-70) %.1f%s (CI %.1f to %.1f)\n',qE,'%',qL,qU)
q   = Ep.lnf(end);
d   = spm_unvec(diag(Cp),Ep);
d   = sqrt(d.lnf(end))*1.64;
qE  = 100*(1 - exp(q));
qL  = 100*(1 - exp(q + d));
qU  = 100*(1 - exp(q - d));
fprintf('preventing fatality when seriously ill %.1f%s (CI %.1f to %.1f)\n',qE,'%',qL,qU)
disp(' ')

disp('The corresponding cumulative (vaccinated vs. unvaccinated) risks are:')
disp(' ')

q      = (mean(exp(Ep.Tin)) + mean(exp(Ep.Tcn))*mean(exp(Ep.ves))) /...
         (mean(exp(Ep.Tin)) + mean(exp(Ep.Tcn)));

infect = mean(erf(exp(Ep.vef)));
mild   = q*infect;
severe = mean(exp(Ep.lnk))*mild;
death  = mean(exp(Ep.lnf))*severe;
fprintf('relative risk of infection %.1f%s\n',      infect*100,'%')
fprintf('relative risk of mild illness %.1f%s\n',   mild*100,  '%')
fprintf('relative risk of severe illness %.1f%s\n', severe*100,'%')
fprintf('relative risk of fatality %.1f%s\n',       death*100, '%')
disp(' ')

fprintf('%s%.1f%s\n\n','For example, vaccination reduces the risk of being infected and subsequently developing a severe illness to ',severe*100,'% of the risk prior to vaccination. The risk then increases slowly until the next vaccination.')

% mean serial interval
%--------------------------------------------------------------------------
fprintf('%s%.1f%s\n','The average serial interval has fallen to ',ST(j),' days')

% exposure time
%--------------------------------------------------------------------------
fprintf('%s%.1f%s\n','The average period between exposure and becoming infectious is ',TIN(j),' days')

% mean incubation (asymptomatic).
%--------------------------------------------------------------------------
fprintf('%s%.1f%s\n\n','The average asymptomatic period is currently ',TIC(j),' days')

% report transmissibility and basic reproduction number
%--------------------------------------------------------------------------
disp('relative transmissibility');
disp(100*TRN(j)/mean(TRN(1:j)))
disp('basic reproduction number');
disp(RT(j))
disp(' ')

%% ancillary predictions
%--------------------------------------------------------------------------
% Z      = spm_SARS_gen(Ep,M,[1,16,31]);
% 
% fprintf('vaccine effectiveness (prevalence of infection) %.1f%s\n', Z(j,3),'%')
% disp(' ')
% 
% [m,d] = max(Z(j:end,2));
% fprintf('peak hospital admissions:  %.0f on %s\n', m,datestr(t(j + d)))
% disp(' ')
% 
% [m,d] = max(Z(j:end,1));
% fprintf('peak (28-day) death rates:  %.0f on %s\n', m,datestr(t(j + d)))
% disp(' ')

%% save figures
%--------------------------------------------------------------------------
spm_figure('GetWin','outcomes (1)');
savefig(gcf,'Fig1_tmp')

spm_figure('GetWin','outcomes (2)');
savefig(gcf,'Fig2_tmp')

spm_figure('GetWin','outcomes (3)');
savefig(gcf,'Fig3_tmp')

spm_figure('GetWin','outcomes (4)');
savefig(gcf,'Fig4_tmp')

spm_figure('GetWin','United Kingdom');
savefig(gcf,'Fig5_tmp')

spm_figure('GetWin','long-term (2)');
savefig(gcf,'Fig6_tmp')

spm_figure('GetWin','long-term (1)');
savefig(gcf,'Fig7_tmp')

% Table
%--------------------------------------------------------------------------
disp(spm_COVID_table(DCM.Ep,DCM.Cp,DCM.M))

return



%% NOTES
% posthoc fitting of age cohort-specific parameters
%==========================================================================
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs
Y   = DCM.M.T; 

xY.y  = spm_vec(Y.Y);
xY.Q  = spm_Ce([Y.n]);

% empirical priors
%--------------------------------------------------------------------------
pE    = Ep;
pC    = spm_zeros(DCM.M.pC);

% augment priors
%--------------------------------------------------------------------------
nN      = numel(Ep.N);
k       = numel(Ep.tra)*2;
pE.mob  = zeros(nN,k);                       % fluctuations in mobility
pC.mob  = ones(nN,k)/8;                      % prior variance

% augment posteriors
%----------------------------------------------------------------------
M.pE      = pE;                    % empirical prior expectation
M.pC      = pC;                    % fix expectations
[Ep,Cp]   = spm_nlsi_GN(M,U,xY);   % new posterior expectation

return

% Notes for precision – increasing precision of recent (64 day) data
%==========================================================================
nY    = numel(Y);
Q     = [];
rdate = datenum(M.end) - 64;
for i = 1:numel(Y)
    j      = Y(i).date > rdate;
    q      = zeros(Y(i).n,nY);
    q(:,i) = 1;
    q(j,i) = 16;
    Q      = [Q; q];
end
nQ    = size(Q,1);
for i = 1:numel(Y)
    xY.Q{i} = sparse(1:nQ,1:nQ,Q(:,i));
end

% fluctuations (adiabatic mean field approximation)
%==========================================================================
for f = 1:numel(fluct)
    
    % augment priors
    %----------------------------------------------------------------------
    pE.(fluct{f})   = zeros(1,16);           % add new prior expectation
    pC.(fluct{f})   =  ones(1,16);           % add new prior covariance
    
    % augment posteriors
    %----------------------------------------------------------------------
    i               = 1:size(Cp,1);          % number of parameters
    C               = Cp;                    % empirical prior covariance
    M.pE            = Ep;                    % empirical prior expectation
    M.pC            = spm_zeros(M.pC);       % fix expectations
    M.pE.(fluct{f}) = zeros(1,16);           % add new prior expectation
    M.pC.(fluct{f}) =  ones(1,16);           % add new prior covariance
    
    [Ep,Cp,Eh,Ff]   = spm_nlsi_GN(M,U,xY);   % new posterior expectation
    Cp(i,i)         = C;                     % new posterior covariance
    F               = F + Ff;                % free energy
    
    % save priors
    %----------------------------------------------------------------------
    M.pE   = pE;
    M.pC   = pC;
    
end

return

%% age-group specific inversion
%--------------------------------------------------------------------------
a = 2;
N = N(a);
i = [Y.age] == a;
Y = Y(i);
S = S(:,i);
for i = 1:numel(Y)
    Y(i).age = 0;
end

% retrieve age-specific priors
%--------------------------------------------------------------------------
pE.Nin = pE.Nin(a,a);
pE.Nou = pE.Nou(a,a);
pC.Nin = pC.Nin(a,a);
pC.Nou = pC.Nou(a,a);
field = fieldnames(pE);
for i = 1:numel(field)
    if size(pE.(field{i}),1) > 1
        pE.(field{i}) = pE.(field{i})(a,:);
        pC.(field{i}) = pC.(field{i})(a,:);
    end
end


%% sensitivity analysis in terms of partial derivatives
%--------------------------------------------------------------------------
i     = spm_fieldindices(Ep,'oth');
nP    = spm_length(Ep);
V     = spm_speye(nP,i);
dYdP  = spm_diff(M.G,Ep,M,4,1,{V});

% plot results
%--------------------------------------------------------------------------
subplot(2,1,1)
plot(dYdP)
subplot(2,1,2)
plot(cumsum(dYdP))
ylabel('First-order sensitivity','FontSize',16), box off
spm_fieldindices(Ep,9)




%% Interventions
%==========================================================================
% Y(:,1)  - Daily deaths (28 days)
% Y(:,2)  - Daily confirmed cases
% Y(:,3)  - Mechanical ventilation
% Y(:,4)  - Reproduction ratio (R)
% Y(:,5)  - Seroprevalence {%}
% Y(:,6)  - PCR testing rate
% Y(:,7)  - Risk of infection (%)
% Y(:,8)  - Prevalence (true) {%}
% Y(:,9)  - Daily contacts
% Y(:,10) - Daily incidence (%)
% Y(:,11) - Prevalence (positivity){%}
% Y(:,12) - Number symptomatic
% Y(:,13) - Mobility (%)
% Y(:,14) - Workplace (%)
% Y(:,15) - Certified deaths
% Y(:,16) - Hospital admissions
% Y(:,17) - Hospital deaths
% Y(:,18) - Non-hospital deaths
% Y(:,19) - Daily incidence (per hundred thousand)
% Y(:,20) - Weekly confirmed cases (per hundred thousand)
% Y(:,21) - Infection fatality ratio (%)
% Y(:,22) - Cumulative first dose
% Y(:,23) - PCR case positivity (%)
% Y(:,24) - Lateral flow tests
% Y(:,25) - Cumulative attack rate
% Y(:,26) - Population immunity (total)
% Y(:,27) - Hospital cases
% Y(:,28) - Incidence of Long Covid
% Y(:,29) - Proportion vaccinated
% Y(:,30) - Cumulative admissions
% Y(:,31) - Vaccine effectiveness (prevalence)
% Y(:,32) - Gross domestic product

% get posterior estimates
%--------------------------------------------------------------------------
clear
DCM = load('DCM_UK.mat','DCM');
DCM = DCM.DCM;

% unpack model and posterior expectations
%--------------------------------------------------------------------------
M   = DCM.M;                                 % model (priors)
Ep  = DCM.Ep;                                % posterior expectation
Cp  = DCM.Cp;                                % posterior covariances
S   = DCM.Y;                                 % smooth timeseries
U   = DCM.U;                                 % indices of outputs
A   = DCM.A;                                 % age cohort

% plot epidemiological trajectories and hold plots
%==========================================================================
spm_figure('GetWin','states'); clf;
%--------------------------------------------------------------------------
T      = datenum(M.end) - datenum(DCM.M.date);
M.T    = T + 360*2;                          % forecast dates
u      = [1];                                % empirical outcome
a      = 0;                                  % age cohort (0 for everyone)
Ep.qua = DCM.Ep.qua + 0;                     % adjusted (log) parameter

[Z,X]  = spm_SARS_gen(Ep,M,u,[],a); % posterior prediction

% and plot
%--------------------------------------------------------------------------
try
    j     = [];
    for i = 1:numel(u)
        j = [j find(U == u(i) & A == a(i))];  
    end
    spm_SARS_plot(Z,X,S(:,j),u)
catch
    spm_SARS_plot(Z,X,[],u)
end
subplot(3,2,1), hold on
plot([T T]/7,[min(Z) max(Z)],'--')


