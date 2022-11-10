function spm_COVID_PV(DCM,i,T)
% FORMAT spm_COVID_PV(DCM,i,T)
% remove ( > T) data from country ( = i)
%--------------------------------------------------------------------------
% i  - country index
% T  - number of days to withhold
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% use priors from parametric empirical Bayes
%--------------------------------------------------------------------------
pE            = DCM{i}.M.pE;
pC            = DCM{i}.M.pC;
data          = DATA_COVID_JHU;
data(i).death = data(i).death(1:end - T);
data(i).cases = data(i).cases(1:end - T);

% invert (using incomplete data) and plot confidence intervals
%--------------------------------------------------------------------------
Y         = [data(i).death, data(i).cases];
[F,Ep,Cp] = spm_COVID(Y,pE,pC);
fig       = ['predictive validity: ' data(i).country];

spm_figure('GetWin',fig); clf
spm_COVID_ci(Ep,Cp,Y,1:3);

% retrieve and overlay withheld data
%--------------------------------------------------------------------------
data = DATA_COVID_JHU;
Y    = [data(i).death, data(i).cases];
NY   = size(Y,1);
t    = (1:NY)/7;
CY   = cumsum(Y(:,1));
i    = (NY - T):NY;
t    = t(i);

spm_figure('GetWin',fig);
subplot(3,2,1), hold on, plot(t,Y(i,1),'.r','MarkerSize',16)
subplot(3,2,3), hold on, plot(t,Y(i,2),'.r','MarkerSize',16)
subplot(2,2,2), hold on, plot(t,CY(i), '.r','MarkerSize',16)

return