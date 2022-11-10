function DCM = DEM_Dispatches
% FORMAT DCM = DEM_COVID_UK
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

% set up and preliminaries
%==========================================================================
% https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata
% https://www.ndm.ox.ac.uk/covid-19/covid-19-infection-survey/results
% https://coronavirus.data.gov.uk/
% https://covid.joinzoe.com/data#levels-over-time

% get figure and data
%--------------------------------------------------------------------------
DCM = DEM_COVID_UK;

% model inversion with Variational Laplace (Gauss Newton)
%--------------------------------------------------------------------------
M   = DCM.M;
Ep  = DCM.Ep;
Cp  = DCM.Cp;
YS  = DCM.Y;
U   = DCM.U;

% posterior predictions
%==========================================================================
spm_figure('GetWin','testing and cases'); clf;
%--------------------------------------------------------------------------
M.T = datenum('01-08-2021','dd-mm-yyyy') - datenum(M.date,'dd-mm-yyyy');

% positive PCR cases and underlying number of new infections per day 
%--------------------------------------------------------------------------
subplot(2,1,1)
spm_SARS_ci(Ep,Cp,YS(:,1),U(1),M); hold on
spm_SARS_ci(Ep,Cp,[],10,M);
title('Daily new infections and tests','FontSize',14), xlabel('number')
UT  = 4e4;

% add window of opportunity 
%--------------------------------------------------------------------------
plot(get(gca,'XLim'),[1,1]*UT,'-.r')
plot(datenum(date)*[1,1],get(gca,'YLim'),'-.b')


% calibrate the probability of testing to about 500,000 tests per day
%--------------------------------------------------------------------------
% spm_figure('GetWin','enhanced testing'); clf;
% 
% Ep.lim(1) = log(0.01);
% Ep.ons(1) = log(datenum('01-Jul-20','dd-mm-yy') - datenum(M.date,'dd-mm-yyyy'));
% Ep.rat(1) = log(4);
% 
% [Z,X]  = spm_SARS_gen(Ep,M,6);
% spm_SARS_plot(Z,X,[],6)

% show the effects on new PCR cases and infections per day
%--------------------------------------------------------------------------
% spm_figure('GetWin','testing and cases');
% subplot(2,1,2), hold off
% spm_SARS_ci(Ep,Cp,[],U(1),M); hold on
% spm_SARS_ci(Ep,Cp,[],10,M);
% title('With enhanced testing capacity','FontSize',14), xlabel('number')
% axis(a)
%  
% add window of opportunity 
%--------------------------------------------------------------------------
% plot(get(gca,'XLim'),[1,1]*UT,'-.r')
% plot(datenum(date)*[1,1],get(gca,'YLim'),'-.b')


% quantify the effect of efficient FTTIS
%==========================================================================
subplot(2,1,2)
spm_SARS_ci(Ep,Cp,YS(:,2),1,M); hold on

Tdate = {'11-Dec-2020'};
M.FTT = .24;
for i = 1:numel(Tdate)
    M.TTT = datenum(Tdate{i},'dd-mmm-yyyy') - datenum(M.date,'dd-mm-yyyy');
    spm_SARS_ci(Ep,Cp,YS(:,2),1,M);
end

% t   = datenum('01-07-2021','dd-mm-yyyy');
% t0  = t - 12*24;
% set(gca,'XLim',[t0,t])
% set(gca,'YLim',[0,500])
% set(gca,'XTick',[t0:28:t])
% datetick('x','mmm-dd','keeplimits','keepticks')


% quantify the effect of efficient FTTIS
%==========================================================================
spm_figure('GetWin','lockdown and hospitals'); clf;
M     = rmfield(M,{'TTT','FTT'});
subplot(2,1,1); spm_SARS_ci(Ep,Cp,YS(:,3),3,M);  hold on
subplot(2,1,2); spm_SARS_ci(Ep,Cp,YS(:,9),13,M); hold on

M.FTT = .24;
for i = 1:numel(Tdate)
    M.TTT = datenum(Tdate{i},'dd-mmm-yyyy') - datenum(M.date,'dd-mm-yyyy');
    subplot(2,1,1); spm_SARS_ci(Ep,Cp,YS(:,3),3,M);  hold on
    subplot(2,1,2); spm_SARS_ci(Ep,Cp,YS(:,9),13,M); hold on
end

return






