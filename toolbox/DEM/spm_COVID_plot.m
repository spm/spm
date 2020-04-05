function spm_COVID_plot(Y,X,Z)
% Graphics for coronavirus simulations
% FORMAT spm_COVID_plot(Y,X,Z)
% Y      - expected timeseries (i.e., new depths and cases)
% X      - latent (marginal ensemble density) states
% Z      - optional empirical data
%
% This auxiliary routine plots the trajectory of outcome variables
% and underlying latent or hidden states, in the form of marginal densities
% over the four factors that constitute the COVID model. if empirical data
% are supplied, they will be superimposed. if more than four expected or
% predicted outcomes are supplied, a threshold will be superimposed;
% reflecting the typical number of beds available per population cell.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_COVID_plot.m 7811 2020-04-05 12:00:43Z karl $

% Plot outcomes
%==========================================================================
% https://www.england.nhs.uk/statistics/statistical-work-areas/critical-care-capacity/
% The NHS also maintains critical care beds for patients who are seriously
% ill and require constant support. Unlike most other categories of
% hospital bed, the total number of critical care beds has increased in
% recent years. In 2011/12 there were around 5,400 critical care beds, by
% 2019/20 this had risen to 5,900 (NHS England 2019b) (Figure 5). Of these,
% around 70 per cent are for use by adults and the remainder for children
% and infants
% https://www.telegraph.co.uk/global-health/science-and-disease/huge-regional-differences-intensive-care-bed-numbers-threaten/
%--------------------------------------------------------------------------
[t,n] = size(Y);
t     = (1:t)/7;
u     = 2000;

% factors and names
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% graphics
%--------------------------------------------------------------------------
subplot(3,2,1)
if n > 2
    plot(t,Y,t,u*t.^0,':m')
else
    plot(t,Y)
end
xlabel('time (weeks)'),ylabel('number of cases/day')
title('Rates (per day)','FontSize',16)
axis square, box off, set(gca,'XLim',[0, t(end)])
legend(str.outcome{1:n}), legend('boxoff')

subplot(3,2,2)
plot(t,cumsum(Y(:,1:min(end,2))))
xlabel('time (weeks)'),ylabel('number of cases')
title('Cumulative cases','FontSize',16)
set(gca,'YLim',[0, 50000]), box off

% marginal densities
%--------------------------------------------------------------------------
for i = 1:numel(X)
    subplot(3,2,2 + i)
    plot(t,X{i})
    xlabel('time (weeks)'),ylabel('probability')
    title(str.factors{i},'FontSize',16), set(gca,'YLim',[0,1]), set(gca,'XLim',[0, t(end)])
    axis square, box off, legend(str.factor{i}), legend('boxoff'), box off 
end

% add empirical data
%--------------------------------------------------------------------------
if nargin > 2
    t = (1:size(Z,1))/7;
    subplot(3,2,1), hold on, plot(t,Z,'.k'), hold off
    subplot(3,2,2), hold on, plot(t,cumsum(Z),'.k'), hold off
end
drawnow
