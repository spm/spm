function spm_CLIMATE_plot(Y,X,U)
% Graphics for coronavirus simulations
% FORMAT spm_CLIMATE_plot(Y,X,Z,U)
% Y      - expected timeseries (i.e., new depths and cases)
% X      - latent (marginal ensemble density) states
% U      - optional indices of outcomes
%
% This auxiliary routine plots the trajectory of outcome variables
% and underlying latent or hidden states, in the form of marginal densities
% over the four factors that constitute the SARS model. if empirical data
% are supplied, they will be superimposed.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_CLIMATE_plot.m 8151 2021-09-13 09:12:37Z karl $

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
global CHOLD, if isempty(CHOLD); CHOLD = 1; end

% plot prior responses
%==========================================================================

% defaults
%==========================================================================
[t,n,m]   = size(Y);
if nargin < 3, U = 1:3;  end

% plot hidden states and outcomes for this region or country
%==========================================================================

% factors and names
%--------------------------------------------------------------------------
[~,~,str] = spm_CLIMATE_priors;

% graphics
%--------------------------------------------------------------------------
subplot(3,1,1), set(gca,'ColorOrderIndex',1);
t  = (1:size(Y,1))/12;
plot(t,Y)

xlabel('time (years)'),ylabel('outcome')
title('Outcomes','FontSize',16)
box off, set(gca,'XLim',[0, t(end)])
legend(str.outcome(U)), legend('boxoff'), box off

% Latent states
%--------------------------------------------------------------------------
k = 1;
l = {};
subplot(3,2,2 + k)
for i = 1:size(X,2)
    
    plot(t,X(:,i)), hold on
    ylabel('state')
    title('Latent states','FontSize',12), set(gca,'XLim',[0, t(end)])
    l = {l{:}, str.states{i}};
    
    if ~rem(i,3)
        legend(l), legend('boxoff'), box off
        k = k + 1;
        l = {};
        if i < size(X,2)
            subplot(3,2,2 + k), set(gca,'ColorOrderIndex',1);
        end
    end
end

drawnow
