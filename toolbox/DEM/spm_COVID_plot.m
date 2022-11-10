function spm_COVID_plot(Y,X,Z,u,U)
% Graphics for coronavirus simulations
% FORMAT spm_COVID_plot(Y,X,Z)
% Y      - expected timeseries (i.e., new depths and cases)
% X      - latent (marginal ensemble density) states
% Z      - optional empirical data
% u      - optional bed capacity threshold
% U      - optional indices of outcomes
%
% This auxiliary routine plots the trajectory of outcome variables
% and underlying latent or hidden states, in the form of marginal densities
% over the four factors that constitute the COVID model. if empirical data
% are supplied, they will be superimposed.
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

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

% defaults
%--------------------------------------------------------------------------
[t,n,m]   = size(Y);
if nargin < 5, U = 1:n;  end
if nargin < 4, u = [];   end
if nargin < 3 || isempty(Z), Z = zeros(0,n,m); end

% remove unpredicted data
%--------------------------------------------------------------------------
Z         = Z(:,1:min(n,end),:);

% deal with multiple groups or stratification
%--------------------------------------------------------------------------
if all(size(X) > 1) && m < 2
    for i = 1:size(X,2)
        j = [i,n];
        spm_COVID_plot(Y(:,j),X(:,i),Z(:,j))
        for j = 1:6, subplot(3,2,j), hold on, end
        if CHOLD, set(gca,'ColorOrderIndex',1); end          
    end
    return
end

% deal with multi-region outcomes
%==========================================================================
if m > 1
    for i = 1:m
        spm_COVID_plot(Y(:,:,i),X(:,i),Z(:,:,i),u,U)
        for j = 1:6, subplot(3,2,j), hold on, end
    end
    return
end

% plot hidden states and outcomes for this region or country
%==========================================================================

% factors and names
%--------------------------------------------------------------------------
[pE,pC,str] = spm_COVID_priors;

% graphics
%--------------------------------------------------------------------------
subplot(3,2,1), if CHOLD, set(gca,'ColorOrderIndex',1); end
t    = (1:t)/7;
if ~isempty(u)
   p = plot(t,Y,t,u*t.^0,':m');
   legend(p,{str.outcome{U},'capacity'})
else
   p = plot(t,Y);
   legend(p,str.outcome{U})
end
xlabel('time (weeks)'),ylabel('number per day')
title('Rates (per day)','FontSize',16)
axis square, box off, set(gca,'XLim',[0, t(end)])
legend('boxoff')

subplot(3,2,2), if CHOLD, set(gca,'ColorOrderIndex',1); end
plot(t,cumsum(Y))
xlabel('time (weeks)'),ylabel('number of cases'), set(gca,'XLim',[0, t(end)])
title('Cumulative cases','FontSize',16), axis square, box off

% marginal densities
%--------------------------------------------------------------------------
for i = 1:numel(X)
    subplot(3,2,2 + i), if CHOLD, set(gca,'ColorOrderIndex',1); end
    plot(t,X{i}(:,2:end)*100)
    xlabel('time (weeks)'),ylabel('proportion (%)')
    title(str.factors{i},'FontSize',16), set(gca,'XLim',[0, t(end)])
    axis square, box off, legend(str.factor{i}(2:end)), legend('boxoff'), box off
end

% add empirical data
%--------------------------------------------------------------------------
t = (1:size(Z,1))/7;
try
    U   = U(ismember(U,1:size(Z,2)));
    subplot(3,2,2), if CHOLD, set(gca,'ColorOrderIndex',1); end; hold on, plot(t,cumsum(Z(:,U)),'.k'), hold off
    subplot(3,2,1), if CHOLD, set(gca,'ColorOrderIndex',1); end; hold on, plot(t,Z(:,U),'.'), hold off
end
drawnow
