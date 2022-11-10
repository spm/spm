function spm_SARS_plot(Y,X,Z,U)
% Graphics for coronavirus simulations
% FORMAT spm_SARS_plot(Y,X,Z,U)
% Y      - expected timeseries (i.e., new depths and cases)
% X      - latent (marginal ensemble density) states
% Z      - optional empirical data (ordered as Y)
% U      - optional indices of outcomes
%
% This auxiliary routine plots the trajectory of outcome variables
% and underlying latent or hidden states, in the form of marginal densities
% over the four factors that constitute the SARS model. if empirical data
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

% plot prior responses
%==========================================================================
if nargin < 2
    
    % default outputs of interest (death, cases and hospitalisation)
    %----------------------------------------------------------------------
    if nargin < 1
        U = [1 2 3 16 27];
    else
        U = Y;
    end
    
    % number of age groups if not specified
    %----------------------------------------------------------------------
    nN      = 4;
    pE      = spm_SARS_priors(nN);
    M.T     = datenum(date) - datenum('01-01-2020','dd-mm-yyyy');
    [Y,X,Z] = spm_SARS_gen(pE,M,U);
    spm_SARS_plot(Y,X,Z,U)
    
    return
end

% defaults
%==========================================================================
[t,n,m]   = size(Y);
if nargin < 4, U = 1:n;  end
if nargin < 3 || isempty(Z), Z = zeros(0,n,m); end

% remove unpredicted data
%--------------------------------------------------------------------------
Z         = Z(:,1:min(n,end),:);

% deal with multiple groups or stratification
%--------------------------------------------------------------------------
if all(size(X) > 1) && m < 2
    for i = 1:size(X,2)
        j = [i,n];
        spm_SARS_plot(Y(:,j),X(:,i),Z(:,j),U)
        for j = 1:6, subplot(3,2,j), hold on, end
        if CHOLD, set(gca,'ColorOrderIndex',1); end          
    end
    return
end

% deal with multi-region outcomes
%==========================================================================
if m > 1
    for i = 1:m
        spm_SARS_plot(Y(:,:,i),X(:,i),Z(:,:,i),U)
        for j = 1:6, subplot(3,2,j), hold on, end
    end
    return
end

% plot hidden states and outcomes for this region or country
%==========================================================================

% factors and names
%--------------------------------------------------------------------------
[pE,pC,str] = spm_SARS_priors;

% graphics
%--------------------------------------------------------------------------
subplot(3,2,1), set(gca,'ColorOrderIndex',1);
t  = (1:t)/7;
p  = plot(t,Y);
nu = numel(U);
un = str.outcome(U);

xlabel('time (weeks)'),ylabel('number per day')
title('Rates (per day)','FontSize',16)
axis square, box off, set(gca,'XLim',[0, t(end)])
legend('off'), legend(p(1:nu),un), legend('boxoff')

subplot(3,2,2), set(gca,'ColorOrderIndex',1);
plot(t,cumsum(Y));
xlabel('time (weeks)'),ylabel('number of cases'), set(gca,'XLim',[0, t(end)])
title('Cumulative cases','FontSize',16), axis square, box off

% marginal densities
%--------------------------------------------------------------------------
k     = 4;
for i = 1:numel(X)
    
    k = k + 1;
    subplot(6,2,k), set(gca,'ColorOrderIndex',1);
    [d,j] = sort(max(X{i}));
    
    % remove redundant states
    %----------------------------------------------------------------------
    if i > 2
        j(end) = [];
    end
    
    if i == 1 || i == 3
        
        plot(t,X{i}(:,j(1:2))*100)
        ylabel('percent')
        title(str.factors{i},'FontSize',12), set(gca,'XLim',[0, t(end)])
        box off, legend(str.factor{i}(j(1:2))), legend('boxoff'), box off
        
        k = k + 1;
        subplot(6,2,k), set(gca,'ColorOrderIndex',1);
        j(1:2) = [];
        plot(t,X{i}(:,j)*100)
        ylabel('percent')
        title(str.factors{i},'FontSize',12), set(gca,'XLim',[0, t(end)])
        box off, legend(str.factor{i}(j)), legend('boxoff'), box off
        
    elseif i == 2
        
        plot(t,X{i}(:,j(1:4))*100)
        ylabel('percent')
        title(str.factors{i},'FontSize',12), set(gca,'XLim',[0, t(end)])
        box off, legend(str.factor{i}(j(1:4))), legend('boxoff'), box off
        
        k = k + 1;
        subplot(6,2,k), set(gca,'ColorOrderIndex',1);
        j(1:4) = [];
        plot(t,X{i}(:,j)*100)
        ylabel('percent')
        title(str.factors{i},'FontSize',12), set(gca,'XLim',[0, t(end)])
        box off, legend(str.factor{i}(j)), legend('boxoff'), box off
        
    else  
        
        plot(t,X{i}(:,j)*100)
        ylabel('percent')
        title(str.factors{i},'FontSize',12), set(gca,'XLim',[0, t(end)])
        box off, legend(str.factor{i}(j)), legend('boxoff'), box off
        
    end
    
end

% add empirical data
%--------------------------------------------------------------------------
t = (1:size(Z,1))/7;
try
    T   = numel(t);
    i   = 1:(T - 8);
    j   = (T - 7):T;
    subplot(3,2,2), set(gca,'ColorOrderIndex',1); hold on
    plot(t,cumsum(Z),'.'), hold off
    subplot(3,2,1), set(gca,'ColorOrderIndex',1); hold on
    plot(t(i),Z(i,:),'.'), plot(t(j),Z(j,:),'.c'), hold off
    legend('off'), legend(p(1:nu),un), legend('boxoff')
end
drawnow
