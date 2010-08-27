function spm_DEM_qH(qH,pH)
% reports on conditional estimates of hyperparameters
% FORMAT spm_DEM_qH(qH,pH);
%
% qH.h    - conditional estimate of log-precision (causes)
% qH.g    - conditional of log-precision (state)
% qH.V    - conditional variance (causes)
% qH.W    - conditional (states)
%
% qH.p    - time-dependent estimates from Laplace scheme
% qH.c    - time-dependent covariances
%
% pH      - option true log-precisions
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_qH.m 4052 2010-08-27 19:22:44Z karl $
 
% unpack conditional covariances
%--------------------------------------------------------------------------
try, qH = qH.qH; end
try, pH = pH.pH; end

% [Re]ML estimates - h
%==========================================================================
ci = spm_invNcdf(1 - 0.05);
h  = spm_vec(qH.h);
c  = spm_vec(qH.V);
c  = sqrt(c)*ci;
subplot(2,2,1)
bar(full(h),'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8)
title({'log-precision';'noise and causes'},'FontSize',16);
axis square
set(gca,'XLim',[0 length(c) + 1])
 
hlabel = {};
for i = 1:length(qH.h)
    for j = 1:length(qH.h{i})
        hlabel{end + 1} = sprintf('h:level %i',i);
    end
end
set(gca,'XTickLabel',hlabel)
 
% conditional covariances
%--------------------------------------------------------------------------
for i = 1:length(c)
    line([i i],[-1 1]*c(i) + h(i),'LineWidth',4,'Color','r');
end

% prior or true means
%--------------------------------------------------------------------------
try
    p     = spm_vec(pH.h);
    for i = 1:length(h)
        line([-1 1]/2 + i,[0 0] + p(i),'LineWidth',4,'Color','k');
    end
end

        
% [Re]ML estimates - g
%==========================================================================
h = spm_vec(qH.g);
c = spm_vec(qH.W);
c = sqrt(c)*spm_invNcdf(1 - 0.05);
if h
    subplot(2,2,2)
    bar(full(h),'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8)
    title({'log-precision';'states'},'FontSize',16);
    axis square
    set(gca,'XLim',[0 length(c) + 1])
    
    % conditional covariances
    %----------------------------------------------------------------------
    for i = 1:length(c)
        line([i i],[-1 1]*c(i) + h(i),'LineWidth',4,'Color','r');
    end
    
    % prior or true means
    %----------------------------------------------------------------------
    try
        p     = spm_vec(pH.g);
        for i = 1:length(h)
            line([-1 1]/2 + i,[0 0] + p(i),'LineWidth',4,'Color','k');
        end
    end
    
    glabel = {};
    for i = 1:length(qH.h)
        for j = 1:length(qH.h{i})
            glabel{end + 1} = sprintf('g:level %i',i);
        end
    end
    set(gca,'XTickLabel',glabel)
else
    glabel = {};
end

% conditional covariance and correlations
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc(qH.C)
title({'covariances among','hyperparameters'},'FontSize',16)
axis square
set(gca,'YTickLabel',{hlabel{:} glabel{:}},'YTick',[1:length(qH.C)])
 
% plot evolution of hyperparameters if supplied
%==========================================================================
subplot(2,2,4)
try
    
    % confidence interval and expectations
    %----------------------------------------------------------------------
    ns = length(qH.p);
    t  = 1:ns;
    for i = 1:ns
        v(:,i) = sqrt(diag(qH.c{i}));
    end
    c  = ci*v;
    h  = spm_cat(qH.p);
 
    % plot
    %----------------------------------------------------------------------
    hold on
    nh    = size(h,1);
    for i = 1:nh
        fill([t fliplr(t)],[(h(i,:) + c(i,:)) fliplr(h(i,:) - c(i,:))],...
            [1 1 1]*.8,'EdgeColor',[1 1 1]/2)
        plot(t,h(i,:))
    end
    set(gca,'XLim',[1 ns])
    title({'dynamics of','hyperparameters'},'FontSize',16)
    xlabel('time')
    axis square
    hold off
 
catch
    
    imagesc(spm_cov2corr(qH.C))
    title({'correlations among','hyperparameters'},'FontSize',16)
    axis square
    
end
