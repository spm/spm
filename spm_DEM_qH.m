function spm_DEM_qH(qH)
% reports on conditional estimates of hyperparameters
% FORMAT spm_DEM_qH(qH);
%
% qH.h    - ReML estimate of log precision (causes)
% qH.g    - ReML estimate of log precision (state)
% qH.V    - conditional variance (causes)
% qH.W    - conditional (states)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_qH.m 3655 2009-12-23 20:15:34Z karl $

% unpack conditional covariances
%--------------------------------------------------------------------------
try, qH = qH.qH; end

% [Re]ML estimates - h
%--------------------------------------------------------------------------
h = spm_vec(qH.h);
c = spm_vec(qH.V);
c = sqrt(c*spm_invNcdf(1 - 0.05));
subplot(2,2,1)
bar(full(h),'c')
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

% [Re]ML estimates - g
%--------------------------------------------------------------------------
h = spm_vec(qH.g);
c = spm_vec(qH.W);
c = sqrt(c*spm_invNcdf(1 - 0.05));
subplot(2,2,2)
bar(full(h),'c')
title({'log-precision';'states'},'FontSize',16);
axis square
set(gca,'XLim',[0 length(c) + 1])

% conditional covariances
%--------------------------------------------------------------------------
for i = 1:length(c)
    line([i i],[-1 1]*c(i) + h(i),'LineWidth',4,'Color','r');
end

glabel = {};
for i = 1:length(qH.h)
    for j = 1:length(qH.h{i})
        glabel{end + 1} = sprintf('g:level %i',i);
    end
end
set(gca,'XTickLabel',glabel)

% conditional covariance and correlations
%--------------------------------------------------------------------------
if length(qH.C) > 1
    subplot(2,2,3)
    imagesc(qH.C)
    title({'covariances among','ln(hyperparameters)'},'FontSize',16)
    axis square
    set(gca,'YTickLabel',{hlabel{:} glabel{:}},'YTick',[1:length(qH.C)])

    subplot(2,2,4)
    imagesc(spm_cov2corr(qH.C))
    title({'correlations among','ln(hyperparameters)'},'FontSize',16)
    axis square
end

