function spm_DEM_qH(qH)
% reports on conditional estimates of hyperparameters
% FORMAT spm_DEM_qH(qH);
%
% qH.h    - ReML estimate
% qH.hi   - ReML estimate - hierarchical form
% qH.C    - covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_qH.m 3058 2009-04-09 18:17:53Z karl $

% unpack conditional covariances
%--------------------------------------------------------------------------
try
    qH = qH.qH;
end
c     = [];
try
    c = sqrt(diag(qH.C))*spm_invNcdf(1 - 0.05);
end

% [Re]ML estimates
%--------------------------------------------------------------------------
h = spm_vec(qH.h);
j = length(qH.h);
subplot(2,1,1)
bar(full(h),'c')
title('hyperparameters');
grid on
axis square
set(gca,'XLim',[0 j + 1])

% conditional covariances
%--------------------------------------------------------------------------
try
    for i = 1:j
        line([i i],[-1 1]*c(i) + h(i),...
            'LineWidth',4,...
            'Color','r');
    end
end

% conditional covariance and correlations
%--------------------------------------------------------------------------
try
    if length(qH.C) > 1
        subplot(2,2,3)
        imagesc(qH.C)
        title({'covariances among','ln(hyperparameters)'})
        axis square

        subplot(2,2,4)
        imagesc(spm_cov2corr(qH.C))
        title({'correlations among','ln(hyperparameters)'})
        axis square
    end
end
