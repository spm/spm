function spm_DEM_qP(qP,pP)
% reports on conditional estimates of parameters
% FORMAT spm_DEM_qP(qP,pP)
%
% qP.P    - conditional expectations
% qP.Pi   - conditional expectations - hierarchical form
% qP.C    - conditional covariance
% 
% pP      - optional priors
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_qP.m 2029 2008-09-02 18:26:23Z karl $


% time-series specification
%--------------------------------------------------------------------------
clf
g     = length(qP.P);                                  % depth of hierarchy

% unpack conditional covariances
%--------------------------------------------------------------------------
ci     = spm_invNcdf(1 - 0.05);
c      = [];
p      = [];
try
    c  = full(sqrt(diag(qP.C)));
end


% loop over levels
%--------------------------------------------------------------------------
for i = 1:(g - 1)
    
    % conditional expectations
    %----------------------------------------------------------------------
    qi    = spm_vec(qP.P{i});
    dk    = 0;
    try
        qi = [qi spm_vec(pP.P{i})];
        dk = -1/8;
    end
    j     = length(qi);
    if j
        subplot(g,1,i)
        bar(qi)        
        title(sprintf('parameters - level %i',i));
        grid on
        axis square
        set(gca,'XLim',[0 j + 1])
        
        % conditional covariances
        %------------------------------------------------------------------
        try
            for k = 1:j
                line([k k] + dk,[-1 1]*ci*c(k) + qi(k),...
                    'LineWidth',4,'Color','r');
            end
            c(1:j) = [];
        end
    end
end

% conditional covariance and correlations
%--------------------------------------------------------------------------
try
    if length(qP.C) == 1;
        return
    end
catch
    return
end

subplot(g,2,g + g - 1)
if exist('pC','var')
    imagesc(spm_cov2corr(pC))
    title({'prior correlations','among parameters'})
else
    imagesc(qP.C)
    title({'conditional covariances','among parameters'})
end
axis square

subplot(g,2,g + g)
imagesc(spm_cov2corr(qP.C))
title({'conditional correlations','among parameters'})
axis square
drawnow
