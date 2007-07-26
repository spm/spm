function spm_DEM_qP(qP,M)
% reports on conditional estimates of parameters
% FORMAT spm_DEM_qP(qP,[M]);
%
% qP.P    - conditional expectations
% qP.Pi   - conditional expectations - hierarchical form
% qP.C    - conditional covariance
%
% M       - model structure for priors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% unpack
%--------------------------------------------------------------------------
try
    M  = qP.M;
    qP = qP.qP;
end


% time-series specification
%--------------------------------------------------------------------------
g     = length(qP.P);                                  % depth of hierarchy

% unpack conditional covariances
%--------------------------------------------------------------------------
ci     = spm_invNcdf(1 - 0.05);
c      = [];
p      = [];
try
    c  = full(sqrt(diag(qP.C)));
end
try
    pC = spm_cat(diag({M.pC}));
    p  = full(sqrt(diag(pC)));
end

% loop over levels
%--------------------------------------------------------------------------
for i = 1:(g - 1)
    
    % conditional expectations
    %----------------------------------------------------------------------
    qi    = spm_vec(qP.P{i});
    try
        pi = spm_vec(M.pE{i});
    end
    j     = length(qi);
    if j
        subplot(g,1,i)
        bar(qi,'c')        
        title(sprintf('parameters - level %i',i));
        grid on
        axis square
        set(gca,'XLim',[0 j + 1])
        
        % conditional covariances
        %------------------------------------------------------------------
        try
            for k = 1:j
                line([k k],[-1 1]*ci*c(k) + qi(k),...
                    'LineWidth',4,'Color',[0 0 0] + 2/8);
            end
            c(1:j) = [];
        end
        
        % prior covariances
        %------------------------------------------------------------------
        try
            for k = 1:j
                line([k k] + 1/2,[-1 1]*ci*p(k) + pi(k),...
                    'LineWidth',8,'Color',[0 0 0] + 6/8);
            end
            p(1:j) = [];
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
