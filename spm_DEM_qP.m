function spm_DEM_qP(qP,pP)
% reports on conditional estimates of parameters
% FORMAT spm_DEM_qP(qP,pP)
%
% qP.P    - conditional expectations
% qP.V   - conditional variance
%
% pP      - optional priors
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_qP.m 3695 2010-01-22 14:18:14Z karl $


% time-series specification
%--------------------------------------------------------------------------
clf
g     = length(qP.P);                                  % depth of hierarchy

% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);

% loop over levels
%--------------------------------------------------------------------------
Label = {};
for i = 1:(g - 1)

    % get lablels
    %----------------------------------------------------------------------
    label = {};
    if isstruct(qP.P{i})
        names = fieldnames(qP.P{i});
        for j = 1:length(names)
            for k = 1:length(spm_vec(getfield(qP.P{i},names{j})))
                label{end + 1} = names{j};
            end
        end
    end

    % conditional expectations (with priors if specified)
    %----------------------------------------------------------------------
    qi     = spm_vec(qP.P{i});
    c      = sqrt(spm_vec(qP.V{i}))*ci;
    j      = find(c);
    qi     = qi(j);
    c      = c(j);
    label  = label(j);
    np     = length(qi);
    try
        pi = spm_vec(pP.P{i});
        pi = pi(j);
    end
    
    if np
        
        % conditional means
        %------------------------------------------------------------------
        subplot(g,1,i)
        bar(qi)
        title(sprintf('parameters - level %i',i));
        axis square
        set(gca,'XLim',[0 np + 1])

        % conditional variances
        %------------------------------------------------------------------
        for k = 1:np
            line([k k], [-1 1]*c(k) + qi(k),'LineWidth',4,'Color','r');
        end

        % prior or true means
        %------------------------------------------------------------------
        try
        for k = 1:np
            line([-1 1]/2 + k,[0 0] + pi(k),'LineWidth',4,'Color','b');
        end
        end

        % labels
        %------------------------------------------------------------------
        for k = 1:length(label)
            text(k + 1/4,qi(k),label{k},'FontSize',12,'FontWeight','Bold','Color','g');
        end
        Label = {Label{:}, label{:}};
    end
end

% conditional (or prior) covariance 
%--------------------------------------------------------------------------
if length(qP.C) == 1;
    return
else
    i  = find(diag(qP.C));
end

subplot(g,2,g + g - 1)
if exist('pC','var')
    imagesc(spm_cov2corr(pC(i,i)))
    title({'prior correlations','among parameters'},'FontSize',16)
else
    imagesc(qP.C(i,i))
    title({'conditional covariances','among parameters'},'FontSize',16)
end
if ~isempty(Label)
    set(gca,'YTickLabel',Label,'YTick',[1:length(Label)])
end
axis square


% and correlations
%--------------------------------------------------------------------------
subplot(g,2,g + g)
imagesc(spm_cov2corr(qP.C(i,i)))
title({'conditional correlations','among parameters'},'FontSize',16)
axis square
drawnow
