function spm_DEM_qP(qP,pP)
% reports on conditional estimates of parameters
% FORMAT spm_DEM_qP(qP,pP)
%
% qP.P   - conditional expectations
% qP.V   - conditional variance
%
% pP     - optional priors
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_qP.m 7322 2018-05-31 09:47:15Z karl $


% unpack conditional covariances
%--------------------------------------------------------------------------
g     = length(qP.P);                                  % depth of hierarchy
ci    = spm_invNcdf(1 - 0.05);

% loop over levels
%--------------------------------------------------------------------------
Label = {};
col   = [1 3/4 3/4];
for i = 1:g
    
    % check for last level
    %----------------------------------------------------------------------
    if isempty(qP.P{i}), break, end

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
    try
        label = label(j);
    end
    try
        pi = spm_vec(pP.P{i});
        pi = pi(j);
    end
    np     = length(qi);
    if np
        
        % use current axes if P = P{1}
        %------------------------------------------------------------------
        if g > 1, subplot(g,1,i), end
                
        % conditional means
        %------------------------------------------------------------------
        bar(qi,'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8)
        title(sprintf('parameters - level %i',i),'FontSize',16);
        axis square
        box off
        set(gca,'XLim',[0 np + 1])

        % conditional variances
        %------------------------------------------------------------------
        for k = 1:np
            line([k k], [-1 1]*c(k) + qi(k),'LineWidth',4,'Color',col);
        end

        % prior or true means
        %------------------------------------------------------------------
        try
            hold on, bar(1:length(qi),pi,1/3), hold off
        end

        % labels
        %------------------------------------------------------------------
        for k = 1:length(label)
            text(k + 1/4,qi(k),label{k},'FontWeight','Bold','Color','r');
        end
        Label = [Label, label];
    end
end

% conditional (or prior) covariance 
%--------------------------------------------------------------------------
try
    if length(qP.C) == 1;
        return
    else
        i  = find(diag(qP.C));
    end
catch
    return
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


% plot evolution of hyperparameters if supplied
%==========================================================================
subplot(g,2,g + g)
try
    
    % confidence interval and expectations
    %----------------------------------------------------------------------
    ns = length(qP.p);
    t  = 1:ns;
    for i = 1:ns
        v(:,i) = sqrt(diag(qP.U*qP.c{i}*qP.U'));
    end
    c  = ci*v;
    p  = qP.U*spm_cat(qP.p);
    i  = find(any(v,2));
    c  = c(i,:);
    p  = p(i,:);
    
 
    % plot
    %----------------------------------------------------------------------
    hold on
    np    = size(p,1);
    for i = 1:np
        fill([t fliplr(t)],[(p(i,:) + c(i,:)) fliplr(p(i,:) - c(i,:))],...
            [1 1 1]*.8,'EdgeColor',[1 1 1]/2)
        plot(t,p(i,:))
    end
    set(gca,'XLim',[1 ns])
    title({'dynamics of parameters','(minus prior)'},'FontSize',16)
    xlabel('time')
    axis square
    hold off
 
catch
    
    % or correlations
    %----------------------------------------------------------------------
    imagesc(spm_cov2corr(qP.C(i,i)))
    title({'conditional correlations','among parameters'},'FontSize',16)
    axis square
    drawnow
    
end



