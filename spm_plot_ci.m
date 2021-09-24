function spm_plot_ci(E,C,x,j,s)
% Plot mean and conditional confidence intervals
% FORMAT spm_plot_ci(E,C,x,j,s)
% E - expectation (structure or array)
% C - variance or covariance (structure or array)
% x - domain
% j - rows of E to plot
% s - string to specify plot type:e.g. '--r' or 'exp', 'log' etc
%
% If E is a row vector with two elements, confidence regions will be
% plotted; otherwise, bar charts with confidence intervals are provided.
%__________________________________________________________________________
% Copyright (C) 2008-2021 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_plot_ci.m 8154 2021-09-24 11:25:10Z karl $


% get axis
%--------------------------------------------------------------------------
ax   = gca;
col  = get(ax,'ColorOrder');
coli = get(ax,'ColorOrderIndex');
coll = col(coli,:);
colf = erf(coll + 1);

% confidence region (CR) plotting
%--------------------------------------------------------------------------
if size(E,1) == 1 && size(E,2) == 2
    E  = E';
    CR = true;
else
    CR = false;
end
    

% unpack expectations into a matrix
%--------------------------------------------------------------------------
if isstruct(E),       E = spm_vec(E);    end
if iscell(E),         E = spm_vec(E);    end

if ~exist('x','var'), x = 1:size(E,2);   end
if ~exist('j','var'), j = 1:size(E,1);   end
if ~exist('s','var'), s = '';            end

if isempty(x),        x = 1:size(E,2);   end
if isempty(j),        j = 1:size(E,1);   end


% order and length of sequence
%--------------------------------------------------------------------------
O     = E;
E     = E(j,:);
[n,N] = size(E);

% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);               % confidence interval
gr    = 0.9;                                 % grey level
if iscell(C)
    
    % try cell array of covariances (from spm_DEM amd spm_LAP)
    %----------------------------------------------------------------------
    try
        for i = 1:N
            c(:,i) = ci*sqrt(diag(C{i}(j,j)));
        end
    catch
        
        % try cell array of variances
        %------------------------------------------------------------------
        c = ci*sqrt(spm_unvec(spm_vec(C),O));
        c = c(j,:);
    end
    
elseif isstruct(C)
    
    % try structure of variances
    %----------------------------------------------------------------------
    c = ci*sqrt(spm_unvec(spm_vec(C),O));
    c = c(j,:);
    
elseif isnumeric(C)
    
    % try matrix of variances
    %----------------------------------------------------------------------
    if all(size(C) == size(O))
        c = ci*sqrt(C(j,:));
    elseif all(size(C') == size(O))
        c = ci*sqrt(C(:,j));
        c = c(:)';
    else
        
        % try covariance matrix
        %------------------------------------------------------------------
        C = diag(C);
        c = ci*sqrt(C);
        c = c(:)';
    end
    
end


% set plot parameters
%--------------------------------------------------------------------------
switch lower(get(ax,'NextPlot'))
    case 'add'
        col   = [1 1/4 1/4];
        width = .9;
    otherwise
        col   = [1 3/4 3/4];
        width = .8;
end

% plot elliptical confidence region
%--------------------------------------------------------------------------
if CR
    [x,y] = ellipsoid(E(1),E(2),1,c(1),c(2),0,32);
    fill(x(16,:)',y(16,:)',[1 1 1]*gr,'EdgeColor',[1 1 1]*.5,'Parent',ax);
    hold(ax,'on');
    plot(ax,E(1),E(2),'.','MarkerSize',16);
    hold(ax,'off'); drawnow
    return
end


% plot bar chart
%--------------------------------------------------------------------------
if N >= 8
    
    % time-series plot
    %======================================================================
    x  = x(:)';
    if strcmpi(s,'exp')
        fill([x fliplr(x)],exp([full(E + c) fliplr(full(E - c))]),...
            colf,'EdgeColor','none','Parent',ax,'FaceAlpha', 0.4);
        hold(ax,'on');
        plot(x,exp(E),'Color',coll);
        set(ax,'ColorOrderIndex',coli + 1);
        
    elseif strcmpi(s,'log')
        fill([x fliplr(x)],log(abs([full(E + c) fliplr(full(E - c))])),...
            colf,'EdgeColor','none','Parent',ax,'FaceAlpha', 0.4);
        hold(ax,'on');
        plot(x,log(abs(E)),'Color',coll);
        set(ax,'ColorOrderIndex',coli + 1);
        
        
    else
        fill([x fliplr(x)],[full(E + c) fliplr(full(E - c))],...
            colf,'EdgeColor','none','Parent',ax,'FaceAlpha', 0.4);
        hold(ax,'on');
        plot(ax,x,E,s,'Color',coll);
        set(ax,'ColorOrderIndex',coli + 1);
        
    end
    
else
    
    % bar
    %======================================================================
    if N == 1
        
        if strcmpi(s,'exp')
            
            % conditional means
            %--------------------------------------------------------------
            bar(ax,exp(E),width,'Edgecolor',colf,'FaceColor',colf);
            hold(ax,'on');
            
            % conditional variances
            %--------------------------------------------------------------
            for k = 1:n
                line([k k],exp([-1 1]*c(k) + E(k)),...
                    'LineWidth',4,'Color',col,'Parent',ax);
            end
            
        elseif strcmpi(s,'log')
            
            % conditional means
            %--------------------------------------------------------------
            bar(ax,log(abs(E)),width,'EdgeColor',colf,'FaceColor',colf);
            hold(ax,'on');
            
            % conditional variances
            %--------------------------------------------------------------
            for k = 1:n
                line([k k],log(abs([-1 1]*c(k) + E(k))),...
                    'LineWidth',4,'Color',col,'Parent',ax);
            end
            
        else
                        
            if n > 1
                
                % conditional means
                %----------------------------------------------------------
                bar(ax,E,width,'EdgeColor',colf,'FaceColor',colf);
                hold(ax,'on');
                
            else
                % conditional means
                %----------------------------------------------------------
                bar(ax,E,'EdgeColor',colf,'FaceColor',colf);
                hold(ax,'on');
                
            end
            
            % conditional variances
            %--------------------------------------------------------------
            for k = 1:n
                line([k k],[-1 1]*c(k) + E(k),...
                    'LineWidth',4,'Color',col,'Parent',ax);
            end
            
        end
        
        box(ax,'off');
        set(ax,'XLim',[0 n + 1]);
        
    else
        
        if strcmpi(s,'exp')
            
            % conditional means (exponential)
            %--------------------------------------------------------------
            h = bar(ax,exp(E)'); hold(ax,'on');
            
            % conditional variances
            %--------------------------------------------------------------
            for m = 1:n
                if ~isempty(get(h(m),'Children'))
                    x = mean(get(get(h(m),'Children'),'XData'));
                else
                    x = get(h(m),'XEndPoints');
                end
                for k = 1:N
                    line([x(k) x(k)],exp([-1 1]*c(m,k) + E(m,k)),...
                        'LineWidth',1,'Color',col,'Parent',ax);
                end
            end
            
        else
            
            % conditional means
            %--------------------------------------------------------------
            h = bar(ax,E); hold(ax,'on');
            
            % conditional variances
            %--------------------------------------------------------------
            for m = 1:N
                if ~isempty(get(h(m),'Children'))
                    x = mean(get(get(h(m),'Children'),'Xdata'));
                else
                    x = get(h(m),'XEndPoints');
                end
                for k = 1:n
                    line([x(k) x(k)],[-1 1]*c(k,m) + E(k,m),...
                        'LineWidth',4,'Color',col,'Parent',ax);
                end
            end
            
        end
    end
    
end
hold(ax,'off');
drawnow
