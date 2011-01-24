function spm_plot_ci(varargin)
% plots mean and conditional confidence intervals
% FORMAT spm_plot_ci(t,E,C,j,s)
% FORMAT spm_plot_ci(t,E,C,j)
% FORMAT spm_plot_ci(t,E,C)
% FORMAT spm_plot_ci(E,C)
%
% t - domain
% E - expectation
% C - variance or covariance
% j - indices of spm_cat(E(:)) to plot
% s - string to specify plot type
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_plot_ci.m 4169 2011-01-24 18:34:20Z karl $

% unpack
%--------------------------------------------------------------------------
if nargin == 5
    t = varargin{1};
    E = varargin{2};
    C = varargin{3};
    j = varargin{4};
    s = varargin{5};
elseif nargin == 4
    t = varargin{1};
    E = varargin{2};
    C = varargin{3};
    j = varargin{4};
elseif nargin == 3
    t = varargin{1};
    E = varargin{2};
    C = varargin{3};
else
    E = varargin{1};
    C = varargin{2};
end

if iscell(E), E = spm_cat(E(:));       end
if ~exist('j','var'), j = 1:size(E,1); end
if ~exist('t','var'), t = 1:size(E,2); end
if ~exist('s','var'), s = '';          end

% order and length of sequence
%--------------------------------------------------------------------------
E     = E(j,:);
[n N] = size(E);

% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);
if iscell(C)
    for i = 1:N
        c(:,i) = ci*sqrt(diag(C{i}(j,j)));
    end
else
    if isvector(C)
        c = ci*sqrt(C);
    else
        C = diag(C);
        c = ci*sqrt(C(j,:));
    end
end

% set plot parameters
%--------------------------------------------------------------------------
switch get(gca,'NextPlot')
    case{lower('add')}
        col   = [1 1/4 1/4];
        width = .9;
    otherwise
        col   = [1 3/4 3/4];
        width = .8;
end

% conditional covariances
%--------------------------------------------------------------------------
if N > 1
    
    % time-series plot
    %======================================================================
    fill([t fliplr(t)],[full(E + c) fliplr(full(E - c))],...
         [1 1 1]*.8,'EdgeColor',[1 1 1]*.5),hold on
    plot(t,E,s)
    
elseif n == 2
    
    % plot in state-space
    %======================================================================    try,  C = C{1};  end
    [x y] = ellipsoid(E(1),E(2),0,c(1),c(2),0,32);
    fill(x',y',[1 1 1]*.8,'EdgeColor',[1 1 1]*.5),hold on
    plot(E(1,1),E(2,1),'.','MarkerSize',16)
    
else
    
    % bar
    %======================================================================
    
    % conditional means
    %----------------------------------------------------------------------
    bar(E,width,'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8), hold on
    box off
    set(gca,'XLim',[0 n + 1])
    
    % conditional variances
    %----------------------------------------------------------------------

    for k = 1:n
        line([k k], [-1 1]*c(k) + E(k),'LineWidth',4,'Color',col);
    end
end
hold off
drawnow
