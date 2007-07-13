function spm_DEM_qU(varargin)
% displays conditional estimates of states (qU)
% FORMAT spm_DEM_qU(qU);
% FORMAT spm_DEM_qU(v{i},x{i},e{i});
%
% qU.v{i}    - causal states (V{1} = y = response)
% qU.x{i}    - hidden states
% qU.e{i}    - prediction error
% qU.C{N}    - conditional covariance - [causal states] for N samples
% qU.S{N}    - conditional covariance - [hidden states] for N samples
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% unpack
%--------------------------------------------------------------------------
if isstruct(varargin{1})
    try
        varargin{1} = varargin{1}.qU;
    end
    V     = varargin{1}.v;
    X     = varargin{1}.x;
    E     = varargin{1}.z;
    try
        C = varargin{1}.C;
        S = varargin{1}.S;
    end
elseif nargin == 3
    V     = varargin{1};
    X     = varargin{2};
    E     = varargin{3};
end

% time-series specification
%--------------------------------------------------------------------------
g     = length(V);          % order of hierarchy
N     = size(V{1},2);       % length of data sequence
dt    = 1;                  % time step
t     = [1:N]*dt;           % time


% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);
s     = [];
c     = [];
try
    for i = 1:N
        c = [c sqrt(diag(C{i}))];
        s = [s sqrt(diag(S{i}))];
    end
end

% loop over levels
%--------------------------------------------------------------------------
for i = 1:g

    if N == 1

        % causal states and error - single observation
        %------------------------------------------------------------------
        subplot(g,2,2*i - 1)
        t = 1:size(V{i},1);
        plot(t,full(V{i}),t,full(E{i}),':')


        % conditional covariances
        %------------------------------------------------------------------
        if i > 1 & size(c,1)
            hold on
            j      = [1:size(V{i},1)];
            y      = ci*c(j,:);
            c(j,:) = [];
            plot(t,full(V{i} + y),'-.',t,full(V{i} - y),'-.')
            hold off
        end

        % title and grid
        %------------------------------------------------------------------
        title(sprintf('causal states - level %i',i));
        xlabel('states')
        grid on
        axis square
        set(gca,'XLim',[t(1) t(end)])

    else

        % causal states and error - time series
        %------------------------------------------------------------------
        subplot(g,2,2*i - 1)
        plot(t,full(V{i}),t,full(E{i}(:,1:N)),':')
        title(sprintf('causal states - level %i',i));
        xlabel('time {bins}')
        ylabel('states (a.u.)')
        grid on
        axis square
        set(gca,'XLim',[t(1) t(end)])

        % conditional covariances
        %------------------------------------------------------------------
        if i > 1 & size(c,1)
            hold on
            j      = [1:size(V{i},1)];
            y      = ci*c(j,:);
            c(j,:) = [];
            plot(t,full(V{i} + y),'-.',t,full(V{i} - y),'-.')
            hold off
        end

        % title and grid
        %------------------------------------------------------------------
        title(sprintf('causal states - level %i',i));
        xlabel('time {bins}')
        ylabel('states (a.u.)')
        grid on
        axis square
        set(gca,'XLim',[t(1) t(end)])
  

        % hidden states
        %------------------------------------------------------------------
        try
 
            subplot(g,2,2*i)
            plot(t,full(X{i}))
            title('hidden states')
            xlabel('time {bins}')
            grid on
            axis square
            set(gca,'XLim',[t(1) t(end)])
            a   = axis;
            
            if length(s)
                hold on
                j      = [1:size(X{i},1)];
                y      = ci*s(j,:);
                s(j,:) = [];
                plot(t,full(X{i} + y),'-.',t,full(X{i} - y),'-.')
                hold off
                axis(a);
            end
        catch
            delete(gca)
        end
    end
end
drawnow
