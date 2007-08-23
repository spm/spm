function spm_DEM_qU(qU,pU)
% displays conditional estimates of states (qU)
% FORMAT spm_DEM_qU(qU,pU);
%
% qU.v{i}    - causal states (V{1} = y = predicted response)
% qU.x{i}    - hidden states
% qU.e{i}    - prediction error
% qU.C{N}    - conditional covariance - [causal states] for N samples
% qU.S{N}    - conditional covariance - [hidden states] for N samples
%
% pU         - optional input for known states
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% unpack
%--------------------------------------------------------------------------
clf
V      = qU.v;
X      = qU.x;
E      = qU.z;
try
    C  = qU.C;
    S  = qU.S;
end
try
    pV = pU.v;
    pX = pU.x;
end
try
    pA = qU.a;
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
        plot(t,full(E{i}),'r:',t,full(V{i}))


        % conditional covariances
        %------------------------------------------------------------------
        if i > 1 & size(c,1)
            hold on
            j      = [1:size(V{i},1)];
            y      = ci*c(j,:);
            c(j,:) = [];
            fill([t fliplr(t)],[full(V{i} + y)' fliplr(full(V{i} - y)')],...
                 [1 1 1]*.8,'EdgeColor',[1 1 1]*.8)
            plot(t,full(E{i}),'r:',t,full(V{i}))
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
        try, 
            plot(t,pV{i},'linewidth',2,'color',[1 1 1]/2)
        end, hold on
        plot(t,full(E{i}(:,1:N)),'r:',t,full(V{i}))
        hold off
        set(gca,'XLim',[t(1) t(end)])
        a   = axis;

        % conditional covariances
        %------------------------------------------------------------------
        if i > 1 & size(c,1)
            hold on
            j      = [1:size(V{i},1)];
            y      = ci*c(j,:);
            c(j,:) = [];
            fill([t fliplr(t)],[full(V{i} + y) fliplr(full(V{i} - y))],...
                        [1 1 1]*.8,'EdgeColor',[1 1 1]*.8)
            plot(t,full(E{i}(:,1:N)),'r:',t,full(V{i}))
            hold off
        end

        % title, action, grid and true casues (if available)
        %------------------------------------------------------------------
        if i == 1
            title('predicted response and error');
        else
            title(sprintf('causal states - level %i',i));
            try, hold on
                plot(t,pV{i},'linewidth',2,'color',[1 1 1]/2)
            end, hold off
            try, hold on
                plot(t,pA{i - 1},'linewidth',1,'color',[1 0 0])
            end, hold off
        end
        xlabel('time {bins}')
        ylabel('states (a.u.)')
        grid on
        axis square
        axis(a)

        % hidden states
        %------------------------------------------------------------------
        try
 
            subplot(g,2,2*i)
            plot(t,full(X{i}))
            try, hold on
                plot(t,pX{i},'linewidth',2,'color',[1 1 1]/2)
            end, hold off
            set(gca,'XLim',[t(1) t(end)])
            a   = axis;
            
            if length(s)
                hold on
                j      = [1:size(X{i},1)];
                y      = ci*s(j,:);
                s(j,:) = [];
                fill([t fliplr(t)],[full(X{i} + y) fliplr(full(X{i} - y))],...
                        [1 1 1]*.8,'EdgeColor',[1 1 1]*.8)
                plot(t,full(X{i}))
                hold off

            end
            
            try, hold on
                plot(t,pX{i},'linewidth',2,'color',[1 1 1]/2)
            end, hold off
            
            % title and grid
            %--------------------------------------------------------------
            title('hidden states')
            xlabel('time {bins}')
            grid on
            axis square
            axis(a);
            
        catch
            delete(gca)
        end
    end
end
drawnow
