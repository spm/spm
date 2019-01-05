function DEM_FEP_Least_Action
%--------------------------------------------------------------------------
% This routine uses a Lorenz system to show that the most likely autonomous
% path (or deterministic path) is uniquely identified by its initial and
% final states. In other words, if we knew the end state of an autonomous
% trajectory, then we would implicitly know the path taken from the initial
% particular state, even if we did not know the external states. This point
% is demonstrated using numerical analyses of the Lorenz system; treating
% the first state and an active state and the third as an external state.
% In this example, 1024 solutions are obtained from the same initial
% particular (i.e., sensory and active) states but sampling from a Gaussian
% distribution over external states. The ensuing trajectories over 128 time
% bins of 1/128 seconds are shown in the left panels. The sample
% distribution over active states is shown as a (scaled) histogram along
% the x-axis. Paths that end within 1/8 of an arbitrary active state (here,
% ?? = -4) are shown in red. The corresponding autonomous (i.e., active)
% paths are shown as a function of time in the right panels. one can repeat
% this analysis for different levels of random fluctuations; e.g.,log
% precisions of 2 and 16. The key thing to observe is that as the amplitude
% of random fluctuations decreases (i.e., its precision increases) the
% paths that begin and end in the same place collapse to a single
% trajectory of least action. This is the most likely or deterministic
% path. Clearly, this behaviour rests upon a diffeomorphic mapping between
% the initial and final states: for example, a final active state of -8 has
% the least two paths of least action (xT in the code below).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_FEP_Least_Action.m 7512 2019-01-05 21:15:16Z karl $

% generative model
%==========================================================================                       % switch for demo
spm_figure('GetWin','DEM'); clf

% flow
%--------------------------------------------------------------------------
dt    = 1/128; 
G.f   = @(x,v,P,G) v(:) + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x*dt;

% set up
%--------------------------------------------------------------------------
T     = 128;                        % length of trajectory
N     = 1024;                       % number of paths
g     = exp(-8);                    % amplitude of random fluctuations
s     = exp(2);                     % amplitude deviations
P     = [10; -8/3; 28];             % Rayleigh parameter
x0    = [1; 1; 25];
for k = 1:N
    
    % random fluctuations (and action)
    %----------------------------------------------------------------------
    U.u      = randn(T,3)*g;
    A(k)     = (1/2)*sum(U.u(:).^2)/(g^2)/T;
    
    % random deviations in external states (and probability)
    %----------------------------------------------------------------------
    dx       = randn*s;
    A(k)     = A(k) + (1/2)*sum(dx.^2)/(s^2);
        
    % integrate timeseries with random initial hidden state
    %----------------------------------------------------------------------
    G.x      = x0;
    G.x(3)   = x0(3) + dx;
    x(:,:,k) = spm_int_L(P,G,U);
    
end

% plot ensemble densities
%--------------------------------------------------------------------------
subplot(2,2,1),cla, hold on
td    = fix(linspace(1,N,32));
for t = td
    plot(x(:,1,t),x(:,3,t),':','MarkerSize',1)
    plot(x(T,1,t),x(T,3,t),'r.','MarkerSize',8)
end
axis square, axis([-20 20 0 60])
title('Paths from initial state','FontSize',16)
xlabel('active state'),ylabel('external state'),drawnow

% find trajectories that start and end at the same place
%==========================================================================


% get surprisal of final (active) state using sample density
%--------------------------------------------------------------------------
nb    = 32;
xN    = squeeze(x(end,1,:));
[n,a] = hist(xN,nb);
n     = 2*nb*n/sum(n);
bar(a,n,1)

% identify and plot trajectories with the same endpoints
%--------------------------------------------------------------------------
xT    = -4;
k     = find(abs(xN - xT) < 1/8);
for t = k(:)'
    plot(x(:,1,t),x(:,3,t),'r')
end
plot([x0(1),x0(1)],[0 48],'r-.')
plot([xT   ,xT   ],[0 48],'r-.')

% paths as a function of time
%--------------------------------------------------------------------------
subplot(2,2,2), hold on
pst   = (1:T)*dt;
for t = k(:)'
    plot(pst,x(:,1,t),'r')
end
title('Autonomous paths','FontSize',16)
xlabel('time (seconds)'),ylabel('active state'),axis square


return

% NB: numerical analysis of the action and surprisal
%--------------------------------------------------------------------------
for i = 1:N
    [d,j] = min(abs(xn - xN(i)));
    L(i)  = -log(n(j));
end
L      = log(spm_softmax(L(:)));
A      = log(spm_softmax(A(:)));

subplot(2,2,2)
plot(L,A,'.',L,L,':')
title('Action and surprisal','FontSize',16)
xlabel('Surprisal of final state'),ylabel('Action of path'),axis square

