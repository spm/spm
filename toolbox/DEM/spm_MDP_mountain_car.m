function spm_MDP_mountain_car(X,V,T)
% Demo for Discrete Markov Decision process (planning)
% FORMAT spm_MDP_mountain_car(X,V,T))
% X    - initial and goal position
% V    - initial and goal velocity
% T    - number of time-steps
%
% This routine uses a Markov decisions process formulation of the mountain
% car problem to illustrate prospective free energy minimization under a
% variational Bayesian learning scheme. The key notion here is that the
% agent represents future states and action (in a pullback sense), where it
% has strong prior beliefs about future states. The intervening states and
% actions are optimized with respect to current sensory data to provide
% predictions about the next sensory state, which action fulfils. The
% result is a planned trajectory through state space that realizes prior
% beliefs in a prospective sense.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_mountain_car.m 6061 2014-06-21 09:02:42Z karl $
 
% set up and preliminaries
%==========================================================================
rng('default');
global eta
eta  = 16;
u    = [-2 -1 0 1 2];                       % levels of control
Nu   = length(u);                           % number of control/actions
Nx   = 32;                                  % number of cells per dimension
Ns   = Nx*Nx;                               % number of states
dt   = 2;
 
% first and final states
%--------------------------------------------------------------------------
try, T;  catch, T = 16;    end
try, X;  catch, X = [0 1]; end
try, V;  catch, V = [0 0]; end
 
X(1)  = X(1);                               % first position
V(1)  = V(1);                               % first velocity
X(T)  = X(2);                               % final position
V(T)  = V(2);                               % final velocity
 
% smoothing matrix to simulate noise
%--------------------------------------------------------------------------
K     = toeplitz(sparse(1,[1 2],[1 1/2],1,Nx));
K     = K + K';
K     = K*diag(1./sum(K,1));

% set-up transition matrix for mountain car problem
%--------------------------------------------------------------------------
x     = linspace(-2,2,Nx);
v     = linspace(-3,3,Nx);
for k = 1:Nu
    P{k}  = sparse(Ns,Ns);
    for i = 1:length(x)
        for j = 1:length(v)
            
            % change in state
            %--------------------------------------------------------------
            ds = spm_fx_mountaincar([x(i);v(j)],0,u(k),[]);
            
            % transition probabilities - position
            %--------------------------------------------------------------
            dx = x - (x(i) + ds(1));
            ii = find(dx > 0,1);
            if ii == 1,
                px = sparse(1,1,1,Nx,1);
            elseif isempty(ii)
                px = sparse(Nx,1,1,Nx,1);
            else
                ii = [ii - 1,ii];
                px = pinv([x(ii); 1 1])*[(x(i) + ds(1)); 1];
                px = sparse(ii,1,px,Nx,1);
            end
            
            % transition probabilities - velocity
            %--------------------------------------------------------------
            dv = v - (v(j) + ds(2));
            ii = find(dv > 0,1);
            if ii == 1,
                pv = sparse(1,1,1,Nx,1);
            elseif isempty(ii)
                pv = sparse(Nx,1,1,Nx,1);
            else
                ii = [ii - 1,ii];
                pv = pinv([v(ii); 1 1])*[(v(j) + ds(2)); 1];
                pv = sparse(ii,1,pv,Nx,1);
            end
            
            % place in P
            %--------------------------------------------------------------
            pv     = K*pv;
            p      = px*pv';
            P{k}(:,(j - 1)*Nx + i) = p(:);
            
        end
    end
end
 
% make into discrete transition probabilities
%--------------------------------------------------------------------------
for i = 1:Nu
    P{i} = P{i}^dt;
end
 
% Generative MDP model – transition probabilities (empirical priors)
%--------------------------------------------------------------------------
for i = 1:Nu
    B{i}   = P{i}';                                 % transition priors
    P{i}   = P{i}*diag(1./max(sum(P{i}),eps));      % push-forward
    B{i}   = B{i}*diag(1./max(sum(B{i}),eps));      % pull-back
end
 
 
% graphics (transition probabilities)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf
 
[d i] = min(abs(x - X(2)));
[d j] = min(abs(v - V(2)));
j     = (j - 1)*Nx + i;
s     = sparse(j,1,1,Ns,1);
k     = [0 4 8];
for i = 1:Nu
    for j = 1:3
        subplot(Nu,3,(i - 1)*3 + j)
        imagesc(x,v,reshape(B{i}^k(j)*s,Nx,Nx)'), axis xy square
        xlabel('position','FontSize',12)
        ylabel('velocity','FontSize',12)
    end
end
drawnow;
 
% initial and final states
%==========================================================================
spm_figure('GetWin','Figure 2');clf
 
 
% first state (s)
%--------------------------------------------------------------------------
[d i] = min(abs(x - X(1)));
[d j] = min(abs(v - V(1)));
X(1)  = x(i);
V(1)  = v(j);
S     = sparse(i,j,1,Nx,Nx);

 
% final state (c)
%--------------------------------------------------------------------------
[d i] = min(abs(x - X(T)));
[d j] = min(abs(v - V(T)));
X(T)  = x(i);
V(T)  = v(j);
C     = sparse(i,j,1,Nx,Nx);

 
% (uniform) cost over control (d)
%--------------------------------------------------------------------------
D     = ones(Nu,1);

% solve
%==========================================================================
MDP.T = T;                         % process depth (the horizon)
MDP.S = S;                         % initial state
MDP.B = P;                         % transition probabilities (priors)
MDP.C = C;                         % terminal cost probabilities (priors)
MDP.D = D;                         % control probabilities (priors)

MDP.lambda = 64;
MDP.plot   = 1;

[Q,R,S]  = spm_MDP(MDP);
 
% set up state space graphically
%--------------------------------------------------------------------------
subplot(2,1,1)
for i = 1:length(x);
    plot(x(i),v,'k','color',[1 1 1]/2), hold on
    axis([-2 2 -2 2])
end
 
% and plot expected and realised trajectories
%==========================================================================
Sx(1) = X(1);
Sv(1) = V(1);
for k = 1:(T - 1)
    
    for t = 1:T
        
        % conditional expectations
        %--------------------------------------------------------------
        q      = reshape(Q(:,k,t),Nx,Nx);
        qx     = sum(q,2);
        qv     = sum(q,1);
        X(t)   = x*qx /sum(qx);
        V(t)   = v*qv'/sum(qv);
        
    end
    
    % plot expectation at this at time k
    %------------------------------------------------------------------
    subplot(2,1,1)
    plot(X,V,':r'), hold on
    axis([-2 2 -2 2])
    
    % current position
    %----------------------------------------------------------------------
    [i j]     = find(reshape(S(:,k + 1),Nx,Nx));
    Sx(k + 1) = x(i);
    Sv(k + 1) = v(j);
    
end
 
subplot(2,1,1)
plot(Sx,Sv,'ok','LineWidth',4)
plot(Sx,Sv,'k','LineWidth',2)
title('Expected and actual trajectory','FontSize',16)
xlabel('position','FontSize',12)
ylabel('velocity','FontSize',12)
axis([-2 2 -2 2])
 


