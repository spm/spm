function spm_MDP_offer
% Demo for active inference with limited offer game
%__________________________________________________________________________
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
% $Id: spm_MDP.m 4804 2012-07-26 13:14:18Z karl $
 
% set up and preliminaries
%==========================================================================

% initial and final states
%==========================================================================
T     = 8;                         % number of offers
Ns    = 5;                         % number of outcomes (hidden states)
Nu    = 2;                         % number of actions (hidden controls)
a     = 8/(2*T);                   % probability of a high offer
b     = 1/(8*T);                       % probability of withdrawn offer

% iinitial state
%--------------------------------------------------------------------------
S     = sparse(1,1,1,Ns,1);
 
% priors over final state
%--------------------------------------------------------------------------
C     = [0 0 0 1 1]';

% (uniform) cost over control (d)
%--------------------------------------------------------------------------
D     = ones(Nu,1);

% transition probabilities (P{1} - decline; P{2} - accept)
%--------------------------------------------------------------------------
P{1}  = [(1 - a - b) 0 0 0 0;
          a          0 0 0 0;
          b          1 1 0 0;
          0          0 0 1 0;
          0          0 0 0 1];
      
P{2}  = [ 0 0 0 0 0;
          0 0 0 0 0;
          0 0 1 0 0;
          1 0 0 1 0;
          0 1 0 0 1];
      

% solve
%==========================================================================
MDP.T = T;                         % process depth (the horizon)
MDP.S = S;                         % initial state
MDP.P = P;                         % transition probabilities (priors)
MDP.C = C;                         % terminal cost probabilities (priors)
MDP.D = D;                         % control probabilities (priors)

[Q,R,S,E] = spm_MDP(MDP);

return
 
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
 


