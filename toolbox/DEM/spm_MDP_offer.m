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
T     = 16;                         % number of offers
Ns    = 5;                          % number of outcomes (hidden states)
Nu    = 2;                          % number of actions (hidden controls)
Pa    = 1/2;                        % probability of a high offer
Pb    = 1/16;                       % probability of withdrawn offer

% transition probabilities (B{1} - decline; B{2} - accept)
%--------------------------------------------------------------------------
for i = 1:T
    
    a       = 1 - (1 - Pa)^(1/T);
    b       = i*Pb/T;
    B{i,1}  = [(1 - a - b) 0 0 0 0;
                a          0 0 0 0;
                b          1 1 0 0;
                0          0 0 1 0;
                0          0 0 0 1];
    
    B{i,2}  = [ 0 0 0 0 0;
                0 0 0 0 0;
                0 0 1 0 0;
                1 0 0 1 0;
                0 1 0 0 1];
end
      

% iinitial state
%--------------------------------------------------------------------------
S     = [1 0 0 0 0]';
 
% priors over final state
%--------------------------------------------------------------------------
C     = [0 0 0 1 4]';

% (uniform) cost over control (d)
%--------------------------------------------------------------------------
D     = [1 1]';


% solve
%==========================================================================
MDP.T = T;                         % process depth (the horizon)
MDP.S = S;                         % initial state
MDP.B = B;                         % transition probabilities (priors)
MDP.C = C;                         % terminal cost probabilities (priors)
MDP.D = D;                         % control probabilities (priors)
MDP.W = 4;                         % log-precision

[Q,R,S,E] = spm_MDP(MDP);

return
 



