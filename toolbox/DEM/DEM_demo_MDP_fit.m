function MDP = DEM_demo_MDP_habits_fit
% Demo of active inference for trust games
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
% inference (with variational Bayes) to model foraging for information in a
% three arm maze.  This demo illustrates variational free energy
% minimisation in the context of Markov decision processes, where the agent
% is equipped with prior beliefs that it will minimise expected free energy
% in the future. This free energy is the free energy of future sensory
% states expected under the posterior predictive distribution. It can be
% regarded as a generalisation of the variational formulation of KL control
% in which information gain or epistemic value is formulated explicitly.
%
% In this example, the agent starts at the centre of a three way maze
% which is baited with a reward in one of the two upper arms. However, the
% rewarded arm changes from trial to trial.  Crucially, the agent can
% identify where the reward (US) is located by accessing a cue (CS) in the
% lower arm. This tells the agent whether the reward is on the left or the
% right upper arm.  This means the optimal policy would first involve
% maximising information gain or epistemic value by moving to the lower arm
% and then claiming the reward this signified. Here, there are eight hidden
% states (four locations times right or left reward), four control states
% (that take the agent to the four locations) and seven outcomes (three
% locations times two cues plus the centre).  The central location has an
% ambiguous or uninformative cue outcome, while the upper arms are rewarded
% probabilistically.
%
% this version  focuses on learning by optimising the parameters of the
% generative model. In particular, it looks at the acquisition of epistemic
% habits  – and how they relate to optimal policies under dynamic
% programming. We start with a series of simulations to illustrate various
% phenomena in electrophysiology and then move on to learning per se.
%
% see also: spm_MPD_game
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_fit.m 6586 2015-10-31 12:02:44Z karl $
 
% set up and preliminaries
%==========================================================================
rng('default')
  
% outcome probabilities: A
%--------------------------------------------------------------------------
% We start by specifying the probabilistic mapping from hidden states
% to outcomes. This necessarily requires one to think carefully about the
% hidden states [centre, right, left, cue] x [reward, loss] and the ensuing
% outcomes
%--------------------------------------------------------------------------
a      = .95;
b      = 1 - a;
A      = [1 1 0 0 0 0 0 0;    % ambiguous starting position (centre)
          0 0 a b 0 0 0 0;    % left arm selected and rewarded
          0 0 b a 0 0 0 0;    % left arm selected and not rewarded
          0 0 0 0 b a 0 0;    % right arm selected and not rewarded
          0 0 0 0 a b 0 0;    % right arm selected and rewarded
          0 0 0 0 0 0 1 0;    % informative cue - reward on right
          0 0 0 0 0 0 0 1];   % informative cue - reward on left
 
% controlled transitions: B{u}
%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions of hidden states
% under each action or control state. Here, there are four actions taking the
% agent directly to each of the four locations. Some of these locations are
% absorbing states; in that once entered, they cannot be left
%--------------------------------------------------------------------------
for i = 1:4
    B{i} = zeros(4,4);
    B{i}([2 3],[2 3]) = eye(2);
    B{i}(i,[1 4])     = 1;
    B{i} = kron(B{i},eye(2));
end
 
% priors: (utility) C
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of  log
% probabilities. Here, the agent does not like to be exposed in the centre
% and, clearly, prefers rewards to losses.
%--------------------------------------------------------------------------
c  = 2;
C  = [0 c -c c -c 0 0]';
 
% now specify prior beliefs about initial state, in terms of counts
%--------------------------------------------------------------------------
d  = kron([1 0 0 0],[1 1])';
 
% allowable policies (of depth T).  These are just sequences of actions
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  3  4  4  4  4
      1  2  3  4  2  3  1  2  3  4];
 
 
 
% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.V = V;                    % allowable policies
mdp.A = A;                    % observation model
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.d = d;                    % prior over initial states
mdp.s = 1;                    % initial state
 
% true initial states
%--------------------------------------------------------------------------
n           = 8;
i           = rand(1,n) > 1/2;
[MDP(1:8)]  = deal(mdp);
[MDP(i).s]  = deal(2);
MDP(1).d    = kron([1 0 0 0],[128 1])';
 
 
% Solve - generate data
%==========================================================================
MDP  = spm_MDP_VB(MDP);
 
% illustrate behavioural responses – single trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1a'); clf
spm_MDP_VB_trial(MDP(1));
 
% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1b'); clf
spm_MDP_VB_game(MDP);

% Invert to recover parameters
%==========================================================================
DCM.MDP   = mdp;                  % MDP model
DCM.field = {'alpha','beta','d'}; % parameter (field) names to optimise
DCM.U     = {MDP.o};              % observable outcomes (stimului)
DCM.Y     = {MDP.u};              % observabale responses (action)

DCM   = spm_dcm_mdp(DCM);



 

