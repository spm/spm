function MDP = DEM_demo_MDP_habits
% Demo of active inference for trust games
%__________________________________________________________________________
%
% This routine uses the Markov decision process formulation of active
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
% (that take the agent to the four locations) and 8 outcomes (four
% locations times two cues times two rewards).  The central location has an
% ambiguous or uninformative cue outcome, while the upper arms are rewarded
% probabilistically with an 90% schedule.
%
% A single trial is simulated followed by an examination of dopaminergic
% responses to conditioned and unconditioned stimuli (cues and rewards). A
% hierarchical version is then implemented, in which the mapping between
% locations in the generative model and the generative process is unknown
% and has to be learned.
%
% see also: spm_MPD_game
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_MDP_habits.m 6523 2015-08-15 20:57:28Z karl $

% set up and preliminaries
%==========================================================================
rng('default')

% observation probabilities
%--------------------------------------------------------------------------
a      = .95;
b      = 1 - a;
A      = [1 1 0 0 0 0 0 0;
          0 0 a b 0 0 0 0;
          0 0 b a 0 0 0 0;
          0 0 0 0 b a 0 0;
          0 0 0 0 a b 0 0;
          0 0 0 0 0 0 1 0;
          0 0 0 0 0 0 0 1];

% transition probabilities: states = [centre,L,R,cue] x [R,L]
%--------------------------------------------------------------------------
for i = 1:4
    B{i} = zeros(4,4);
    B{i}([2 3],[2 3]) = eye(2);
    B{i}(i,[1 4])     = 1;
    B{i} = kron(B{i},eye(2));
end

% priors: (utility)
%--------------------------------------------------------------------------
c  = 2;
C  = c*[-2 -2 1 -1 -1 1 0 0]';

% prior beliefs about initial state
%--------------------------------------------------------------------------
d  = kron([1 0 0 0],[1 1])' + 1;

% allowable policies (of depth T)
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  3  4  4  4  4
      1  2  3  4  2  3  1  2  3  4];

% true initial states
%--------------------------------------------------------------------------
s     = [2 1 2 1 1 1 1 1 1 1 1 2 1 1 1 1];

% MDP Structure
%==========================================================================
for i = 1:length(s)
    
    MDP(i).V = V;                   % allowable policies
    MDP(i).s = s(i);                % true initial state
    MDP(i).A = A;                   % observation model
    MDP(i).B = B;                   % transition probabilities (priors)
    MDP(i).C = C;                   % terminal cost probabilities (priors)
    MDP(i).d = d;                   % initial state probabilities (priors)
    
    MDP(i).alpha  = 8;              % gamma hyperparameter
    MDP(i).beta   = 1;              % gamma hyperparameter

end

% Solve - an example game: a run of reds then an oddball
%==========================================================================
MDP  = spm_MDP_VB(MDP);

% illustrate behavioural respsonses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_game(MDP(1:16));

% illustrate phase-precession and responses to chosen option - 1st trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP(1),[4 6;3 3]);

% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_LFP(MDP(1:8));

% illustrate oddball responses (MMN) - CS
%--------------------------------------------------------------------------
MMN  = MDP(1:2);
MMN(1).s = [2 8 6];
MMN(2).s = [2 8 6];
MMN(1).d = kron([1 0 0 0],[1 1])' + 1;
MMN(2).d = kron([1 0 0 0],[8 1])' + 1;
MMN(1).h = [];
MMN(2).h = [];

MMN(1)  = spm_MDP_VB(MMN(1));
MMN(2)  = spm_MDP_VB(MMN(2));
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MMN,[1 2;1 1]);

% illustrate oddball responses (MMN) - US
%--------------------------------------------------------------------------
MMN(1).s = [2 8 6];
MMN(2).s = [2 8 5];
MMN(1).d = kron([1 0 0 0],[1 4])' + 1;
MMN(2).d = kron([1 0 0 0],[1 4])' + 1;
MMN(1).h = [];
MMN(2).h = [];

MMN(1)  = spm_MDP_VB(MMN(1));
MMN(2)  = spm_MDP_VB(MMN(2));
spm_figure('GetWin','Figure 5'); clf
spm_MDP_VB_LFP(MMN,[5 6;3 3]);


% illustrate habit leaning and repetition supression
%==========================================================================

% true initial states
%--------------------------------------------------------------------------
clear MDP
s     = [ones(1,24) 2*ones(1,128)];
d     = kron([1 0 0 0],[1 1])' + 1;

% MDP Structure
%==========================================================================
for i = 1:length(s)
    
    MDP(i).V = V;                   % allowable policies
    MDP(i).s = s(i);                % true initial state
    MDP(i).A = A;                   % observation model
    MDP(i).B = B;                   % transition probabilities (priors)
    MDP(i).C = C;                   % terminal cost probabilities (priors)
    MDP(i).d = d;                   % initial state probabilities (priors)
    
    MDP(i).alpha  = 8;              % gamma hyperparameter
    MDP(i).beta   = 1;              % gamma hyperparameter

end

OPTIONS.habit = 0;
MDP  = spm_MDP_VB(MDP,OPTIONS);


spm_figure('GetWin','Figure 6'); clf
spm_MDP_VB_game(MDP);




[B0,BV] = spm_MDP_DP(MDP(1));


return


