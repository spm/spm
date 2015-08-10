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
% (that take the agent to the four locations) and 16 outcomes (four
% locations times two cues times two rewards).  The central location has an
% ambiguous or uninformative cue outcome, while the upper arms are rewarded
% probabilistically with an 80% schedule.
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
% $Id: DEM_demo_MDP_habits.m 6517 2015-08-10 11:21:53Z karl $

% set up and preliminaries
%==========================================================================
% rng('default')

% observation probabilities
%--------------------------------------------------------------------------
a      = .9;
A{1,1} = [.5 .5; .5 .5];
A{2,2} = [a (1 - a); (1 - a) a];
A{3,3} = [(1 - a) a; a (1 - a)];
A{4,4} = [1 0; 0 1];
A      = spm_cat(A);

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
C  = c*[-1 -1 1 -1 -1 1 0 0]';

% prior beliefs about initial state
%--------------------------------------------------------------------------
d  = kron([1 0 0 0],[4 4])' + 1;

% allowable policies (of depth T)
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  3  4  4  4  4
      1  2  3  4  2  3  1  2  3  4];

% true initial states
%--------------------------------------------------------------------------
nt    = 16;                         % number of trials
s     = [1 1 0 1 1 1 1 1 1 1 1 0 1 1 1 1];
for i = 1:nt
    p    = rand > 1/2;
    p    = s(i);
    S{i} = kron([1 0 0 0],[p (1 - p)])';
end


% MDP Structure
%==========================================================================
for i = 1:length(S)
    
    MDP(i).V = V;                   % allowable policies
    MDP(i).S = S{i};                % true initial state
    MDP(i).A = A;                   % observation model
    MDP(i).B = B;                   % transition probabilities (priors)
    MDP(i).C = C;                   % terminal cost probabilities (priors)
    MDP(i).d = d;                   % initial state probabilities (priors)
    
    MDP(i).alpha  = 64;             % gamma hyperparameter
    MDP(i).beta   = 4;              % gamma hyperparameter

end

% Solve - an example game: a run of reds with an oddball
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
OPTIONS.plot  = gcf;
OPTIONS.habit = 0;
MDP           = spm_MDP_VB(MDP,OPTIONS);

spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP)

spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_LFP(MDP(1:2),[4 6;3 3])


return


