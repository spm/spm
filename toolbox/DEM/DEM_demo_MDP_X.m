function MDP = DEM_demo_MDP_X
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
% $Id: DEM_demo_MDP_X.m 6652 2015-12-21 10:51:54Z karl $
 
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
a      = .98;
b      = 1 - a;
A{1}(:,:,1) = [...
    1 1 1 0;    % cue neutral
    0 0 0 1;    % cue left
    0 0 0 0];   % cue right
A{1}(:,:,2) = [...
    1 1 1 0;    % cue neutral
    0 0 0 0;    % cue left
    0 0 0 1];   % cue right

A{2}(:,:,1) = [...
    1 0 0 1;    % reward neutral
    0 a b 0;    % reward positive
    0 b a 0];   % reward negative
A{2}(:,:,2) = [...
    1 0 0 1;    % reward neutral
    0 b a 0;    % reward positive
    0 a b 0];   % reward negative


% controlled transitions: B{u}
%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions of hidden states
% under each action or control state. Here, there are four actions taking the
% agent directly to each of the four locations.
%--------------------------------------------------------------------------
B{1}(:,:,1)  = [1 0 0 1; 0 1 0 0;0 0 1 0;0 0 0 0];
B{1}(:,:,2)  = [0 0 0 0; 1 1 0 1;0 0 1 0;0 0 0 0];
B{1}(:,:,3)  = [0 0 0 0; 0 1 0 0;1 0 1 1;0 0 0 0];
B{1}(:,:,4)  = [0 0 0 0; 0 1 0 0;0 0 1 0;1 0 0 1];

% contextthe point that
%--------------------------------------------------------------------------
B{2}         = eye(2);

% priors: (utility) C
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of  log
% probabilities over outcomes. Here, the agent prefers rewards to losses.
%--------------------------------------------------------------------------
c     = 3;
C{1}  = [0  0  0;
         0  0  0;
         0  0  0];

C{2}  = [0  0  0;
         0  c  c;
         0 -c -c];
 
% now specify prior beliefs about initial state, in terms of counts
%--------------------------------------------------------------------------
d{1} = [1 0 0 0]';
d{2} = [1 1]';


% allowable policies (of depth T).  These are just sequences of actions
%--------------------------------------------------------------------------
V(:,:,1) = [1  1  1  1  2  3  4  4  4  4
            1  2  3  4  2  3  1  2  3  4];
V(:,:,2) = 1;
 

 
% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.V = V;                    % allowable policies
mdp.A = A;                    % observation model
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred outcomes
mdp.d = d;                    % prior over initial states
mdp.s = [1 1]';               % true initial state

mdp.Aname = {'exteroceptive','interoceptive'};
mdp.Bname = {'position','context'};
 
% true initial states – with context change at trial 12
%--------------------------------------------------------------------------
i              = [1,3];          % change context in a couple of trials
[MDP(1:32)]    = deal(mdp);      % create structure array
[MDP(i).s]     = deal([1 2]');   % deal context changes


% Solve - an example game: a run of reds then an oddball
%==========================================================================
MDP  = spm_MDP_VB_X(MDP);
 
% illustrate behavioural responses – single trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1));

% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_game(MDP);

% illustrate phase-precession and responses to chosen option - 1st trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_LFP(MDP(1),[2 3;3 3],1);

% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MDP(1:8));
 
% illustrate oddball responses (MMN) - CS
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
spm_MDP_VB_LFP(MDP([2,16]),[1;2],2);
subplot(4,1,1), title('Repetition suppression and DA transfer','FontSize',16)
 
spm_figure('GetWin','Figure 6'); clf
v  = spm_MDP_VB_LFP(MDP([2,16]),[1;2],2);
t  = (1:16)*16 + 80;
subplot(2,1,1),plot(t,v{1}{2,1},'b-.',t,v{2}{2,1},'b:',t,v{2}{2,1} - v{1}{2,1})
xlabel('Time (ms)'),ylabel('LFP'),title('Difference waveform (MMN)','FontSize',16)
legend({'oddball','standard','MMN'}), grid on, axis square

