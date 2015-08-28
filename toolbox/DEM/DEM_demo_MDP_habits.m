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
% $Id: DEM_demo_MDP_habits.m 6538 2015-08-28 12:54:40Z karl $

% set up and preliminaries
%==========================================================================
rng('default')

% allowable policies (of depth T)
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  3  4  4  4  4
      1  2  3  4  2  3  1  2  3  4];
  
% outcome probabilities
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
C  = [-1 -1 c -c -c c 0 0]';

% prior beliefs about initial state
%--------------------------------------------------------------------------
d  = kron([1 0 0 0],[1 1])';

% prior beliefs about policies
%--------------------------------------------------------------------------
e  = [8 + zeros(size(V,2),1); 1];


% MDP Structure array
%==========================================================================
mdp.V = V;                    % allowable policies
mdp.A = A;                    % observation model
mdp.B = B;                    % transition probabilities (priors)
mdp.C = C;                    % terminal cost probabilities (priors)
mdp.d = d;                    % initial state probabilities (priors)
mdp.s = 1;                    % initial state

% true initial states – with context change at trial 12
%--------------------------------------------------------------------------
i           = [1,3];          % change context in a few trials
[MDP(1:24)] = deal(mdp);      % create structure array
[MDP(i).s]  = deal(2);        % change context
MDP(12).s   = [1 7 5];        % unexpected state transition
MDP(23).u   = 4;              % go to cue (CS)
MDP(24).u   = 2;              % go straight to reward (US)


% Solve - an example game: a run of reds then an oddball
%==========================================================================
OPTIONS.habit = 0;
MDP  = spm_MDP_VB(MDP,OPTIONS);

% illustrate behavioural respsonses – single trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1a'); clf
spm_MDP_VB_trial(MDP(1));

% illustrate behavioural respsonses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1b'); clf
spm_MDP_VB_game(MDP);

% illustrate phase-precession and responses to chosen option - 1st trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP(1),[4 6;3 3]);

% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_LFP(MDP(1:8));

% illustrate oddball responses (MMN) - US
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MDP([11,12]),[5 6;3 3]);
subplot(4,1,1), title('Mismatch negativity','FontSize',16)

% illustrate oddball responses (MMN) - CS
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
u = spm_MDP_VB_LFP(MDP([2,11]),[1 2;1 1]);
subplot(4,1,1), title('Repetition suppression','FontSize',16)

spm_figure('GetWin','Figure 5a');clf, t = (1:16)*64;
subplot(2,2,1),plot(t,u{1}{2,1},'b',t,u{2}{2,1},'b-.'), grid on, axis square
xlabel('Time (ms)'),ylabel('LFP'),title('Early and late responses', 'FontSize',16)
subplot(2,2,2),plot(t,u{2}{2,1} - u{1}{2,1}),           grid on, axis square ij
xlabel('Time (ms)'),ylabel('LFP'),title('Difference waveform (MMN)','FontSize',16)
subplot(2,2,1), legend({'oddball','standard'})

% illustrate transfer of dopamine responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 6'); clf
spm_MDP_VB_LFP(MDP([23,24]),[3 4;2 2]);
subplot(4,1,1), title('Transfer of dopamine responses','FontSize',16)



% illustrate reversal learning - after trial 32
%==========================================================================
clear MDP,   c = 4;

[MDP(1:64)]    = deal(mdp);
[MDP(32:64).s] = deal(2);
[MDP.C]        = deal([-1 -1 c -c -c c 0 0]');

OPTIONS.habit  = 0;
spm_figure('GetWin','Figure 7'); clf
spm_MDP_VB_game(spm_MDP_VB(MDP,OPTIONS));

OPTIONS.habit  = 1;
spm_figure('GetWin','Figure 8'); clf
spm_MDP_VB_game(spm_MDP_VB(MDP,OPTIONS));


% the effect of prior exposure
%--------------------------------------------------------------------------
clear MDP,   c = 3; 
[MDP(1:16)]    = deal(mdp);
[MDP(4:16).s]  = deal(2);
[MDP.C]        = deal([-1 -1 c -c -c c 0 0]');
OPTIONS.habit  = 0;
OPTIONS.plot   = 0;

d     = 1:4:32;
for i = 1:length(d)
    MDP(1).d = kron([1 0 0 0],[d(i) 1])' + 1;
    MDP      = spm_MDP_VB(MDP,OPTIONS);
    Q        = spm_MDP_VB_game(MDP);
    ext(i)   = find([Q.R(9,4:end) 1],1);
end

spm_figure('GetWin','Figure 9'); clf
subplot(2,1,1); bar(d,ext,'c'), axis square
xlabel('Previous exposures'), ylabel('Trials until reversal')
title('Reversal learning','FontSize',16)



% illustrate devaluation
%==========================================================================

% devalue (i.e. switch) preferences after after habitization (trial 32)
%--------------------------------------------------------------------------
clear MDP; OPTIONS.habit = 1; c = 2;

[MDP(1:48)]     = deal(mdp);
[MDP.C]         = deal([-1 -1 c -c -c  c 0 0]');
[MDP(32:end).C] = deal([-1 -1 -c c  c -c 0 0]');
MDP(1).e        = e;
MDP(1).d        = kron([1 0 0 0],[8 8])';

spm_figure('GetWin','Figure 10'); clf
spm_MDP_VB_game(spm_MDP_VB(MDP,OPTIONS));

% repeat but before habit formation (1216)
%--------------------------------------------------------------------------
[MDP(8:end).C]  = deal([-1 -1 -c c  c -c 0 0]');
MDP(1).e        = e;
MDP(1).d        = kron([1 0 0 0],[8 8])';

spm_figure('GetWin','Figure 11'); clf
spm_MDP_VB_game(spm_MDP_VB(MDP,OPTIONS));



% illustrate epistemic habit leaning
%==========================================================================

% true initial states
%--------------------------------------------------------------------------
clear MDP; OPTIONS.habit = 1;
i            = rand(1,32) > 1/2;
[MDP(1:32)]  = deal(mdp);
[MDP(i).s]   = deal(2);
MDP(1).e     = e;

% habitual (non-sequential) policy
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 12a'); clf
MDP = spm_MDP_VB(MDP,OPTIONS); spm_MDP_VB_game(MDP);
h   = MDP(end).c;
h   = h*diag(1./sum(h));

subplot(3,3,7); image(64*(1 - h)), axis square
xlabel('Hidden state'), xlabel('Hidden state')
title('Epistemic habit','FontSize',16)

% get equivalent dynamic programming solutions
%--------------------------------------------------------------------------
[B0,BV] = spm_MDP_DP(MDP(1));

subplot(3,3,8); image(64*(1 - spm_softmax(B0))), axis square
xlabel('Hidden state'), xlabel('Hidden state')
title('Active inference','FontSize',16)

subplot(3,3,9); image(64*(1 - BV)), axis square
xlabel('Hidden state'), xlabel('Hidden state')
title('Dynamic programming','FontSize',16)


% now repeat with unambiguous outcomes
%--------------------------------------------------------------------------
A      = [1 0 0 0 0 0 0 0;
          0 1 0 0 0 0 0 0
          0 0 a b 0 0 0 0;
          0 0 b a 0 0 0 0;
          0 0 0 0 b a 0 0;
          0 0 0 0 a b 0 0;
          0 0 0 0 0 0 1 0;
          0 0 0 0 0 0 0 1];
      
[MDP.A]  = deal(A);

% habitual (non-sequential) policy
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 12b'); clf
MDP = spm_MDP_VB(MDP(1:32),OPTIONS); spm_MDP_VB_game(MDP);
h   = MDP(end).c;
h   = h*diag(1./sum(h));

subplot(3,3,7); image(64*(1 - h)), axis square
xlabel('Hidden state'), xlabel('Hidden state')
title('Epistemic habit','FontSize',16)

% get equivalent dynamic programming solutions
%--------------------------------------------------------------------------
[B0,BV] = spm_MDP_DP(MDP(1));

subplot(3,3,8); image(64*(1 - spm_softmax(B0))), axis square
xlabel('Hidden state'), xlabel('Hidden state')
title('Active inference','FontSize',16)

subplot(3,3,9); image(64*(1 - BV)), axis square
xlabel('Hidden state'), xlabel('Hidden state')
title('Dynamic programming','FontSize',16)

return


