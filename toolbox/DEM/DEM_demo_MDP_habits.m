function MDP = DEM_demo_MDP_habits
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
% (that take the agent to the four locations) and five outcomes (two
% locations times two cues plus the centre).  The central location has an
% ambiguous or uninformative cue outcome, while the upper arms are rewarded
% probabilistically.
%
% This version  focuses on learning by optimising the parameters of the
% generative model. In particular, it looks at the acquisition of epistemic
% habits  – and how they relate to optimal policies under dynamic
% programming. We start with a series of simulations to illustrate various
% phenomena in electrophysiology and then move on to learning per se.
%
% see also: spm_MPD_game
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_habits.m 6707 2016-01-31 13:16:28Z karl $
 
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
A      = [1 1 0 0 0 0 0 0;    % ambiguous starting position (centre)
          0 0 a b 0 0 0 0;    % reward - right
          0 0 b a 0 0 0 0;    % no reward - right
          0 0 0 0 b a 0 0;    % reward - left
          0 0 0 0 a b 0 0;    % no reward - left
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
% probabilities. Here, the agent prefers rewards to losses.
%--------------------------------------------------------------------------
c  = 3;
C  = [0 0  0  0  0 0 0;
      0 c -c  c -c 0 0;
      0 c -c  c -c 0 0]';
 
% now specify prior beliefs about initial state, in terms of counts
%--------------------------------------------------------------------------
d  = kron([1 0 0 0],[8 8])';
 
% allowable policies (of depth T).  These are just sequences of actions
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  3  4  4  4  4
      1  2  3  4  2  3  1  2  3  4];
 
% prior beliefs about policies
%--------------------------------------------------------------------------
Np = size(V,2);
e  = [4 + zeros(Np,1); 1];
 
 
% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.V = V;                    % allowable policies
mdp.A = A;                    % observation model
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.d = d;                    % prior over initial states
mdp.s = 1;                    % true initial state
 
% true initial states – with context change at trial 12
%--------------------------------------------------------------------------
i           = [1,3];          % change context in a couple of trials
[MDP(1:32)] = deal(mdp);      % create structure array
[MDP(i).s]  = deal(2);        % deal context changes
% [MDP.C]   = deal(C - C);    % for epistemic simulation
MDP(12).o   = [1 6 7];        % unexpected outcome


 
% Solve - an example game: a run of reds then an oddball
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

% illustrate phase-precession and responses to chosen option - 1st trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP(1),[4 6;3 3]);

% place cells
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
subplot(2,2,1),spm_MDP_VB_place_cell(MDP(1:6),[3 6;3 3]);
subplot(2,2,2),spm_MDP_VB_place_cell(MDP(1:6),[7 8;2 2]);

% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MDP(1:8));

% illustrate oddball responses (P300) - US
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
spm_MDP_VB_LFP(MDP([11,12]),[8;3]);
subplot(4,1,1), title('Violation response (P300)','FontSize',16)
 
% illustrate oddball responses (MMN)  - CS and dopamine transfer
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 6a'); clf
i  = find(ismember(spm_cat({MDP.u}'),[4 2],'rows'));
spm_MDP_VB_LFP(MDP([i(1),i(end)]),[1;1]);
subplot(4,1,1), title('Repetition suppression and DA transfer','FontSize',16)
 
spm_figure('GetWin','Figure 6b');clf
v  = spm_MDP_VB_LFP(MDP([i(1),i(end)]),[1;1]);
t  = (1:16)*16 + 80;
subplot(2,1,1),plot(t,v{1}{2,1},'b-.',t,v{2}{2,1},'b:',t,v{2}{2,1} - v{1}{2,1})
xlabel('Time (ms)'),ylabel('LFP'),title('Difference waveform (MMN)','FontSize',16)
legend({'oddball','standard','MMN'}), grid on, axis square

% return

% illustrate reversal learning - after trial 32
%==========================================================================
clear MDP
[MDP(1:64)]    = deal(mdp);
[MDP(32:64).s] = deal(2);
 
% just learn the context
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7'); clf
spm_MDP_VB_game(spm_MDP_VB(MDP));

 
% the effect of prior exposure on reversal learning
%--------------------------------------------------------------------------
clear MDP
[MDP(1:16)]    = deal(mdp);
[MDP(4:16).s]  = deal(2);
OPTIONS.plot   = 0;
 
d     = linspace(32,64,4);
for i = 1:length(d)
    MDP(1).d   = kron([1 0 0 0],[d(i) 8])';
    M          = spm_MDP_VB(MDP,OPTIONS);
    Q          = spm_MDP_VB_game(M);
    ext(i)     = sum(Q.O(4:end) == 3);
end
 
spm_figure('GetWin','Figure 8'); clf
subplot(2,1,1); bar(d,ext,'c'), axis square
xlabel('Previous exposures'), ylabel('Trials until reversal')
title('Reversal learning','FontSize',16)

 
% illustrate devaluation: enable habit learning from now on
%==========================================================================
OPTIONS.habit   = 1;
mdp.e           = e;
N               = 64;

% devalue (i.e. switch) preferences after habitisation (trial 32)
%--------------------------------------------------------------------------
clear MDP; 
 
[MDP(1:N)]      = deal( mdp);
[MDP(48:end).C] = deal(-mdp.C);
 
spm_figure('GetWin','Figure 9'); clf
spm_MDP_VB_game(spm_MDP_VB(MDP,OPTIONS));
 
% repeat but now devalue before habit formation (at trial 8)
%--------------------------------------------------------------------------
[MDP(16:end).C]  = deal(-mdp.C);
 
spm_figure('GetWin','Figure 10'); clf
spm_MDP_VB_game(spm_MDP_VB(MDP,OPTIONS));
 
 
 
% illustrate epistemic habit leaning
%==========================================================================
 
% true initial states
%--------------------------------------------------------------------------
clear MDP;
i            = rand(1,N) > 1/2;
[MDP(1:N)]   = deal(mdp);
[MDP(i).s]   = deal(2);

 
% habitual (non-sequential) policy
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 11'); clf
M   = spm_MDP_VB(MDP,OPTIONS); spm_MDP_VB_game(M);
h   = M(end).c;
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
A = [1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0
     0 0 a b 0 0 0 0;
     0 0 b a 0 0 0 0;
     0 0 0 0 b a 0 0;
     0 0 0 0 a b 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1];

C = [0 0 0  0 0  0 0 0;
     0 0 c -c c -c 0 0;
     0 0 c -c c -c 0 0]';
      
[MDP.A]  = deal(A);
[MDP.C]  = deal(C);
 
% habitual (non-sequential) policy
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 12'); clf
M = spm_MDP_VB(MDP,OPTIONS); spm_MDP_VB_game(M);
h   = M(end).c;
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

function spm_MDP_VB_place_cell(MDP,UNITS)

% place cell plotting subroutine
%--------------------------------------------------------------------------
col = {'r','g','b','m','c'};
for t = 1:length(MDP)
    for j = 1:size(UNITS,2)
        
        [~,v] = spm_MDP_VB_LFP(MDP(t),UNITS(:,j));
        qu    = spm_cat(v);
        qe    = randn(2,size(qu,1))/4;
        L     = [0 0 -1 -1 1 1  0  0;
                 0 0  1  1 1 1 -1 -1];
        for i = 1:size(qu,1)
            X(:,i) = L(:,MDP(t).s(ceil(i/16))) + qe(:,i);
        end
        X     = spm_conv(X,0,3);
        plot(X(1,:),X(2,:),'r:'), hold on
        for i = 1:size(qu,1)
            if qu(i) > .80
                plot(X(1,i),X(2,i),'.','MarkerSize',16,'Color',col{j})
            end
        end
        
    end
end
axis([-2 2 -2 2]), axis square, hold off
title('Place field responses','fontsize',16)
