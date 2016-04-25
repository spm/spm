function MDP = DEM_MDP_decision
% Demo of active inference for visual salience
%__________________________________________________________________________
%
% This routine illustrates the treatment of signal detection paradigms in
% the context of active inference and Markov decision processes. This is
% probably one of the simplest paradigms to model; in which there are just
% too  hidden states generating ambiguous stimuli - and the agent move from
% an undecided (hidden) state to a definitive choice. The A tensor in this
% instanceen codes ambiguity (perceptual noise), while the B matrix encodes
% the behaviour-dependent transitions among decision states. Finally,
% the C matrix  encodes  prior costs or preferences. In this instance, the
% agent does not want to be wrong – and prefers to be right.
%
% in what follows, we simulate a single trial to illustrate the underlying
% Bayesian belief updating and associated behavioural and physiological
% responses. We then consider multiple trials under different levels of
% ambiguity and cost. The dependent measures in this instance include the
% behavioural accuracy, reaction times (assuming 250 ms time bins) and the
% uncertainty about the cause of sensory cues and control – as measured by
% the entropy of posterior beliefs prior to making a choice.
%
% see also: DEM_demo_MDP_rule.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_MDP_decision.m 6779 2016-04-25 09:40:35Z karl $

% set up and preliminaries
%==========================================================================
% rng('default')

% [D]efault prior beliefs about initial states (in terms of counts): D
%--------------------------------------------------------------------------
D{1} = [1 1]';        % rule:   {'left','right'}
D{2} = [0 0 1]';      % report: {'left','right','undecided'}

% [A]mbiguous mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
Ns    = zeros(1,Nf);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
r      = 2;                            % odds ratio (ambiguity)
for f1 = 1:Ns(1)                       % rule
    for f2 = 1:Ns(2)                   % report
        
        % A{1} what: {'left','right'}
        %------------------------------------------------------------------
        A{1}(1:2,f1,f2) = 1;
        A{1}(f1, f1,f2) = r;
        
        
        % A{2} feedback: {'null','right','wrong'}
        %------------------------------------------------------------------
        if f2 == 3,
            A{2}(1,f1,f2) = 1; % undecided
        elseif f1 == f2
            A{2}(2,f1,f2) = 1; % right
        else
            A{2}(3,f1,f2) = 1; % wrong
        end
        
    end
end
Ng    = numel(A);
for g = 1:Ng
    No(g) = size(A{g},1);
    A{g}  = double(A{g});
end

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% control states B(2): {'left','right','undecided'}
%--------------------------------------------------------------------------
B{2}(:,:,1) = [1 0 1;
               0 1 0;
               0 0 0];
B{2}(:,:,2) = [1 0 0;
               0 1 1;
               0 0 0];
B{2}(:,:,3) = [1 0 0;
               0 1 0;
               0 0 1];

% allowable policies (specified as the next action) U
%--------------------------------------------------------------------------
U(1,1,:)  = [1 1]';         % choose left
U(1,2,:)  = [1 2]';         % choose right
U(1,3,:)  = [1 3]';         % keep watching

% priors: (utility) C; the agent expects to avoid mistakes
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end
% and expects itself to make a decision after the fifth observation
%--------------------------------------------------------------------------
c    = 1;
C{2} = [ 0  0  0  0  0  0  0;
         c  c  c  c  c  c  c;
        -4 -4 -4 -4 -4 -4 -4];

% MDP Structure
%--------------------------------------------------------------------------
T     = size(C{2},2);
mdp.T = T;                      % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states


mdp.Aname = {'cue', 'feedback'};
mdp.Bname = {'rule','decision'};
mdp.alpha = 16;
mdp.tau   = 1;

mdp  = spm_MDP_check(mdp);


% illustrate a single trial
%==========================================================================
MDP  = spm_MDP_VB_X(mdp);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);

% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],1);

return

% illustrate a sequence of trials
%==========================================================================
clear MDP

% create structure array
%--------------------------------------------------------------------------
N     = 32;
for i = 1:N
    MDP(i) = mdp;
end


% Solve an example sequence undercover initial states
%==========================================================================
MDP   = spm_MDP_VB_X(MDP);




% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_game(MDP);


subplot(2,1,2)
subplot(4,1,3),
plot(1:N,Fs,1:N,F),  xlabel('trial'), spm_axis tight
title('States ','Fontsize',16), legend({'Surprise','Total Free energy'})
subplot(4,1,4), spm_MDP_plot_moves(MDP)
plot(1:N,Fu,1:N,Fq), xlabel('trial'), spm_axis tight
title('Action ','Fontsize',16), legend({'saccades','hits'})




