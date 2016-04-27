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
% $Id: DEM_MDP_decision.m 6786 2016-04-27 19:38:30Z karl $

% set up and preliminaries
%==========================================================================

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
a      = exp(1/4);              % odds ratio (ambiguity)
for f1 = 1:Ns(1)                % rule
    for f2 = 1:Ns(2)            % report
        
        % A{1} what: {'left','right'}
        %------------------------------------------------------------------
        A{1}(1:2,f1,f2) = 1;
        A{1}(f1, f1,f2) = a;
        
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
T    = 6;
C{2}(1,1:T) = 0;                % null
C{2}(2,1:T) = 3;                % correct
C{2}(3,1:T) = -6;               % wrong


% MDP Structure
%--------------------------------------------------------------------------
mdp.T = T;                      % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.s = [1 3]'*ones(1,T);
mdp.o = [1 1]'*ones(1,T);
mdp.u = [1 3]'*ones(1,T);

mdp.Aname = {'cue', 'feedback'};
mdp.Bname = {'rule','decision'};
mdp.alpha = 2;

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


% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');

P  = squeeze(MDP.P);
P  = (1 - cumprod([1 P(3,1:end - 1)])).*P(1,:);
subplot(3,2,3), bar(P,'c')
xlabel('epoch'),ylabel('P(correct choice)')
title('Expected behaviour','FontSize',16)



% illustrate choice behaviour over different levels ambiguity and reward
%==========================================================================
a = (1:32)/64;                      % range of perceptual ambiguity
c = [1 2];                          % levels of reward
t = (2:T)/4;                        % time (250 ms bins)
for i = 1:numel(a)
    for j = 1:numel(c)
        
        % perceptual precision
        %------------------------------------------------------------------
        for f1 = 1:Ns(1)                       % rule
            for f2 = 1:Ns(2)                   % report
                MDP.A{1}(1:2,f1,f2) = 1;
                MDP.A{1}(f1, f1,f2) = exp(a(i));
            end
        end
        
        % choice precision
        %------------------------------------------------------------------
        MDP.C{2}(1,:) = C{2}(1,:)*c(j);       % wait
        MDP.C{2}(2,:) = C{2}(2,:)*c(j);       % right
        MDP.C{2}(3,:) = C{2}(3,:)*c(j);       % wrong
        
        % solve active inference scheme
        %------------------------------------------------------------------
        MDP     = spm_MDP_VB_X(MDP);
        
        % accuracy and reaction time
        %------------------------------------------------------------------
        P       = squeeze(MDP.P);
        P       = (1 - cumprod([1 P(3,1:end - 1)])).*P(1,:);
        RT(i,j) = t*P(:)/sum(P);
        AC(i,j) = P(end);
        
        % uncertainty
        %------------------------------------------------------------------
        P       = squeeze(MDP.X{1});
        HS(i,j) = mean(sum(-P.*log(P + eps)));
        P       = squeeze(MDP.P);
        HU(i,j) = mean(sum(-P.*log(P + eps)));
        
    end
end

% plot behavioural responses and uncertainty
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf

subplot(2,2,1),plot(a,RT),xlabel('log odds ratio'),ylabel('seconds')
title('reaction time',     'FontSize',16), axis square, legend({'low reward','high reward'})
subplot(2,2,2),plot(a,AC),xlabel('log odds ratio'),ylabel('proportion')
title('choice accuracy',   'FontSize',16), axis square
subplot(2,2,3),plot(a,HS),xlabel('log odds ratio'),ylabel('nats')
title('state uncertainty', 'FontSize',16), axis square
subplot(2,2,4),plot(a,HU),xlabel('log odds ratio'),ylabel('nats')
title('choice uncertainty','FontSize',16), axis square



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

% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
MDP = spm_MDP_VB_X(MDP);
spm_MDP_VB_game(MDP);






