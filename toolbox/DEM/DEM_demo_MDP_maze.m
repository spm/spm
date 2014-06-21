function MDP = DEM_demo_MDP_maze
% Demo of active inference for trust games
%__________________________________________________________________________
%
% This routine uses the Markov decision process formulation of active
% inference (with variational Bayes) to model a simple trust game. In trust
% games, one plays an opponent who can either cooperate or defect. The
% payoff contingencies depend upon the joint choices of you and your
% opponent, which in turn depend upon your inferences about the nature of
% the opponent (pro-social or non-social). This example illustrates single
% round games with a special focus on Bayesian belief updating between
% games. This is illustrated in terms of evidence accumulation about
% the nature of the opponent by using the posterior marginal distributions
% following one game as the prior distribution over beliefs about the
% opponent in the next. This accumulation is shown in the final figures.
%
% In this example, there are nine states. The first is a starting state
% and the subsequent eight states model the four combinations of
% cooperation and defection (between you and your opponent) under the
% prior beliefs that the opponent is either pro-social or non-social. 
% Initially, these prior beliefs are uninformative but are subsequently 
% informed through experience. prior beliefs about behaviour are based on
% relative entropy or KL divergence in the usual way - which requires the
% specification of utility functions over states based upon standard payoff
% tables in these sorts of games. It is interesting to see how precision
% or confidence in beliefs about choices, fluctuates with beliefs about
% the nature of one's opponent.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_maze.m 6063 2014-06-21 11:05:30Z karl $

% set up and preliminaries
%==========================================================================
% rng('default')

% observation probabilities
%--------------------------------------------------------------------------
A{1,1} = [.5 .5; .5 .5; 0 0; 0 0];
A{2,2} = [0 0; 0 0; .8 .2; .2 .8];
A{3,3} = [0 0; 0 0; .2 .8; .8 .2];
A{4,4} = [1 0; 0 1; 0 0; 0 0];
A      = spm_cat(A);

% transition probabilities
%--------------------------------------------------------------------------
for i = 1:4
    B{i} = zeros(4,4);
    B{i}([2 3],[2 3]) = eye(2);
    B{i}(i,[1 4])     = 1;
    B{i} = kron(B{i},eye(2));
end

% priors: softmax(utility)
%--------------------------------------------------------------------------
C    = spm_softmax(kron(ones(4,1),[0; 0; 0;0]));

% prior beliefs about initial state
%--------------------------------------------------------------------------
D    = kron([1 0 0 0],[1 1]/2)';

% true initial state
%--------------------------------------------------------------------------
S    = kron([1 0 0 0],[1 0])';


% allowable policies (of depth T)
%--------------------------------------------------------------------------
V = [];
for i = 1:4
    v(1,1) = i;
    for j = 1:4
        v(2,1) = j;
        v(3,1) = 1;
        V  = [V v];
    end
end

 
% MDP Structure
%==========================================================================
MDP.N = 8;                          % number of variational iterations
MDP.S = S;                          % true initial state
MDP.A = A;                          % observation model
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
MDP.D = D;                          % initial state probabilities (priors)
MDP.V = V;                          % allowable policies

MDP.alpha = 8;                      % gamma hyperparameters
MDP.beta  = 1/2;

% Solve - an example game
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
MDP.plot = gcf;
MDP      = spm_MDP_game(MDP);

return


% now iterate repeated games accumulating posterior beliefs
%==========================================================================
MDP.plot = 0;
NG       = 64;

for i = 1:NG
    
    % solve and marginalise over posterior beliefs about hidden states
    %----------------------------------------------------------------------
    MDP    = spm_MDP_game(MDP);
    Q(:,i) = spm_softmax(k);
    O(:,i) = MDP.O(:,end);
    W(:,i) = MDP.W(:,end);
    P(:,i) = MDP.P(:,1);
    
    % update prior beliefs about initial state (pro-social or
    % non-social) and associated utility functions
    %----------------------------------------------------------------------
    a      = find(MDP.U(:,1));
    p      = MDP.O(:,2)'*MDP.A*MDP.B{a}(:,[1 6]);
    p      = p(:)/sum(p);
    k      = k + log(p);
    p      = spm_softmax(k);
    MDP.D  = kron(p,[1 0 0 0 0]');
    MDP.C  = [pp*p(1); pn*p(2)];
    
end


% graphics
%==========================================================================
spm_figure('GetWin','Figure 2'); clf

% posterior beliefs about hidden states (prosocial versus nonsocial)
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(1:NG,Q)
title('Beliefs about other','FontSize',16)
xlabel('Number of games','FontSize',12)
ylabel('True and posterior expectations','FontSize',12)
spm_axis tight, axis square
legend({'prosocial','nonsocial'})

subplot(2,2,2)
plot(1:NG,P)
title('Beliefs about control','FontSize',16)
xlabel('Number of games','FontSize',12)
ylabel('True and posterior expectations','FontSize',12)
spm_axis tight, axis square
legend({'cooperate','defect'})

subplot(2,2,3)
imagesc(O)
title('Outcomes','FontSize',16)
xlabel('Number of games','FontSize',12)
ylabel('Observed outcome','FontSize',12)
axis square

subplot(2,2,4)
plot(W)
title('Precision','FontSize',16)
xlabel('Number of games','FontSize',12)
ylabel('Expected precision','FontSize',12)
axis square

