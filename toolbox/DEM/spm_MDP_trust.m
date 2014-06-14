function spm_MDP_trust
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
% $Id: spm_MDP_trust.m 6044 2014-06-14 10:22:46Z karl $

% set up and preliminaries
%==========================================================================
rng('default')
%
% Payoffs (reward):
% _________________________________________________________________________
%                                 Trustee
% _________________________________________________________________________
% Investor              Cooperate (fT high) Defect (fT low)
% _________________________________________________________________________
% Cooperate            A (e.g.=26)         C (e.g.= 10)
% (fI high)            a (e.g.=26)         b (e.g.= 42)
%
% Defect               B (e.g.=21)         D (e.g.= 18)
% (fI low)             c (e.g.= 7)         d (e.g.= 10)
% _________________________________________________________________________

U{1} = [26 10;
        21 18]/8;                   % self payoff (utility)
     
U{2} = [26 42;
         7 10]/8;                   % other payoff (utility)
    
         
% initial state - encoding the actual type of trustee
%--------------------------------------------------------------------------
S    = [0 1]';                      % indicator - [prosocial nonsocial]
S    = kron(S,[1 0 0 0 0]');


% prior beliefs about initial state
%--------------------------------------------------------------------------
k    = [1 1];
p    = spm_softmax(k(:));
k    = log(p);
D    = kron(p,[1 0 0 0 0]');

         
% investor's payoffs or prior beliefs (softmax(utility))
%--------------------------------------------------------------------------
a    = 1/2;
pp   = [0; spm_softmax(spm_vec((1 - a)*U{1} + a*U{2}))];
pn   = [0; spm_softmax(spm_vec(        U{1}         ))];

C    = [pp*p(1); pn*p(2)];

% investor's belief (based on a prosocial and nonsocial trustee)
%--------------------------------------------------------------------------
cp   = spm_softmax(((1 - a)*U{2}(1,:) + a*U{1}(1,:))');
dp   = spm_softmax(((1 - a)*U{2}(2,:) + a*U{1}(2,:))');
cn   = spm_softmax((        U{2}(1,:)              )');
dn   = spm_softmax((        U{2}(2,:)              )');
 
% transition probabilities (B{1} - (c)ooperate; B{2} - (d)efect)
%--------------------------------------------------------------------------
B{1} = ...                        % cooperate:
   [0     0 0 0 0  0 0 0 0 0;     % start - prosocial trustee
    cp(1) 1 0 0 0  0 0 0 0 0;     % cc - prosocial - state 2
    0     0 1 0 0  0 0 0 0 0;     % dc - prosocial - state 3
    cp(2) 0 0 1 0  0 0 0 0 0;     % cd - prosocial - state 4
    0     0 0 0 1  0 0 0 0 0;     % dd - prosocial - state 5the
    
    0 0 0 0 0  0     0 0 0 0;     % cc - nonsocial
    0 0 0 0 0  cn(1) 1 0 0 0;     % cc - nonsocial
    0 0 0 0 0  0     0 1 0 0;     % dc - nonsocial
    0 0 0 0 0  cn(2) 0 0 1 0;     % cd - nonsocial
    0 0 0 0 0  0     0 0 0 1];    % dd - nonsocial
    
B{2} = ...                        % defect:
   [0     0 0 0 0  0 0 0 0 0;     % start - nonsocial trustee
    0     1 0 0 0  0 0 0 0 0;     % ...
    dp(1) 0 1 0 0  0 0 0 0 0;
    0     0 0 1 0  0 0 0 0 0;
    dp(2) 0 0 0 1  0 0 0 0 0;
    
    0 0 0 0 0  0     0 0 0 0;
    0 0 0 0 0  0     1 0 0 0;
    0 0 0 0 0  dn(1) 0 1 0 0;
    0 0 0 0 0  0     0 0 1 0;
    0 0 0 0 0  dn(2) 0 0 0 1];


% observation probabilities
%--------------------------------------------------------------------------
A    = kron([1 1],speye(5,5));

% allowable policies - (of depth T); here, simply defect will cooperate
%--------------------------------------------------------------------------
V    = [1 2;
        1 1];

 
% MDP Structure
%==========================================================================
MDP.N = 8;                          % number of variational iterations
MDP.T = 2;                          % process depth (one-shot game)
MDP.S = S;                          % true initial state
MDP.A = A;                          % observation model
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
MDP.D = D;                          % initial state probabilities (priors)
MDP.V = V;                          % allowable policies

MDP.alpha = 2;                      % gamma hyperparameters
MDP.beta  = 1/2;

% Solve - an example game
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
MDP.plot = gcf;
MDP      = spm_MDP_game(MDP);


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

