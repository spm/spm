function spm_MDP_trust
% Demo for active inference with limited offer game
%__________________________________________________________________________
%
% This demonstration routine uses variational Bayes to minimise the free
% energy to model decision-making. The particular focus here is on
% decisions that are time-sensitive, requiring an explicit representation
% of future states. The example considered here represents a limited offer
% game, where a low offer can be converted to a high offer, which may or
% may not occur. Furthermore, offers may be withdrawn. The objective is
% to understand model choices about accepting or declining the current
% offer in terms of active inference, under prior beliefs about future
% states. The model is specified in a fairly general way in terms of
% probability transition matrices and beliefs about future states. The
% particular inversion scheme used here is spm_MDP_select, which uses a
% mean-field approximation between hidden control and hidden states. It is
% assumed that the agent believes that it will select a particular action
% (accept or decline) at a particular time.
%
% We run an exemplar game, examine the distribution of time to acceptance
% as a function of different beliefs (encoded by parameters of the
% underlying Markov process) and demonstrate how the model can be used to
% produce trial-specific changes in uncertainty – or how one can use
% behaviour to identify the parameters used by a subject.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_trust.m 5318 2013-03-08 17:38:06Z karl $

% set up and preliminaries
%==========================================================================
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
    
         
% initial state – encoding the actual type of trustee
%--------------------------------------------------------------------------
S    = [0 1]';                    % indicator - [prosocial nonsocial]
S    = kron(S,[1 0 0 0 0]');


% prior beliefs about initial state
%--------------------------------------------------------------------------
k    = [1 1];
p    = softmax(k(:));
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

% allowable policies – sequences of control (of depth T)
%--------------------------------------------------------------------------
V    = [1 2;
        1 1];

 
% MDP Structure
%==========================================================================
MDP.N = 8;                          % number of variational iterations
MDP.K = 0;                          % no memory
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
    Q(:,i) = softmax(k);
    O(:,i) = MDP.O(:,end);
    P(:,i) = MDP.P(:,end);
    W(:,i) = MDP.W(:,end);
    
    % update prior beliefs about initial state (context)
    %----------------------------------------------------------------------
    a      = find(MDP.U);
    p      = MDP.O(:,2)'*MDP.A*MDP.B{a}(:,[1 6]);
    p      = p(:)/sum(p);
    
    % update prior beliefs about initial state (context)
    %----------------------------------------------------------------------
    k      = k + log(p);
    p      = softmax(k);
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
legend({'prosocial','nnonsocial'})

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

