function [Q,R,S,E] = spm_MDP(MDP)
% solves the active inference problem for Markov decision processes
% FROMAT [Q,R,S,E] = spm_MDP(MDP)
%
% MDP.T         - process depth (the horizon)
% MDP.S(N,1)    - initial state
% MDP.P{M}(N,N) - transition probabilities among hidden states (priors)
% MDP.C(N,1)    - terminal cost probabilities (prior N over hidden states)
% MDP.D(M,1)    - control probabilities (prior over M control states)
%
% MDP.G         - transition probabilities used to generate outcomes 
%                 (default: the prior transition probabilities)
% MDP.A(N,N)    - Likelihood of outcomes given hidden states
%                 (default: an identity mapping from states to outcomes)
%
% produces:
%
% Q(N,K,T) - an array of conditional (posterior) expectations over N hidden
%            states and time 1,...,T at time K
% R(M,K,T) - an array of conditional expectations over M control
%            states and time 1,...,T at time K.
% S(N,T)   - a sparse matrix of ones, encoding the state at time T
% E(M,T)   - a sparse matrix of ones, encoding the action at time T
%
% This routine provides solutions of active inference (minimisation of
% variational free energy)using a generative model based upon a Markov 
% decision process. This model and inference scheme is formulated
% in discrete space and time. This means that the generative model (and
% process) are  finite state  machines or hidden Markov models whose
% dynamics are given by transition probabilities among states. For
% simplicity, we assume an isomorphism between hidden states and outcomes,
% where the likelihood corresponds to a particular outcome conditioned upon
% hidden states. Similarly, for simplicity, this routine assumes that action
% and hidden controls are isomorphic. If the dynamics of transition
% probabilities of the true process are not provided, this routine will use
% the equivalent probabilities from the generative model.
%
% The transition probabilities are a cell array of probability transition
% matrices corresponding to each (discrete) the level of the control state.
%
% Mote that the conditional expectations are functions of time but also
% contain expectations about fictive states over time at each time point.
% To create time dependent transition probabilities, one can specify a
% function in place of the transition probabilities under different levels
% of control.
%
% partially observed Markov decision processes can be modelled by
% specifying a likelihood (as part of a generative model) and absorbing any
% probabilistic mapping between (isomorphic) hidden states and outcomes
% into the transition probabilities G.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP.m 5043 2012-11-07 21:57:56Z karl $
 
% set up and preliminaries
%==========================================================================
 
% plotting defaults
%--------------------------------------------------------------------------
try MDP.nograph, PLOT = 0; catch, PLOT = 1; end
if PLOT, spm_figure('GetWin','MDP'); clf,   end

 
% generative model and initial states
%--------------------------------------------------------------------------
T     = MDP.T;            % process depth (the horizon)
P     = MDP.P;            % transition probabilities (priors)
S     = spm_vec(MDP.S);   % initial state
C     = spm_vec(MDP.C);   % terminal cost probabilities (priors)
D     = spm_vec(MDP.D);   % control probabilities (priors)
Ns    = length(S);        % number of hidden states
Nu    = length(P);        % number of hidden controls
 
% generative process (assume the true process is the same as the model)
%--------------------------------------------------------------------------
try
    G = MDP.G;
catch
    G = MDP.P;
end
 
% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A = MDP.A;
catch
    A = speye(Ns,Ns);
end
 
 
% get log-transforms and ensure normalization
%--------------------------------------------------------------------------
A     = A*diag(1./sum(A));
lnA   = log(A + eps);
for i = 1:Nu
    B{i}   = P{i}';
    P{i}   = P{i}*diag(1./max(sum(P{i}),eps)); % Push-forward
    B{i}   = B{i}*diag(1./max(sum(B{i}),eps)); % Pull-back
    lnB{i} = log(B{i} + exp(-2));
end
C     = C/sum(C);
D     = C/sum(D);
lnD   = log(D + eps);
 
% effector action (e) and sufficient statistics of proposal density
%--------------------------------------------------------------------------
a      = ones(Ns,T);                           % state probability
b      = ones(Nu,T);                           % control probability
a(:,1) = S;                                    % first state
a(:,T) = C;                                    % final state
a      = a*diag(1./sum(a));
b      = b*diag(1./sum(b));
 
% posterior expectations (sates Q, control R) and action (E)
%----------------------------------------------------------------------
Q(:,1,:) = a;
R(:,1,:) = b;
E        = sparse(Nu,T);
 
% solve
%==========================================================================
for k = 1:(T - 1)
    
    % forward and backward passes at this time point
    %----------------------------------------------------------------------
    for i = 1:16
        for t = [(T - 1):-1:max(k,2) max(k,2):(T - 1)]
            
            % get data likelihood if available at this time
            %--------------------------------------------------------------
            try
                at = lnA'*S(:,t);
            catch
                at = sparse(Ns,1);
            end
            
            % and accumulate empirical priors
            %--------------------------------------------------------------
            for j = 1:Nu
                at = at + b(j,t + 1) * lnB{j} *a(:,t + 1);
                at = at + b(j,t    ) * lnB{j}'*a(:,t - 1);
            end
            for j = 1:Nu
                bt(j,1) = lnD(j);
                bt(j,1) = bt(j,1) + a(:,t - 1)'* lnB{j} *a(:,t);
            end
            
            % update sufficient statistics of hidden states and control
            %--------------------------------------------------------------
            at     = exp(at - max(at));
            bt     = exp(bt - max(bt));
            a(:,t) = at/sum(at);
            b(:,t) = bt/sum(bt);
            
        end
        
        % plot
        %------------------------------------------------------------------
        if PLOT
            subplot(2,1,1)
            imagesc(1 - a)
            title('Expected hidden states','FontSize',16)
            xlabel('time','FontSize',12)
            ylabel('hidden state','FontSize',12)
            
            subplot(2,2,3)
            imagesc(1 - b); axis square
            title('Expected control states','FontSize',16)
            xlabel('time','FontSize',12)
            ylabel('control state','FontSize',12)
            drawnow
        end
    end
    
    
    % sampling of next state (outcome)
    %======================================================================
    for i = 1:Nu
        F(i) = P{i}(:,find(S(:,k)))'*lnA*a(:,k + 1);
    end
    
    % next action (the action and minimises expected free energy)
    %----------------------------------------------------------------------
    Pu       = exp(F(:) - max(F));
    Pu       = Pu/sum(Pu);
    [~, i]   = max(Pu.*rand(Nu,1));
    
    % next state (assuming G mediates uncertainty modelled the likelihood)
    %----------------------------------------------------------------------
    Ps       = G{i}*S(:,k);
    Ps       = Ps/sum(Ps);
    [~, j]   = max(Ps.*rand(Ns,1));
    
    
    % save action, state and posterior expectations (states Q, control R)
    %----------------------------------------------------------------------
    Q(:,k + 1,:) = a;
    R(:,k + 1,:) = b;
    S(:,k + 1)   = sparse(j,1,1,Ns,1);
    E(:,k + 1)   = sparse(i,1,1,Nu,1);
    
    % plot actual action
    %----------------------------------------------------------------------
    if PLOT
        subplot(2,2,4)
        imagesc(1 - E); axis square
        title('Selected action','FontSize',16)
        xlabel('time','FontSize',12)
        ylabel('action','FontSize',12)
    end
    
end
