function [Q,R,S,U,P] = spm_MDP_select(MDP)
% solves the active inference problem for action selection
% FROMAT [Q,R,S,U,P] = spm_MDP(MDP)
%
% MDP.T           - process depth (the horizon)
% MDP.S(N,1)      - initial state
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - terminal probabilities (prior N over hidden states)
%
% optional:
%
% MDP.W           - log-precision of beliefs about transitions (default: 16)
% MDP.G{M}(N,N)   - transition probabilities used to generate outcomes
%                   (default: the prior transition probabilities)
% MDP.A(N,N)      - Likelihood of outcomes given hidden states
%                   (default: an identity mapping from states to outcomes)
% MDP.B{T,M}(N,N) - transition probabilities for each time point
% MDP.G{T,M}(N,N) - transition probabilities for each time point
%                   (default: MDP.B{T,M} = MDP.B{M})
%
% MDP.plot        - switch to suppress graphics
%
% produces:
%
% Q(N,T)   - an array of conditional (posterior) expectations over N hidden
%            states and time 1,...,T
% R(M,T)   - an array of conditional expectations over M control
%            states and time 1,...,T
% S(N,T)   - a sparse matrix of ones, encoding the state at time 1,...,T
% U(M,T)   - a sparse matrix of ones, encoding the action at time 1,...,T
% P(M,T)   - probability of emitting action 1,...,M at time 1,...,T
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
 
% This particular scheme is designed for situations in which a particular
% action is to be selected at a particular time. The first control state
% is considered the baseline or default behaviour. Subsequent control
% states are assumed to be selected under the constraint that only one
% action can be emitted once. This provides a particular simplification for
% the generative model, allowing a fairly exhaustive model of potential
% outcomes – eschewing a mean field approximation over successive
% control states. In brief, the agent simply represents the current state
% and states in the immediate and distant future. 
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
% Partially observed Markov decision processes can be modelled by
% specifying a likelihood (as part of a generative model) and absorbing any
% probabilistic mapping between (isomorphic) hidden states and outcomes
% into the transition probabilities G.
%
% See also spm_MDP which uses multiple future states and a the mean
% field, approximation for control states – but allows for different actions
% at all times (as in control problems).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP.m 5067 2012-11-18 22:05:39Z karl $
 
% set up and preliminaries
%==========================================================================
 
% plotting and precision defaults
%--------------------------------------------------------------------------
try PLOT   = MDP.plot;   catch, PLOT   = 1;     end
try lambda = MDP.lambda; catch, lambda = 1;     end
try W      = MDP.W;      catch, W      = 1;   end
 
if PLOT, spm_figure('GetWin','MDP'); clf, end
 
% generative model and initial states
%--------------------------------------------------------------------------
T     = MDP.T;            % process depth (the horizon)
Ns    = size(MDP.B{1},1); % number of hidden states
Nb    = size(MDP.B,1);    % number of time-dependent probabilities
Nu    = size(MDP.B,2);    % number of hidden controls
 
 
% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A = MDP.A + exp(-16);
catch
    A = speye(Ns,Ns) + exp(-16);
end
A     = A*diag(1./sum(A));
lnA   = log(A);
 
 
% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:T
    for j = 1:Nu
        if i == 1 || Nb == T
            B{i,j}   = MDP.B{i,j} + exp(-16);
            B{i,j}   = B{i,j}*diag(1./sum(B{i,j}));
            lnB{i,j} = log(B{i,j})*W;
        else
            B{i,j}   = B{1,j};
            lnB{i,j} = lnB{1,j};
        end
    end
end
 
% terminal cost probabilities (priors)
%--------------------------------------------------------------------------
C     = spm_vec(MDP.C) + exp(-16);
C     = C/sum(C);
 
 
% generative process (assume the true process is the same as the model)
%--------------------------------------------------------------------------
try
    G = MDP.G;
catch
    G = MDP.B;
end
Ng    = size(G,1);
for i = 1:T
    for j = 1:Nu
        if i == 1 || Ng == T
            G{i,j}   = G{i,j}*diag(1./sum(G{i,j}));
        else
            G{i,j}   = G{1,j};
        end
    end
end
 
% initial states and outcomes
%--------------------------------------------------------------------------
S        = spm_vec(MDP.S);   % 
s        = find(S);
S        = sparse(s,1,1,Ns,T);
U        = sparse(Nu,T);
 
 
% solve
%==========================================================================
for t = 1:(T - 1)
    
    % sufficient statistics of hidden states
    %----------------------------------------------------------------------
    a      = zeros(Ns,3);
    a(:,3) = C;
    
    % transition probabilities from current to final state
    %----------------------------------------------------------------------
    for j = 2:Nu
        for k = t:(T - 1)
            H = 1;
            for f = t:(T - 1)
                if f == k
                    H = B{f,j}*H;
                else
                    H = B{f,1}*H;
                end
            end
            lnH{k,j}  = max(log(H),-16);
        end
    end
    
    
    % forward and backward passes at this time point
    %----------------------------------------------------------------------
    for i = 1:8
        
        % current state a(:,1)
        %------------------------------------------------------------------
        a(:,1) = lnA'*S(:,t);
        a(:,1) = spm_softmax(a(:,1));
                        
        % policy (b)
        %------------------------------------------------------------------
        for j = 2:Nu
            for k  = t:(T - 1)
                b(j,k) =              a(:,3)'*lnH{k,j}*a(:,1);
                if k == t
                    b(j,k) = b(j,k) + a(:,2)'*lnB{t,j}*a(:,1);
                else
                    b(j,k) = b(j,k) + a(:,2)'*lnB{t,1}*a(:,1);
                end
            end
        end
        j      = 2:Nu;
        k      = t:(T - 1);
        b(j,k) = spm_unvec(spm_softmax(spm_vec(b(j,k))),b(j,k));
        b(1,:) = 1 - sum(b(j,:),1);
        
        
        % expected utility policy (p)
        %------------------------------------------------------------------
        for k = t:(T - 1)
            for j = 2:Nu
                K(:,k,j - 1) =   exp(lnH{k,j})*a(:,1);
                p(j,k) = log(a(:,3))'*exp(lnH{k,j})*a(:,1);
            end
        end
        j      = 2:Nu;
        k      = t:(T - 1);
        p(j,k) = spm_unvec(spm_softmax(spm_vec(p(j,k)),32),p(j,k));
        p(1,:) = 1 - sum(p(j,:),1);
        
        % next state a(:,2)
        %------------------------------------------------------------------
        a(:,2) = zeros(Ns,1);
        for j  = 1:Nu
            a(:,2) = a(:,2) + b(j,t)*lnB{t,j}*a(:,1);
        end
        a(:,2) = spm_softmax(a(:,2));
        
        
        % graphics to inspect update scheme
        %==================================================================
        if PLOT > 2
            
            % posterior beliefs about hidden states
            %--------------------------------------------------------------
            spm_plot_states(a,b)
            
            % pause if requested
            %--------------------------------------------------------------
            if PLOT > 3, pause, end
            
        end
        
    end
        
    
    % graphics to inspect update scheme
    %======================================================================
    if PLOT > 1
        
        % posterior beliefs about hidden states
        %------------------------------------------------------------------
        spm_plot_states(a,b)
        
    end
    
    % sampling of next state (outcome)
    %======================================================================
    for i = 1:Nu
        F(i) = B{t,i}(:,s)'*lnA*a(:,2);
    end
    % F        = log(p(:,t));
    
    % next action (the action and minimises expected free energy)
    %----------------------------------------------------------------------
    Pu       = spm_softmax(F(:),lambda);
    i        = find(rand < cumsum(Pu),1);
    
    % next state (assuming G mediates uncertainty modelled the likelihood)
    %----------------------------------------------------------------------
    Ps       = G{t,i}(:,s);
    s        = find(rand < cumsum(Ps),1);
    
    
    % save action, state and posterior expectations (states Q, control R)
    %----------------------------------------------------------------------
    P(:,t)     = Pu;
    Q(:,t)     = full(a(:,2));
    R(:,t)     = full(b(:,t));
    S(s,t + 1) = 1;
    U(i,t)     = 1;
    
    
    % plot
    %======================================================================
    if PLOT > 0
        
        % posterior beliefs about hidden states
        %------------------------------------------------------------------
        spm_plot_states(a,b)
        
        % states sampled (outcome)
        %------------------------------------------------------------------
        subplot(2,2,3)
        if size(S,1) > 512
            spy(S,16)
        else
            imagesc(1 - S)
        end
        axis square
        title('Sampled state','FontSize',16)
        xlabel('time','FontSize',12)
        ylabel('action','FontSize',12)
        
        % action sampled (selected)
        %------------------------------------------------------------------
        subplot(2,2,4)
        if size(U,1) > 512
            spy(U,16)
        else
            imagesc(1 - U)
        end
        axis square
        title('Selected action','FontSize',16)
        xlabel('time','FontSize',12)
        ylabel('action','FontSize',12)
        drawnow
        
    end
    
end
 
 
function spm_plot_states(a,b)
 
% posterior beliefs about hidden states
%--------------------------------------------------------------------------
subplot(2,2,1)
bar(a)
axis square
title('Expected states','FontSize',16)
xlabel('time','FontSize',12)
ylabel('hidden state','FontSize',12)
legend({'now','next','last'})
 
% posterior beliefs about control states
%--------------------------------------------------------------------------
subplot(2,2,2)
if size(b,1) > 512
    spy(b > 1/8,16)
else
    imagesc(1 - b)
end
axis square
title('Expected control states','FontSize',16)
xlabel('time','FontSize',12)
ylabel('control state','FontSize',12)
drawnow
