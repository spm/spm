function [Q,R,S,U,P] = spm_MDP(MDP)
% solves the active inference problem for Markov decision processes
% FROMAT [Q,R,S,U,P] = spm_MDP(MDP)
%
% MDP.T           - process depth (the horizon)
% MDP.S(N,1)      - initial state
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - terminal cost probabilities (prior N over hidden states)
% MDP.D(M,1)      - control probabilities (prior over M control states)
%
% optional:
%
% MDP.W           - log-precision of beliefs about transitions (default: 1)
% MDP.G{M}(N,N)   - transition probabilities used to generate outcomes
%                   (default: the prior transition probabilities)
% MDP.A(N,N)      - Likelihood of outcomes given hidden states
%                   (default: an identity mapping from states to outcomes)
% MDP.B{T,M}(N,N) - transition probabilities for each time point
% MDP.G{T,M}(N,N) - transition probabilities for each time point
%                   (default: MDP.B{T,M} = MDP.B{M})
%
% MDP.plot        -  swtich to suppress graphics
%
% produces:
%
% Q(N,K,T) - an array of conditional (posterior) expectations over N hidden
%            states and time 1,...,T at time 1,...,K
% R(M,K,T) - an array of conditional expectations over M control
%            states and time 1,...,T at time 1,...,K
% S(N,T)   - a sparse matrix of ones, encoding the state at time 1,...,T
% U(M,T)   - a sparse matrix of ones, encoding the action at time 1,...,T
% P(M,T)   - probabaility of emitting action 1,...,M at time 1,...,T
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
%
% See also spm_MDP_game
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP.m 6061 2014-06-21 09:02:42Z karl $

% set up and preliminaries
%==========================================================================

% plotting and precision defaults
%--------------------------------------------------------------------------
try PLOT   = MDP.plot;   catch, PLOT   = 1; end
try lambda = MDP.lambda; catch, lambda = 1; end
try W      = MDP.W;      catch, W      = 1; end

% get figure
%--------------------------------------------------------------------------
if PLOT, spm_figure('GetWin','MDP'); clf,   end

% generative model and initial states
%--------------------------------------------------------------------------
P0    = exp(-32);             % smallest probability
T     = MDP.T;                % process depth (the horizon)
Ns    = size(MDP.B{1},1);     % number of hidden states
Nb    = size(MDP.B,1);        % number of time-dependent probabilities
Nu    = size(MDP.B,2);        % number of hidden controls


% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A = MDP.A + P0;
catch
    A = speye(Ns,Ns) + P0;
end
A     = A*diag(1./sum(A));

% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:T
    for j = 1:Nu
        if i == 1 || Nb == T
            B{i,j}   = MDP.B{i,j} + P0;
            H{i,j}   = B{i,j}';
            B{i,j}   = B{i,j}*diag(1./sum(B{i,j}));
            H{i,j}   = H{i,j}*diag(1./sum(H{i,j}));
            lnB{i,j} = log(B{i,j})*W;
            lnH{i,j} = log(H{i,j})*W;
        else
            B{i,j}   = B{1,j};
            lnB{i,j} = lnB{1,j};
            lnH{i,j} = lnH{1,j};
        end
    end
end

% terminal cost probabilities (priors)
%--------------------------------------------------------------------------
try
    C = spm_vec(MDP.C) + P0;
catch
    C = ones(Ns,1);
end

% control probabilities (priors)
%--------------------------------------------------------------------------
try
    D = spm_vec(MDP.D) + P0;
catch
    D = ones(Nu,1);
end

% state probabilities (priors)
%--------------------------------------------------------------------------
try
    E = spm_vec(MDP.E) + P0;
catch
    E = ones(Ns,1);
end

C     = C/sum(C);
D     = D/sum(D);
E     = E/sum(E);
lnA   = log(A);
lnD   = log(D)*W;
lnE   = log(E)*W;

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


% effector action (e) and sufficient statistics of proposal density
%--------------------------------------------------------------------------
a        = zeros(Ns,T);                          % state probability
b        = ones(Nu,T);                           % control probability
a(:,T)   = C;                                    % final state
b        = b*diag(1./sum(b));

% posterior expectations (states Q, control R)
%--------------------------------------------------------------------------
Q(:,1,:) = a;
R(:,1,:) = b;

%  initial state and posterior (states Q, control R) and action (E)
%--------------------------------------------------------------------------
s     = find(spm_vec(MDP.S));
S     = sparse(s,1,1,Ns,T);
U     = sparse(Nu,T);

% solve
%==========================================================================
for k = 1:(T - 1)
    
    % forward and backward passes at this time point
    %----------------------------------------------------------------------
    for i = 1:8
        for t = (T - 1):-1:k
            
            % get data likelihood if available at this time
            %--------------------------------------------------------------
            if t > k
                at = lnE;
            else
                at = lnA'*S(:,t) + lnE;
            end
            
            % and accumulate empirical priors
            %--------------------------------------------------------------
            for j = 1:Nu
                if t > 1
                    at  = at     + b(j,t - 1) *lnH{t,j}'*a(:,t - 1);
                end
                at      = at     + b(j,t    ) *lnB{t,j}'*a(:,t + 1);
                bt(j,1) = lnD(j) + a(:,t    )'*lnB{t,j}'*a(:,t + 1);
            end
            
            % update sufficient statistics of hidden states and control
            %--------------------------------------------------------------
            a(:,t) = spm_softmax(at);
            b(:,t) = spm_softmax(bt);
            
            % graphics to inspect update scheme
            %==============================================================
            if PLOT > 2
                
                % posterior beliefs about hidden states
                %----------------------------------------------------------
                spm_plot_states(a,b)
                
                % pause if requested
                %----------------------------------------------------------
                if PLOT > 3, pause, end
                
            end
            
            
        end
        
        % graphics to inspect update scheme
        %==================================================================
        if PLOT > 1
            
            % posterior beliefs about hidden states
            %--------------------------------------------------------------
            spm_plot_states(a,b)
            
        end
    end
    
    
    % sampling of next state (outcome)
    %======================================================================
    for i = 1:Nu
        F(i,1) = B{k,i}(:,s)'*lnA*a(:,k + 1);
    end
    
    % next action (the action and minimises expected free energy)
    %----------------------------------------------------------------------
    Pu       = spm_softmax(F,lambda);
    i        = find(rand < cumsum(Pu),1);
    
    % next state (assuming G mediates uncertainty modelled the likelihood)
    %----------------------------------------------------------------------
    Ps       = spm_softmax(G{k,i}(:,s),lambda);
    s        = find(rand < cumsum(Ps),1);
    
    
    % save action, state and posterior expectations (states Q, control R)
    %----------------------------------------------------------------------
    P(:,k)     = Pu;
    Q(:,k,:)   = a;
    R(:,k,:)   = b;
    S(s,k + 1) = 1;
    U(i,k)     = 1;
    
    
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
if size(a,1) > 512
    spy(a > 1/32,16)
else
    imagesc(1 - a)
end
axis square
title('Expected hidden states','FontSize',16)
xlabel('time','FontSize',12)
ylabel('hidden state','FontSize',12)

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


