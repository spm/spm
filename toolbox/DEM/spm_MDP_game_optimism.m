function [MDP] = spm_MDP_game_optimism(MDP,OPTION)
% aaction selection using active inference (with optimism bias)
% FORMAT [MDP] = spm_MDP_game_optimism(MDP,OPTION)
%
% MDP.T           - process depth (the horizon)
% MDP.N           - number of variational iterations (default 4)
% MDP.S(N,1)      - true initial state
%
% MDP.A(O,N)      - Likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - terminal cost probabilities (prior over outcomes)
% MDP.D(N,1)      - initial prior probabilities (prior over hidden states)
%
% MDP.V(T,P)      - P allowable policies (control sequences over T times)
%
% optional:
% MDP.s(1 x T)    - vector of true states  - for deterministic solutions
% MDP.o(1 x T)    - vector of observations - for deterministic solutions
% MDP.a(1 x T)    - vector of action       - for deterministic solutions
% MDP.w(1 x T)    - vector of precisions   - for deterministic solutions
%
% MDP.B{T,M}(N,N) - model transition probabilities for each time point
% MDP.G{T,M}(N,N) - true  transition probabilities for each time point
%                   (default: MDP.G{T,M} = MDP.G{M} = MDP.B{M})
%
% MDP.plot        - switch to suppress graphics: (default: [0])
% MDP.alpha       - upper bound on precision (Gamma hyperprior – shape [8])
% MDP.beta        - precision over precision (Gamma hyperprior - rate  [1])
%
% produces:
%
% MDP.P(M,T)   - probability of emitting an action 1,...,M at time 1,...,T
% MDP.Q(N,T)   - an array of conditional (posterior) expectations over
%                N hidden states and time 1,...,T
% MDP.O(O,T)   - a sparse matrix of ones encoding outcomes at time 1,...,T
% MDP.S(N,T)   - a sparse matrix of ones encoding states at time 1,...,T
% MDP.U(M,T)   - a sparse matrix of ones encoding action at time 1,...,T
% MDP.W(1,T)   - posterior expectations of precision
% MDP.d        - simulated dopamine responses
%
% OPTION       - {'Free Energy' | 'KL Control' | 'Expected Utility'};
%
% This routine provides solutions of active inference (minimisation of
% variational free energy) using a generative model based upon a Markov
% decision process. This model and inference scheme is formulated
% in discrete space and time. This means that the generative model (and
% process) are  finite state machines or hidden Markov models whose
% dynamics are given by transition probabilities among states and the
% likelihood corresponds to a particular outcome conditioned upon
% hidden states. For simplicity, this routine assumes that action
% and hidden controls are isomorphic. If the dynamics of transition
% probabilities of the true process are not provided, this routine will use
% the equivalent probabilities from the generative model.
%
% This implementation equips agents with the prior beliefs that they will
% maximise expected free energy: expected free energy is the free energy
% of future outcomes under the posterior predictive distribution. This can
% be interpreted in several ways – most intuitively as minimising the KL
% divergence between predicted and preferred outcomes (specified as prior
% beliefs) – while simultaneously minimising the (predicted) entropy of 
% outcomes conditioned upon hidden states. Expected free energy therefore
% combines KL optimality based upon preferences or utility functions with
% epistemic value or information gain.
%
% This particular scheme is designed for any allowable policies or control
% sequences specified in MDP.V. Constraints on allowable policies can limit
% the numerics or combinatorics considerable. For example, situations in
% which one action can be selected at one time can be reduced to T polices
% – with one (shift) control being emitted at all possible time points.
% This specification of polices simplifies the generative model, allowing a
% fairly exhaustive model of potential outcomes – eschewing a mean field
% approximation over successive control states. In brief, the agent simply
% represents the current state and states in the immediate and distant
% future.
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
% probabilistic mapping between hidden states and outcomes
% into the transition probabilities G.
%
% See also:spm_MDP, which uses multiple future states and a mean field
% approximation for control states – but allows for different actions
% at all times (as in control problems).
%
% See also: spm_MDP_game_KL, which uses a very similar formulation but just
% maximises the KL divergence between the posterior predictive distribution
% over hidden states and those specified by preferences or prior beliefs.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_game_optimism.m 6077 2014-06-30 16:55:03Z spm $

% set up and preliminaries
%==========================================================================

% options and precision defaults
%--------------------------------------------------------------------------
try, PLOT  = MDP.plot;  catch, PLOT  = 0;             end
try, alpha = MDP.alpha; catch, alpha = 8;             end
try, beta  = MDP.beta;  catch, beta  = 8;             end
try, N     = MDP.N;     catch, N     = 4;             end
try, T     = MDP.T;     catch, T     = size(MDP.V,1); end

% options
%--------------------------------------------------------------------------
if nargin < 2, OPTION = 'Free Energy'; end
MDP.OPT = OPTION;


% set up figure if necessary
%--------------------------------------------------------------------------
if PLOT
    if ishandle(PLOT)
        figure(PLOT); clf
        PLOT = 2;
    else
        spm_figure('GetWin','MDP'); clf
    end
end

% generative model and initial states
%--------------------------------------------------------------------------
Ns    = size(MDP.B{1},1);         % number of hidden states
Nb    = size(MDP.B,1);            % number of time-dependent probabilities
Nu    = size(MDP.B,2);            % number of hidden controls
p0    = exp(-32);                 % smallest probability

% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A  = MDP.A + p0;
    No = size(MDP.A,1);           % number of outcomes
catch
    A  = speye(Ns,Ns) + p0;
    No = Ns;
end
A     = A*diag(1./sum(A));        % normalise
lnA   = log(A);                   % log probabilities
HA    = sum(A.*lnA);              % negentropy of observations

% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:T
    for j = 1:Nu
        if i == 1 || Nb == T
            B{i,j}   = MDP.B{i,j} + p0;
            B{i,j}   = B{i,j}*diag(1./sum(B{i,j}));
        else
            B{i,j}   = B{1,j};
        end
    end
end

% terminal probabilities (priors)
%--------------------------------------------------------------------------
try
    C = MDP.C + p0;
    
    % convert preferences over hidden states to priors over outcomes
    %----------------------------------------------------------------------
    if size(C,1) ~= No
        C = A*C;
    end
    
    % asume no preferences if only final outceoms are specifed
    %----------------------------------------------------------------------
    if size(C,2) ~= T
        C = C(:,end)*ones(1,T);
    end
    
catch
    C = ones(No,T);
end
C     = C*diag(1./sum(C));
lnC   = log(C);

% intital probabilities (priors)
%--------------------------------------------------------------------------
try
    D = spm_vec(MDP.D) + p0;
catch
    D = ones(Ns,1);
end
D     = D/sum(D);
lnD   = log(D);

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
            G{i,j} = G{i,j} + p0;
            G{i,j} = G{i,j}*diag(1./sum(G{i,j}));
        else
            G{i,j} = G{1,j};
        end
    end
end

% policies and their expectations
%--------------------------------------------------------------------------
V      = MDP.V;
Np     = size(V,2);                % number of allowable policies
w      = 1:Np;                     % indices of allowable policies
CA     = A'*C(:,T);                % for plotting preferred hidden states
CA     = CA/sum(CA);


% initial states and outcomes
%--------------------------------------------------------------------------
[p q]  = max(A*MDP.S(:,1));        % initial outcome (index)
s      = find( MDP.S(:,1));        % initial state   (index)
o      = sparse(1,1,q,1,T);        % observations    (index)
S      = sparse(s,1,1,Ns,T);       % states sampled  (1 in K vector)
O      = sparse(q,1,1,No,T);       % states observed (1 in K vector)
U      = zeros(Nu,T);              % action selected (1 in K vector)
P      = zeros(Nu,T);              % posterior beliefs about control
x      = zeros(Ns,T);              % expectations of hidden states
u      = zeros(Np,T);              % expectations of hidden states
a      = zeros(1, T);              % action (index)
W      = zeros(1, T);              % posterior precision


% solve
%==========================================================================
gamma  = [];                       % simulated dopamine responses
for t  = 1:T
    
    
    % allowable policies
    %----------------------------------------------------------------------
    if t > 1
        
        % retain allowable policies (that are consistent with last action)
        %------------------------------------------------------------------
        j = ismember(V(t - 1,:),a(t - 1));
        V = V(:,j);
        w = w(j);
        
        % update policy expectations
        %------------------------------------------------------------------
        u(w,t) = u(w,t - 1)/sum(u(w,t - 1));
        x(:,t) = B{t - 1,a(t - 1)}*x(:,t - 1);
        
    else
        
        u(:,t) = ones(Np,1)/Np;
        x(:,t) = spm_softmax(lnD);
        
    end
    
    
    % Variational iterations (assuming precise inference about past action)
    %======================================================================
    Q     = zeros(size(V,2),Ns);
    for i = 1:N
        
        % current state (x)
        %------------------------------------------------------------------
        if t > 1
            v  = log(B{t - 1,a(t - 1)}*x(:,t - 1));
        else
            v  = lnD;
        end
        v      = v + lnA(o(t),:)' + W(t)*Q'*u(w,t);
        x(:,t) = x(:,t)/2 + spm_softmax(v)/2;
        
        
        % value of policies x current state (under allowable policies)
        %==================================================================
        Q     = zeros(size(V,2),Ns);
        for k = 1:size(V,2)
            
            % path integral of expected free energy
            %--------------------------------------------------------------
            Bj    = 1;
            for j = t:T
                
                % transition probaility from current state
                %----------------------------------------------------------
                Bj  = B{j,V(j,k)}*Bj;
                ABj = A*Bj;
                
                % predicted entropy and divergence
                %----------------------------------------------------------
                switch OPTION
                    case{'Free Energy','FE'}
                        Q(k,:) = Q(k,:) + HA*Bj + (lnC(:,j) - log(ABj*x(:,t)))'*ABj;
                        
                    case{'KL Control','KL'}
                        Q(k,:) = Q(k,:) + (lnC(:,j) - log(ABj*x(:,t)))'*ABj;
                        
                    case{'Expected Utility','EU','RL'}
                        Q(k,:) = Q(k,:) + lnC(:,j)'*ABj;
                        
                    otherwise
                        disp(['unkown optiion: ' OPTION])
                end 
            end    
        end
        
        
        % precision (W)
        %------------------------------------------------------------------
        if isfield(MDP,'w')
            W(t) = MDP.w(t);
        else
            W(t) = alpha/(beta - u(w,t)'*Q*x(:,t));
        end

        % policy (u)
        %------------------------------------------------------------------
        u(w,t)   = spm_softmax(W(t)*Q*x(:,t));
        
        
        % simulated dopamine responses (precision as each iteration)
        %------------------------------------------------------------------
        gamma(end + 1,1) = W(t);
        
    end
    
    % posterior expectations (control)
    %======================================================================
    for j = 1:Nu
        for k = t:T
            P(j,k) = sum(u(w(ismember(V(k,:),j)),t));
        end
    end
    
    % next action (the action that minimises expected free energy)
    %------------------------------------------------------------------
    try
        a(t) = MDP.a(t);
    catch
        try
            a(t) = find(rand < cumsum(P(:,t)),1);
        catch
            error('there are no more allowable policies')
        end
    end
    
    % save action
    %------------------------------------------------------------------
    U(a(t),t) = 1;
    
    
    % sampling of next state (outcome)
    %======================================================================
    if t < T
        
        % next sampled state
        %------------------------------------------------------------------
        try
            s(t + 1) = MDP.s(t + 1);
        catch
            s(t + 1) = find(rand < cumsum(G{t,a(t)}(:,s(t))),1);
        end
        
        % next obsverved state
        %------------------------------------------------------------------
        try
            o(t + 1) = MDP.o(t + 1);
        catch
            o(t + 1) = find(rand < cumsum(A(:,s(t + 1))),1);
        end
        
        % save outcome and state sampled
        %------------------------------------------------------------------
        W(1,t + 1)        = W(t);
        O(o(t + 1),t + 1) = 1;
        S(s(t + 1),t + 1) = 1;
        
    end
    
    
    % plot
    %======================================================================
    if PLOT > 0
        
        % posterior beliefs about hidden states
        %------------------------------------------------------------------
        subplot(4,2,1)
        imagesc(1 - [x CA*max(max(x))/max(CA)])
        if size(x,1) > 128
            hold on, spm_spy(x,16,1), hold off
        end
        title('Inferred states (and utility)','FontSize',14)
        xlabel('Time','FontSize',12)
        ylabel('Hidden state','FontSize',12)
        
        
        % posterior beliefs about control states
        %==================================================================
        subplot(4,2,2)
        
        % make previous plots dotted lines
        %------------------------------------------------------------------
        if T > 2
            h     = get(gca,'Children'); hold on
            for i = 1:length(h)
                set(h(i),'LineStyle',':');
            end
            plot(P')
            title('Inferred action','FontSize',14)
            xlabel('Time','FontSize',12)
            ylabel('Control state','FontSize',12)
            spm_axis tight
        else
            bar(P)
            title('Inferred action','FontSize',14)
            xlabel('Contol state','FontSize',12)
            ylabel('Posterior expectation','FontSize',12)
        end
        
        
        % policies
        %------------------------------------------------------------------
        subplot(4,2,3)
        imagesc(MDP.V')
        title('Allowable policies','FontSize',14)
        ylabel('Policy','FontSize',12)
        xlabel('Time','FontSize',12)
        
        % expectations over policies
        %------------------------------------------------------------------
        subplot(4,2,4)
        image((1 - u)*64)
        title('Posterior probability','FontSize',14)
        ylabel('Policy','FontSize',12)
        xlabel('Time','FontSize',12)
        
        % true state (outcome)
        %------------------------------------------------------------------
        subplot(4,2,5)
        if size(S,1) > 128
            spm_spy(S,16)
        else
            imagesc(1 - S)
        end
        title('True states','FontSize',14)
        ylabel('State','FontSize',12)
        
        % sample (observation)
        %------------------------------------------------------------------
        subplot(4,2,7)
        if size(O,1) > 128
            spm_spy(O,16,1)
        else
            imagesc(1 - O)
        end
        title('Observed states','FontSize',14)
        xlabel('Time','FontSize',12)
        ylabel('State','FontSize',12)
        
        
        % action sampled (selected)
        %------------------------------------------------------------------
        subplot(4,2,6)
        if size(U,1) > 128
            spm_spy(U,16,1)
        else
            imagesc(1 - U)
        end
        title('Selected action','FontSize',14)
        ylabel('Action','FontSize',12)
        
        % expected action
        %------------------------------------------------------------------
        subplot(4,2,8)
        plot((1:length(gamma))/N,gamma)
        title('Expected precision (confidence)','FontSize',14)
        xlabel('Time','FontSize',12)
        ylabel('Precision','FontSize',12)
        spm_axis tight
        drawnow
        
    end
    
end

% deconvolve to simulate dopamine responses
%--------------------------------------------------------------------------
da     = pinv( tril(toeplitz(exp(-((1:length(gamma)) - 1)'/N))) )*gamma;

% assemble results and place in NDP structure
%--------------------------------------------------------------------------
MDP.P  = P;              % probability of action at time 1,...,T - 1
MDP.Q  = x;              % conditional expectations over N hidden states
MDP.O  = O;              % a sparse matrix, encoding outcomes at 1,...,T
MDP.S  = S;              % a sparse matrix, encoding the states
MDP.U  = U;              % a sparse matrix, encoding the action
MDP.W  = W;              % posterior expectations of precision
MDP.d  = gamma;          % simulated dopamine responses
MDP.da = da;             % simulated dopamine responses (deconvolved)


