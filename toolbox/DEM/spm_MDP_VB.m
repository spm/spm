function [MDP] = spm_MDP_VB(MDP,OPTION,W)
% action selection using active inference
% FORMAT [MDP] = spm_MDP_VB(MDP,OPTION,W)
%
% MDP.N           - number of variational iterations (default 4)
% MDP.S(N,1)      - true initial state
%
% MDP.A(O,N)      - Likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - prior preferences (prior over future states)
% MDP.D(N,1)      - prior probabilities (prior over initial states)
%
% MDP.V(T - 1,P)  - P allowable policies (control sequences)
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
% MDP.gamma       - initial precision
% MDP.lamba       - precision update rate
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
% MDP.da       - simulated dopamine responses (deconvolved)
% MDP.KLx      - updating as scored with KL (state estimation)
% MDP.KLu      - updating as scored with KL (policy selection)
%
% OPTION       - {'Free Energy' | 'KL Control' | 'Expected Utility'};
% W            - optional fixed precision
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
% $Id: spm_MDP_VB.m 6451 2015-05-26 09:26:03Z karl $

% set up and preliminaries
%==========================================================================

% options and precision defaults
%--------------------------------------------------------------------------
try, PLOT   = MDP.plot;   catch, PLOT   = 0;             end
try, alpha  = MDP.alpha;  catch, alpha  = 8;             end
try, beta   = MDP.beta;   catch, beta   = 4;             end
try, g      = MDP.gamma;  catch, g      = 1;             end
try, lambda = MDP.lambda; catch, lambda = 0;             end
try, N      = MDP.N;      catch, N      = 4;             end

% options and number of outcomes
%--------------------------------------------------------------------------
T       = size(MDP.V,1) + 1;
if nargin < 2, OPTION = 'Free Energy'; end
MDP.OPT = OPTION;
if nargin > 2
    MDP.w = zeros(1,T) + W;
end


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
Ns     = size(MDP.B{1},1);         % number of hidden states
Nb     = size(MDP.B,1);            % number of time-dependent probabilities
Nu     = size(MDP.B,2);            % number of hidden controls
p0     = exp(-16);                 % smallest probability

% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A  = MDP.A + p0;
    No = size(MDP.A,1);           % number of outcomes
catch
    A  = speye(Ns,Ns) + p0;
    No = Ns;
end
A      = A*diag(1./sum(A));        % normalise
lnA    = log(A);                   % log probabilities
H      = sum(A.*lnA)';             % negentropy of observations

% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:(T - 1)
    for j = 1:Nu
        if i == 1 || Nb == (T - 1)
            B{i,j} = MDP.B{i,j} + p0;
            B{i,j} = B{i,j}*diag(1./sum(B{i,j}));
        else
            B{i,j} = B{1,j};
        end
        lnB{i,j} = log(B{1,j});
    end
end

% terminal probabilities (priors)
%--------------------------------------------------------------------------
try
    C = MDP.C + p0;
    
    % asume constant preferences if only final states are specified
    %----------------------------------------------------------------------
    if size(C,2) ~= T
        C = C(:,end)*ones(1,T);
    end
    
catch
    C = ones(Ns,T);
end
C     = C*diag(1./sum(C));
lnC   = log(C);

% initial probabilities (priors)
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
for i = 1:(T - 1)
    for j = 1:Nu
        if i == 1 || Ng == (T - 1)
            G{i,j} = G{i,j} + p0;
            G{i,j} = G{i,j}*diag(1./sum(G{i,j}));
        else
            G{i,j} = G{1,j};
        end
    end
end

% pparameters (concentration parameters)
%--------------------------------------------------------------------------
try
    c  = MDP.c;
catch
    c  = 0;
    for j = 1:Nu
        c = c + B{1,j}/Nu;
    end
end
B0     = psi(c) - ones(Ns,1)*psi(sum(c));

% policies and their expectations
%--------------------------------------------------------------------------
V      = MDP.V;                    % allowable policies (T - 1,Np)
Np     = size(V,2);                % number of allowable policies
R      = ones(1,Np);               % policies in play

% initial states and outcomes
%--------------------------------------------------------------------------
[p q]  = max(A*MDP.S(:,1));        % initial outcome (index)
s      = find( MDP.S(:,1));        % initial state   (index)
o      = sparse(1,1,q,1,T);        % observations    (index)
S      = sparse(s,1,1,Ns,T);       % states sampled  (1 in K vector)
O      = sparse(q,1,1,No,T);       % states observed (1 in K vector)
U      = zeros(Nu,T - 1);          % action selected (1 in K vector)
P      = zeros(Nu,T - 1);          % posterior beliefs about control
x      = zeros(Ns,T,Np + 1);       % expectations of hidden states
X      = zeros(Ns,T + 1);          % expectations of hidden states
u      = zeros(Np + 1,T - 1);      % expectations of hidden states
a      = zeros(1, T - 1);          % action (index)
W      = zeros(1, T);              % posterior precision

X(:,T + 1) = spm_softmax(H + C(:,end));

% solve
%==========================================================================
KLx    = zeros(1,T);               % state updates
KLu    = zeros(1,T);               % policy updates
wn     = zeros(T*N,1);             % simulated DA responses
b      = alpha/g;                  % expected rate parameter
for t  = 1:T
    
    % Variational updates (hidden states) under habitual policies
    %======================================================================
    k     = Np + 1;
    for i = 1:4
        
        % past and future states
        %------------------------------------------------------------------
        for j = 1:T
            v     = 0;
            if j <= t, v = lnA(o(j),:)';         end
            if j > 1,  v = v + B0 *x(:,j - 1,k); end
            if j < T,  v = v + B0'*x(:,j + 1,k); end
            x(:,j,k) = spm_softmax(v);
        end
        
    end
    
    % learning
    %======================================================================
    if t > 1
        c = c + x(:,t,k)*x(:,t - 1,k)';
    end
    B0    = psi(c) - ones(Ns,1)*psi(sum(c));
    
    % Variational updates (hidden states) under sequential policies
    %======================================================================
    
    % retain allowable policies (that are consistent with last action)
    %----------------------------------------------------------------------
    if t > 1
        R = R & ismember(V(t - 1,:),a(t - 1));
    end
    w     = find(R);
    for k = w
        
        xk    = spm_vec(x(:,:,k)) + p0;
        for i = 1:4
            for j = 1:T
                v     = 0;
                if j <= t, v = lnA(o(j),:)';                            end
                if j > 1,  v = v + lnB{j - 1,V(j - 1,k)} *x(:,j - 1,k); end
                if j < T,  v = v + lnB{j,    V(j,    k)}'*x(:,j + 1,k); end
                x(:,j,k) = spm_softmax(v);
            end
        end
        
        % KL update – states
        %------------------------------------------------------------------
        kx     = spm_vec(x(:,:,k)) + p0;
        KLx(t) = KLx(t) + kx'*(log(kx) - log(xk));
        
    end
    
    % expected (negative) free energy of policies (Q)
    %======================================================================
    Q     = zeros(Np + 1,1);
    w     = [w Np + 1];
    for k = w
        
        % path integral of expected free energy
        %------------------------------------------------------------------
        for j = (t + 1):T
            
            switch OPTION
                case{'Free Energy','FE'}
                    v = lnC(:,j) - log(x(:,j,k) + p0) + H;
                    
                case{'KL Control','KL'}
                    v = lnC(:,j) - log(x(:,j,k) + p0);
                    
                case{'Expected Utility','EU','RL'}
                    v = lnC(:,j);
                    
                otherwise
                    disp(['unkown option: ' OPTION])
            end
            Q(k)   = Q(k) + v'*x(:,j,k);
            
        end
    end
    
    
    % variational updates - policies and precision
    %======================================================================
    for i = 1:N
        n = (t - 1)*N + i;
        
        % policy (u)
        %------------------------------------------------------------------
        v        = W(t)*Q(w);
        u(w,t)   = spm_softmax(v);
        un(w,n)  = u(w,t);
        
        % precision (W)
        %------------------------------------------------------------------
        if isfield(MDP,'w')
            W(t) = MDP.w(t);
        else
            b    = lambda*b + (1 - lambda)*(beta - u(w,t)'*Q(w));
            W(t) = alpha/b;
        end
        
        % simulated dopamine responses (precision as each iteration)
        %------------------------------------------------------------------
        wn(n,1)  = W(t);
        
    end
    
    % KL update – policies
    %----------------------------------------------------------------------
    if t > 1
        KLu(t) = u(:,t)'*(log(u(:,t) + p0) - log(u(:,t - 1) + p0));
    end
    
    % Baysian model averaging of hidden states over policies
    %----------------------------------------------------------------------
    for i = 1:T
        X(:,i) = squeeze(x(:,i,:))*u(:,t);
    end
    
    
    % action selection and sampling of next state (outcome)
    %======================================================================
    if t < T
        
        % posterior expectations (control)
        %==================================================================
        lnx   = log(X(:,t + 1));
        for j = 1:Nu
            xu     = B{t,j}*X(:,t);
            P(j,t) = (H + lnx - log(xu))'*xu;
        end
        P(:,t) = spm_softmax(P(:,t));
        
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
    if PLOT > 0 && (T > 3 || t == T)
        
        
        % posterior beliefs about hidden states
        %------------------------------------------------------------------
        subplot(4,2,1)
        imagesc(1 - X),hold on
        if size(X,1) > 128
            spm_spy(X,16,1)
        end
        plot(s,'.c','MarkerSize',16), hold off
        title('Hidden states (and utility)','FontSize',14)
        xlabel('trial','FontSize',12)
        ylabel('Hidden state','FontSize',12)
        
        % posterior beliefs about control states
        %==================================================================
        subplot(4,2,2)
        imagesc(1 - P), hold on
        plot(a,'.c','MarkerSize',16), hold off
        title('Inferred and selected action','FontSize',14)
        xlabel('trial','FontSize',12)
        ylabel('action','FontSize',12)
        
        % policies
        %------------------------------------------------------------------
        subplot(4,2,3)
        imagesc(MDP.V')
        set(gca,'YLim',[0 (Np + 1)] + 1/2)
        title('Allowable policies','FontSize',14)
        ylabel('policy','FontSize',12)
        xlabel('trial','FontSize',12)
        
        % expectations over policies
        %------------------------------------------------------------------
        subplot(4,2,4)
        imagesc(1 - un)
        title('Posterior probability','FontSize',14)
        ylabel('Policy','FontSize',12)
        xlabel('updates','FontSize',12)
        
        % sample (observation)
        %------------------------------------------------------------------
        subplot(4,2,5)
        if size(O,1) > 128
            spm_spy(O,16,1)
        else
            imagesc(1 - O)
        end
        title('Observed states','FontSize',14)
        xlabel('trial','FontSize',12)
        ylabel('outcome','FontSize',12)
        
        % expected action
        %------------------------------------------------------------------
        subplot(4,2,6)
        plot(wn,'k')
        title('Expected precision (dopamine)','FontSize',14)
        xlabel('updates','FontSize',12)
        ylabel('Precision','FontSize',12)
        drawnow
        
        % learned transition matrix
        %------------------------------------------------------------------
        subplot(4,1,4)
        image(spm_softmax(B0)*64)
        title('learned transitions','FontSize',14)
        xlabel('states','FontSize',12)
        ylabel('states','FontSize',12)
        axis square
        drawnow
        
    end
end

% simulated dopamine responses
%--------------------------------------------------------------------------
da  = gradient(wn) + wn/16;
if PLOT > 0
    subplot(4,2,6), hold on
    bar(4*da,'c'), plot(wn,'k'), hold off
    spm_axis tight
end


% assemble results and place in NDP structure
%--------------------------------------------------------------------------
MDP.P   = P;              % probability of action at time 1,...,T - 1
MDP.Q   = x;              % conditional expectations over N hidden states
MDP.O   = O;              % a sparse matrix, encoding outcomes at 1,...,T
MDP.S   = S;              % a sparse matrix, encoding the states
MDP.U   = U;              % a sparse matrix, encoding the action
MDP.W   = W;              % posterior expectations of precision
MDP.c   = c;              % concentration parameters of transitions
MDP.da  = da;             % simulated dopamine responses (deconvolved)
MDP.KLx = KLx;            % updating as scored with KL (state estimation)
MDP.KLu = KLu;            % updating as scored with KL (policy selection)

return



