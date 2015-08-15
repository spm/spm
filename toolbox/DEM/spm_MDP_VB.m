function [MDP] = spm_MDP_VB(MDP,OPTIONS)
% action selection using active inference
% FORMAT [MDP] = spm_MDP_VB(MDP,OPTIONS)
%
% MDP.S(N,1)      - true initial state
% MDP.V(T - 1,P)  - P allowable policies (control sequences)
%
% MDP.A(O,N)      - Likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - prior preferences   (prior over future states)
% MDP.D(N,1)      - prior probabilities (prior over initial states)
%
% MDP.a(O,N)      - concentration parameters for A
% MDP.b{M}(N,N)   - concentration parameters for B
% MDP.h(N,N)      - concentration parameters for H
% MDP.d(N,1)      - concentration parameters for D
%
% optional:
% MDP.s(1 x T)    - vector of true states  - for deterministic solutions
% MDP.o(1 x T)    - vector of observations - for deterministic solutions
% MDP.u(1 x T)    - vector of action       - for deterministic solutions
% MDP.w(1 x T)    - vector of precisions   - for deterministic solutions
%
% MDP.alpha       - upper bound on precision (Gamma hyperprior – shape [8])
% MDP.beta        - precision over precision (Gamma hyperprior - rate  [1])
%
% OPTIONS.plot    - switch to suppress graphics: (default: [0])
% OPTIONS.scheme  - {'Free Energy' | 'KL Control' | 'Expected Utility'};
% OPTIONS.habit   - switch to suppress habit learning: (default: [1])
%
%
% produces:
%
% MDP.P(M,T)      - probability of emitting action 1,...,M at time 1,...,T
% MDP.Q(N,T)      - an array of conditional (posterior) expectations over
%                   N hidden states and time 1,...,T
% MDP.X           - and Bayesian model averages over policies
% MDP.R           - conditional expectations over policies
% MDP.O(O,T)      - a sparse matrix encoding outcomes at time 1,...,T
% MDP.S(N,T)      - a sparse matrix encoding states at time 1,...,T
% MDP.U(M,T)      - a sparse matrix encoding action at time 1,...,T
% MDP.W(1,T)      - posterior expectations of precision
%
% MDP.un  = un;   - simulated neuronal encoding of hidden states
% MDP.xn  = Xn;   - simulated neuronal encoding of policies
% MDP.wn  = wn;   - simulated neuronal encoding of precision
% MDP.dn  = dn;   - simulated dopamine responses (deconvolved)
% MDP.rt  = rt;   - simulated dopamine responses (deconvolved)
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
% $Id: spm_MDP_VB.m 6523 2015-08-15 20:57:28Z karl $


% deal with a sequence of trials
%==========================================================================

% options
%--------------------------------------------------------------------------
try, OPTIONS.scheme; catch, OPTIONS.scheme = 'Free Energy'; end
try, OPTIONS.habit;  catch, OPTIONS.habit  = 1;             end
try, OPTIONS.plot;   catch, OPTIONS.plot   = 0;             end

if length(MDP) > 1
    for i = 1:length(MDP)
        
        % transfer concentration parameters
        %------------------------------------------------------------------
        if i > 1
            try,  MDP(i).a = OUT(i - 1).a; end
            try,  MDP(i).b = OUT(i - 1).b; end
            try,  MDP(i).h = OUT(i - 1).h; end
            try,  MDP(i).d = OUT(i - 1).d; end
        end
        
        % integrate this trial
        %------------------------------------------------------------------
        OPTS = OPTIONS;
        if i < length(MDP)
            OPTS.plot = 0;
        else
            OPTS.plot = 'plot';
        end
        OUT(i) = spm_MDP_VB(MDP(i),OPTS);
        
    end
    
    MDP = OUT;
    
    % summary statistics -over trials (set up figure if necessary)
    %----------------------------------------------------------------------
    if OPTIONS.plot
        if ishandle(OPTIONS.plot)
            figure(OPTIONS.plot); clf
        else
            spm_figure('GetWin','MDP'); clf
        end
        spm_MDP_VB_game(MDP)
    end
    return
end


% set up and preliminaries
%==========================================================================

% options and precision defaults
%--------------------------------------------------------------------------
try, alpha = MDP.alpha;  catch, alpha  = 8;             end
try, beta  = MDP.beta;   catch, beta   = 2;             end

% generative model and initial states
%--------------------------------------------------------------------------
T      = size(MDP.V,1) + 1;        % number of transitions
Ns     = size(MDP.B{1},1);         % number of hidden states
Nu     = size(MDP.B,2);            % number of hidden controls
p0     = exp(-8);                  % smallest probability

% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A  = MDP.A + p0;
    No = size(MDP.A,1);            % number of outcomes
catch
    A  = speye(Ns,Ns) + p0;
    No = Ns;
end
A      = A*diag(1./sum(A));        % normalise

% parameters (concentration parameters) - A
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    qA = psi(MDP.a) - ones(Ns,1)*psi(sum(MDP.a));
else
    qA = log(A);
end

% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:Nu
    
    B{i} = MDP.B{i} + p0;
    B{i} = B{i}*diag(1./sum(B{i}));
    
    % parameters (concentration parameters) - B
    %----------------------------------------------------------------------
    if isfield(MDP,'b')
        qB{i} = psi(MDP.b{i}) - ones(Ns,1)*psi(sum(MDP.b{i}));
        sB{i} = spm_softmax(qB{i});
        rB{i} = spm_softmax(log(B{i})');
    else
        qB{i} = log(B{i});
        sB{i} = B{i};
        rB{i} = spm_softmax(log(B{i})');
    end
end

% initial probabilities (priors)
%--------------------------------------------------------------------------
try
    d = MDP.d;
catch
    d = ones(Ns,1);
end
qD    = psi(d) - ones(Ns,1)*psi(sum(d));

% parameters (concentration parameters) - Habits
%--------------------------------------------------------------------------
try
    h = MDP.h;
catch
    h = [];
end
if isempty(h)
    h = 1;
    for j = 1:Nu
        h = h + B{j}*2;
    end
    MDP.h = h;
end
qH    = psi(h) - ones(Ns,1)*psi(sum(h));

% terminal probabilities (priors)
%--------------------------------------------------------------------------
try
    C = MDP.C;
catch
    C = zeros(Ns,1);
end
C     = log(spm_softmax(C));

% asume constant preferences if only final states are specified
%--------------------------------------------------------------------------
if size(C,2) ~= T
    C = C(:,end)*ones(1,T);
end

% OPTIONS
%--------------------------------------------------------------------------
switch OPTIONS.scheme
    case{'Free Energy','FE'}
        hA = sum(spm_softmax(qA).*qA)';
        
    case{'KL Control','KL','Expected Utility','EU','RL'}
        hA = 0;

    otherwise
        disp(['unkown option: ' OPTIONS])
end

% policies and their expectations
%--------------------------------------------------------------------------
V  = MDP.V;                         % allowable policies (T - 1,Np)
Np = size(V,2);                     % number of allowable policies
N  = 8;                             % iterations of precision

% initial states and outcomes
%--------------------------------------------------------------------------
try
    s  = MDP.s(1);                  % initial state   (index)
catch
    s  = find(MDP.S(:,1));          % initial state   (index)
end
q  = find(rand < cumsum(A(:,s)),1); % initial outcome (index)

o  = sparse(1,1,q,1,T);             % observations    (index)
S  = sparse(s,1,1,Ns,T);            % states sampled  (1 in K vector)
O  = sparse(q,1,1,No,T);            % states observed (1 in K vector)
U  = zeros(Nu,T - 1);               % action selected (1 in K vector)
P  = zeros(Nu,T - 1);               % posterior beliefs about control
x  = zeros(Ns,T,Np + 1) + 1/Ns;     % expectations of hidden states | policy
X  = zeros(Ns,T + 1);               % expectations of hidden states
u  = zeros(Np + 1,T - 1);           % expectations of policy
a  = zeros(1, T - 1);               % action (index)
W  = zeros(1, T);                   % posterior precision

% add utility to expected sates (for display purposes)
%--------------------------------------------------------------------------
X(:,T + 1) = spm_softmax(hA + C(:,end));


% solve
%==========================================================================
Ni     = 16;                        % number of VB iterations
rt     = zeros(1,T);                % reaction times
xn     = zeros(Ni,Ns,T,T,Np);       % state updates
un     = zeros(Np + 1,T*N);         % policy updates
wn     = zeros(T*N,1);              % simulated DA responses
qbeta  = beta;                      % expected rate parameter
for t  = 1:T
    
    % processing time and reset
    %----------------------------------------------------------------------
    tstart = tic;
    x      = spm_softmax(log(x)/2);
    
    
    % Variational updates (hidden states) under sequential policies
    %======================================================================
    
    % retain allowable policies (that are consistent with last action)
    %----------------------------------------------------------------------
    if t > 1
        p = find([u(1:Np,t - 1)' 1] > exp(-4));
    else
        p = 1:(Np + 1);
    end
    
    for k = p
        
        % gradient descent on free energy
        %------------------------------------------------------------------
        for i = 1:Ni
            for j = 1:T
                
                % current state
                %----------------------------------------------------------
                xj   = x(:,j,k);
                qx   = log(xj);
                
                % evaluate free energy and gradients (v = dFdx)
                %----------------------------------------------------------
                v    = 0;
                if j <= t, v = v - qA(o(j),:)';              end
                if j == 1, v = v + qx - qD;                  end
                if k > Np
                    if j > 1, v = v + qx - qH* x(:,j - 1,k); end
                    if j < T, v = v      - qH'*x(:,j + 1,k); end
                else
                    if j > 1, v = v + qx - log(sB{V(j - 1,k)}*x(:,j - 1,k)); end
                    if j < T, v = v + qx - log(rB{V(j    ,k)}*x(:,j + 1,k)); end
                end
                
                % update
                %----------------------------------------------------------
                x(:,j,k) = spm_softmax(qx - v/4);
                F(j,k)   = xj'*v;
                
                % record neuronal activity
                %----------------------------------------------------------
                xn(i,:,j,t,k) = x(:,j,k) - xj;

            end
            
            % convergence
            %--------------------------------------------------------------
            if i > 1
                dF = F0 - sum(F(:,k));
                if dF > 1/128
                    F0 = F0 - dF; 
                else
                    break
                end
            else
                F0 = sum(F(:,k));
            end
            
        end 
    end
    
    
    % expected (negative) free energy of policies (Q)
    %======================================================================
    Q     = zeros(Np + 1,1);
    if OPTIONS.habit
        Q(Np + 1) = 0;
    else
        Q(Np + 1) = -4;
    end
    for k = p
        
        % path integral of expected free energy
        %------------------------------------------------------------------
        for j = 1:T
            if j > t
                Q(k) = Q(k) + x(:,j,k)'*(C(:,j) + hA - log(x(:,j,k)));
            else
                Q(k) = Q(k) - F(j,k);
            end       
        end
    end
    
    
    % variational updates - policies and precision
    %======================================================================
    for i = 1:N
        
        % policy (u)
        %------------------------------------------------------------------
        v      = W(t)*Q(p);
        u(p,t) = spm_softmax(v);

        % precision (W)
        %------------------------------------------------------------------
        if isfield(MDP,'w')
            try
                W(t) = MDP.w(t);
            catch
                W(t) = MDP.w;
            end
        else
            v     = beta - u(p,t)'*Q(p) - qbeta;
            qbeta = qbeta + v*3/4;
            W(t)  = alpha/qbeta;
        end
        
        
        % simulated dopamine responses (precision as each iteration)
        %------------------------------------------------------------------
        n       = (t - 1)*N + i;
        wn(n,1) = W(t);
        un(p,n) = u(p,t);
        
    end
    
    % Baysian model averaging of hidden states over policies
    %----------------------------------------------------------------------
    for i = 1:T
        X(:,i) = squeeze(x(:,i,:))*u(:,t);
    end
    
    % processing time
    %----------------------------------------------------------------------
    rt(t) = toc(tstart);
    
    
    % action selection and sampling of next state (outcome)
    %======================================================================
    if t < T
    
        % posterior expectations (control)
        %==================================================================
        v     = log(A*X(:,t + 1));
        for j = 1:Nu
            qo     = A*B{j}*X(:,t);
            P(j,t) = (v - log(qo))'*qo;
        end
        
        % accommodate absorbing states (where action is undefined)
        %------------------------------------------------------------------
        P(:,t) = spm_softmax(8*P(:,t));
                
        % next action (the action that minimises expected free energy)
        %------------------------------------------------------------------
        try
            a(t) = MDP.u(t);
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
            s(t + 1) = find(rand < cumsum(B{a(t)}(:,s(t))),1);
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

end

% learning
%==========================================================================
for t = 1:T
    
    % mapping from hidden states to outcomes - A
    %----------------------------------------------------------------------
    if isfield(MDP,'a')
        MDP.a = MDP.a + O(:,t)*X(:,t)';
    end
    
    % mapping from hidden states to hidden states - B(u)
    %----------------------------------------------------------------------
    if isfield(MDP,'b') && t > 1
        for k = 1:Np
            v        = V(t - 1,k);
            MDP.b{v} = MDP.b{v} + u(k,t - 1)*x(:,t,k)*x(:,t - 1,k)';
        end
    end
    
    % mapping from hidden states to hidden states - H
    %----------------------------------------------------------------------
    if isfield(MDP,'h') && t > 1
        k     = Np + 1;
        MDP.h = MDP.h + x(:,t,k)*x(:,t - 1,k)';
    end
    
end

% initial hidden states - D
%--------------------------------------------------------------------------
if isfield(MDP,'d')
    MDP.d = MDP.d + X(:,1);
end

% simulated dopamine responses
%--------------------------------------------------------------------------
dn    = gradient(wn) + wn/64;

% Baysian model averaging (selection) of hidden states over policies
%--------------------------------------------------------------------------
Xn    = zeros(Ni,Ns,T,T);
for i = 1:T
    for k = 1:Np
        Xn(:,:,:,i) = Xn(:,:,:,i) + xn(:,:,:,i,k)*u(k,i);
    end
end


% assemble results and place in NDP structure
%--------------------------------------------------------------------------
MDP.P   = P;              % probability of action at time 1,...,T - 1
MDP.Q   = x;              % conditional expectations over N hidden states
MDP.X   = X;              % Bayesian model averages
MDP.R   = u;              % conditional expectations over policies
MDP.O   = O;              % a sparse matrix, encoding outcomes at 1,...,T
MDP.S   = S;              % a sparse matrix, encoding the states
MDP.U   = U;              % a sparse matrix, encoding the action
MDP.W   = W;              % posterior expectations of precision
MDP.C   = C;              % utility

MDP.un  = un;             % simulated neuronal encoding of hidden states
MDP.xn  = Xn;             % simulated neuronal encoding of policies
MDP.wn  = wn;             % simulated neuronal encoding of precision
MDP.dn  = dn;             % simulated dopamine responses (deconvolved)
MDP.rt  = rt;             % simulated dopamine responses (deconvolved)

% plot
%======================================================================
if OPTIONS.plot
    if ishandle(OPTIONS.plot)
        figure(OPTIONS.plot); clf
    else
        spm_figure('GetWin','MDP'); clf
    end
    spm_MDP_VB_trial(MDP)
end

return
