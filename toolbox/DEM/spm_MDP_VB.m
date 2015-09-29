function [MDP] = spm_MDP_VB(MDP,OPTIONS)
% active inference and learning using variational Bayes
% FORMAT [MDP] = spm_MDP_VB(MDP,OPTIONS)
%
% MDP.S(N,1)      - true initial state
% MDP.V(T - 1,P)  - P allowable policies (control sequences)
%
% MDP.A(O,N)      - likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - prior preferences   (prior over future states)
% MDP.D(N,1)      - prior probabilities (prior over initial states)
%
% MDP.a(O,N)      - concentration parameters for A
% MDP.b{M}(N,N)   - concentration parameters for B
% MDP.c(N,N)      - concentration parameters for H
% MDP.d(N,1)      - concentration parameters for D
% MDP.e(P,1)      - concentration parameters for u
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
% and hidden controls are isomorphic.
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
% the numerics or combinatorics considerably. For example, situations in
% which one action can be selected at one time can be reduced to T polices
% – with one (shift) control being emitted at all possible time points.
% This specification of polices simplifies the generative model, allowing a
% fairly exhaustive model of potential outcomes – eschewing a mean field
% approximation over successive control states. In brief, the agent encodes
% beliefs about hidden states in the past and in the future conditioned
% on each policy (and a non-sequential state-state policy called a
% habit). These conditional expectations are used to evaluate the (path
% integral) of free energy that then determines the prior over policies.
% This prior is used to create a predictive distribution over outcomes,
% which specifies the next action.
%
% In addition to state estimation and policy selection, the scheme also
% updates model parameters; including the state transition matrices,
% mapping to outcomes and the initial state. This is useful for learning
% the context. In addition, by observing its own behaviour, the agent will
% automatically learn habits. Finally, by observing policies chosen over
% trials, the agent develops prior expectations or beliefs about what it 
% will do. If these priors (over policies – that include the habit) render
% some policies unlikely (using an Ockham's window), they will not be
% evaluated.
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
% $Id: spm_MDP_VB.m 6564 2015-09-29 08:10:22Z karl $
 
 
% deal with a sequence of trials
%==========================================================================
 
% options
%--------------------------------------------------------------------------
try, OPTIONS.scheme; catch, OPTIONS.scheme = 'Free Energy'; end
try, OPTIONS.habit;  catch, OPTIONS.habit  = 0;             end
try, OPTIONS.plot;   catch, OPTIONS.plot   = 0;             end
 
 
% if there are multiple trials ensure that parameters are updated
%--------------------------------------------------------------------------
if length(MDP) > 1
    
    OPTS      = OPTIONS;
    OPTS.plot = 0;
    for i = 1:length(MDP)
        
        % update concentration parameters
        %------------------------------------------------------------------
        if i > 1
            try,  MDP(i).a = OUT(i - 1).a; end
            try,  MDP(i).b = OUT(i - 1).b; end
            try,  MDP(i).c = OUT(i - 1).c; end
            try,  MDP(i).d = OUT(i - 1).d; end
            try,  MDP(i).e = OUT(i - 1).e; end
        end
        
        % solve this trial
        %------------------------------------------------------------------
        OUT(i) = spm_MDP_VB(MDP(i),OPTS);
        
    end
    MDP = OUT;
    
    % plot summary statistics - over trials
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
 
% numbers of transitions, policies and states
%--------------------------------------------------------------------------
T   = size(MDP.V,1) + 1;            % number of transitions
Np  = size(MDP.V,2);                % number of allowable policies
Ns  = size(MDP.B{1},1);             % number of hidden states
Nu  = size(MDP.B,2);                % number of hidden controls
V   = MDP.V;                        % allowable policies (T - 1,Np)
p0  = exp(-8);                      % smallest probability
q0  = 1/16;                         % smallest probability
Nh  = Np + 1;                       % index of habit
 
% parameters of generative model and policies
%==========================================================================
 
% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A  = MDP.A + p0;
    No = size(MDP.A,1);            % number of outcomes
catch
    A  = speye(Ns,Ns) + p0;
    No = Ns;
end
A      = spm_norm(A);              % normalise
 
% parameters (concentration parameters): A
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    qA = MDP.a + q0;
    qA = psi(qA) - ones(Ns,1)*psi(sum(qA));
else
    qA = log(A);
end
 
% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:Nu
    
    B{i} = MDP.B{i} + p0;
    B{i} = spm_norm(B{i});
    
    % parameters (concentration parameters): B
    %----------------------------------------------------------------------
    if isfield(MDP,'b')
        pB{i} = MDP.b{i} + q0;
        sB{i} = spm_norm(pB{i} );
        rB{i} = spm_norm(pB{i}');
    else
        sB{i} = spm_norm(B{i} );
        rB{i} = spm_norm(B{i}');
    end

end
 
% parameters (concentration parameters) - Habits
%--------------------------------------------------------------------------
if ~isfield(MDP,'c')
    MDP.c = 0;
    for j = 1:Nu
        MDP.c = MDP.c + MDP.B{j};
    end
end
pH    = MDP.c + q0;
sH    = spm_norm(pH );
rH    = spm_norm(pH');

 
% priors over initial hidden states - concentration parameters
%--------------------------------------------------------------------------
try
    d = MDP.d + q0;
catch
    d = ones(Ns,1);
end
qD    = psi(d) - ones(Ns,1)*psi(sum(d));
 
% priors over policies - concentration parameters
%--------------------------------------------------------------------------
try
    e = MDP.e + q0;
catch
    e(1:Np,1) = 8;
    e(Nh)     = 1;
end
if ~OPTIONS.habit
    e(Nh)     = q0;
end
qE    = psi(e) - ones(Nh,1)*psi(sum(e));
 
% prior preferences (log probabilities) : C
%--------------------------------------------------------------------------
try
    C = MDP.C;
catch
    C = zeros(Ns,1);
end
C     = log(spm_softmax(C));
 
% assume constant preferences, if only final states are specified
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
 
    
% precision defaults
%--------------------------------------------------------------------------
try, alpha = MDP.alpha;  catch, alpha  = 2; end
try, beta  = MDP.beta;   catch, beta   = 1; end
 
 
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
x  = zeros(Ns,T,Nh) + 1/Ns;         % expectations of hidden states | policy
X  = zeros(Ns,T + 1);               % expectations of hidden states
u  = zeros(Nh,T - 1);               % expectations of policy
a  = zeros(1, T - 1);               % action (index)
 
 
% initialise priors over policies and display utility with expected states
%--------------------------------------------------------------------------
u(:,1)     = spm_softmax(qE);
X(:,T + 1) = spm_softmax(hA + C(:,end));
 
% expected rate parameter
%--------------------------------------------------------------------------
qbeta  = beta;                      % initialise rate parameter
W      = zeros(1, T) + alpha/qbeta; % posterior precision
 
% solve
%==========================================================================
Ni     = 16;                        % number of VB iterations
rt     = zeros(1,T);                % reaction times
xn     = zeros(Ni,Ns,T,T,Np);       % state updates
un     = zeros(Nh,T*Ni);            % policy updates
wn     = zeros(T*Ni,1);             % simulated DA responses
for t  = 1:T
    
    % processing time and reset
    %----------------------------------------------------------------------
    tstart = tic;
    x      = spm_softmax(log(x)/2);
    
    % Variational updates (hidden states) under sequential policies
    %======================================================================
    
    % retain allowable policies within Ockham's window
    %----------------------------------------------------------------------
    if t > 1
        p = find([u(1:Np,t - 1)' 1] > 1/32);
    else
        p = 1:(Nh);
    end
    
    F     = zeros(T,Nh) + 16;
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
                if j <= t, v = v - qA(o(j),:)'; end
                if j == 1, v = v + qx - qD;     end
                if k > Np
                    if j > 1, v = v + qx - log(sH*x(:,j - 1,k)); end
                    if j < T, v = v + qx - log(rH*x(:,j + 1,k)); end
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
                if dF > 1/128, F0 = F0 - dF; else,  break,  end
            else
                F0 = sum(F(:,k));
            end
            
        end 
    end
    
    
    % (negative path integral of) free energy of policies (Q)
    %======================================================================
    Q     = zeros(T,Nh);
    for k = p
        for j = (t + 1):T
            Q(j,k) = x(:,j,k)'*(C(:,j) + hA - log(x(:,j,k)));
        end
    end
    
    
    % variational updates - policies and precision
    %======================================================================
    F     = sum(F,1)';
    Q     = sum(Q,1)';
    for i = 1:Ni
        
        % policy (u)
        %------------------------------------------------------------------
        qu = spm_softmax(qE(p) +W(t)*Q(p) - F(p));
        pu = spm_softmax(qE(p) +W(t)*Q(p));
        
        % precision (W) with free energy gradients (v = -dF/dw)
        %------------------------------------------------------------------
        if isfield(MDP,'w')
            try
                W(t) = MDP.w(t);
            catch
                W(t) = MDP.w;
            end
        else
            v     = qbeta - beta + (qu - pu)'*Q(p);
            qbeta = qbeta - v/2;
            W(t)  = alpha/qbeta;
        end
        
        % simulated dopamine responses (precision as each iteration)
        %------------------------------------------------------------------
        u(p,t)  = qu;
        n       = (t - 1)*Ni + i;
        wn(n,1) = W(t);
        un(p,n) = qu;
        
    end
 
    
    % Bayesian model averaging of hidden states over policies
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
    
        % posterior expectations about action
        %==================================================================
        v     = log(A*X(:,t + 1));
        for j = 1:Nu
            qo     = A*B{j}*X(:,t);
            P(j,t) = (v - log(qo))'*qo;
        end
        
        % action selection
        %------------------------------------------------------------------
        P(:,t) = spm_softmax(8*P(:,t));
                
        % next action
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
        
        % next observed state
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
    
    % mapping from hidden states to outcomes: a
    %----------------------------------------------------------------------
    if isfield(MDP,'a')
        i        = MDP.a > 0;
        da       = O(:,t)*X(:,t)';
        MDP.a(i) = MDP.a + da(i);
    end
    
    % mapping from hidden states to hidden states: b(u)
    %----------------------------------------------------------------------
    if isfield(MDP,'b') && t > 1
        for k = 1:Np
            v           = V(t - 1,k);
            i           = MDP.b{v} > 0;
            db          = u(k,t - 1)*x(:,t,k)*x(:,t - 1,k)';
            MDP.b{v}(i) = MDP.b{v}(i) + db(i);
        end
    end
    
    % mapping from hidden states to hidden states - habit: c
    %----------------------------------------------------------------------
    if isfield(MDP,'c') && t > 1
        k        = Nh;
        i        = MDP.c > 0;
        dc       = x(:,t,k)*x(:,t - 1,k)';
        MDP.c(i) = MDP.c(i) + dc(i) - (MDP.c(i) - 1)/32;
    end
    
end
 
% initial hidden states: d
%--------------------------------------------------------------------------
if isfield(MDP,'d')
    i        = MDP.d > 0;
    MDP.d(i) = MDP.d(i) + X(i,1) - (MDP.d(i) - 1)/8;
end
 
% policies: e
%--------------------------------------------------------------------------
if isfield(MDP,'e')
    i         = MDP.e > 0;
    K(1:Np,1) = 8;        % forget all policies (time constant of 8)
    K(Nh)     = 16;       % habits decay more slowly
    MDP.e(i)  = MDP.e(i) + W(t)*u(i,T) - (MDP.e(i) - 1)./K(i);
end
 
% simulated dopamine responses
%--------------------------------------------------------------------------
dn    = gradient(wn) + wn/64;
 
% Bayesian model averaging of expected hidden states over policies
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
%==========================================================================
if OPTIONS.plot
    if ishandle(OPTIONS.plot)
        figure(OPTIONS.plot); clf
    else
        spm_figure('GetWin','MDP'); clf
    end
    spm_MDP_VB_trial(MDP)
end


function A = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A  = A*diag(1./sum(A,1));
