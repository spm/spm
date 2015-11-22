function [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
% active inference and learning using variational Bayes
% FORMAT [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
%
% MDP.S(N,1)      - true initial state
% MDP.V(T - 1,P)  - P allowable policies (control sequences)
%
% MDP.A(O,N)      - likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - prior preferences   (prior over future outcomes)
% MDP.D(N,1)      - prior probabilities (prior over initial states)
%
% MDP.a(O,N)      - concentration parameters for A
% MDP.b{M}(N,N)   - concentration parameters for B
% MDP.c(N,N)      - concentration parameters for H
% MDP.d(N,1)      - concentration parameters for D
% MDP.e(P,1)      - concentration parameters for u
%
% optional:
% MDP.s(1,T)      - vector of true states  - for deterministic solutions
% MDP.o(1,T)      - vector of observations - for deterministic solutions
% MDP.u(1,T)      - vector of action       - for deterministic solutions
% MDP.w(1,T)      - vector of precisions   - for deterministic solutions
%
% MDP.alpha       - upper bound on precision (Gamma hyperprior – shape [1])
% MDP.beta        - precision over precision (Gamma hyperprior - rate  [1/2])
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
%
% MDP.un          - simulated neuronal encoding of hidden states
% MDP.xn          - simulated neuronal encoding of policies
% MDP.wn          - simulated neuronal encoding of precision (tonic)
% MDP.dn          - simulated dopamine responses (phasic)
% MDP.rt          - simulated reaction times
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
% $Id: spm_MDP_VB_X.m 6606 2015-11-22 18:57:07Z karl $


% deal with a sequence of trials
%==========================================================================

% options
%--------------------------------------------------------------------------
try, OPTIONS.plot;    catch, OPTIONS.plot    = 0; end
try, OPTIONS.gamma_u; catch, OPTIONS.gamma_u = 0; end

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
            try,  MDP(i).d = OUT(i - 1).d; end
        end
        
        % solve this trial
        %------------------------------------------------------------------
        OUT(i) = spm_MDP_VB_X(MDP(i),OPTS);
        
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
V   = MDP.V;                        % allowable policies (T - 1,Np)

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
Nf  = numel(MDP.B);                 % number of hidden state factors
Ng  = numel(MDP.A);                 % number of outcome factors
T   = size(MDP.V,1) + 1;            % number of transitions
Np  = size(MDP.V,2);                % number of allowable policies
for f = 1:Nf
    Nu(f) = size(MDP.B{f},3);       % number of hidden controls
    Ns(f) = size(MDP.B{f},1);       % number of hidden states
end
for g = 1:Ng
    No(g) = size(MDP.A{g},1);       % number of outcomes
end
p0  = exp(-8);                      % smallest probability
q0  = 1/16;                         % smallest probability


% parameters of generative model and policies
%==========================================================================

% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
for g = 1:Ng
    
    A{g} = spm_norm(MDP.A{g} + p0);
    
    % parameters (concentration parameters): A
    %----------------------------------------------------------------------
    if isfield(MDP,'a')
        qA{g} = spm_psi(MDP.a{g} + q0);
    else
        qA{g} = log(A{g});
    end
    
    % entropy
    %----------------------------------------------------------------------
    H{g} = spm_ent(qA{g});
    
end

% transition probabilities (priors)
%--------------------------------------------------------------------------
for f = 1:Nf
    for j = 1:Nu(f)
        
        % controlable transition probabilities
        %------------------------------------------------------------------
        B{f}(:,:,j)      = spm_norm(MDP.B{f}(:,:,j) + p0);
        
        % parameters (concentration parameters): B
        %------------------------------------------------------------------
        if isfield(MDP,'b')
            sB{f}(:,:,j) = spm_norm((MDP.b{f}(:,:,j) + q0) );
            rB{f}(:,:,j) = spm_norm((MDP.b{f}(:,:,j) + q0)');
        else
            sB{f}(:,:,j) = spm_norm(B{f}(:,:,j) );
            rB{f}(:,:,j) = spm_norm(B{f}(:,:,j)');
        end
        
    end
end


% priors over initial hidden states - concentration parameters
%--------------------------------------------------------------------------
for f = 1:Nf
    if isfield(MDP,'d')
        qD{f} = spm_psi(MDP.d{f} + q0);
    elseif isfield(MDP,'D')
        qD{f} = log(spm_norm(MDP.D{f} + q0));
    else
        qD{f} = spm_psi(ones(Ns(f),1));
    end
end


% prior preferences (log probabilities) : C
%--------------------------------------------------------------------------
for g = 1:Ng
    if isfield(MDP,'C')
        Vo{g} = MDP.C{g};
    else
        Vo{g} = zeros(No(g),1);
    end
    
    % assume constant preferences, if only final states are specified
    %--------------------------------------------------------------------------
    if size(Vo{g},2) ~= T
        Vo{g} = Vo{g}(:,end)*ones(1,T);
    end
    Vo{g}     = log(spm_softmax(Vo{g})); 
end

% precision defaults
%--------------------------------------------------------------------------
try, alpha = MDP.alpha;  catch, alpha = 16; end
try, beta  = MDP.beta;   catch, beta  = 1;  end

% initialise
%--------------------------------------------------------------------------
Ni    = 16;                         % number of VB iterations
rt    = zeros(1,T);                 % reaction times
wn    = zeros(T*Ni,1);              % simulated DA responses
for f = 1:Nf
    
    % initialise priors over states
    %----------------------------------------------------------------------
    try
        s(f,1) = MDP.s(f,1);
    catch
        s(f,1) = 1;
    end
    
    % initialise posteriors over states
    %----------------------------------------------------------------------
    xn{f} = zeros(Ni,Ns(f),T,T,Np);
    x{f}  = zeros(Ns(f),T,Np) + 1/Ns(f);
    X{f}  = zeros(Ns(f),T); 
    for k = 1:Np
        x{f}(:,1,k) = spm_softmax(qD{f});
    end
    
end

% initialise posteriors over polices and action
%--------------------------------------------------------------------------
P  = zeros([Nu (T - 1)]) - 16;
un = zeros(Np,T*Ni);
u  = zeros(Np,T - 1);
a  = zeros(f, T - 1);


% initial outcome (index)
%--------------------------------------------------------------------------
for g = 1:Ng
    try
        o(g,1) = MDP.o(g,1);
    catch
        ind    = num2cell(s(:,1));
        o(g,1) = find(rand < cumsum(A{g}(:,ind{:})),1);
    end
end

% expected rate parameter
%--------------------------------------------------------------------------
p     = 1:Np;                       % allowable policies
qbeta = beta;                       % initialise rate parameters
gu    = zeros(1,T)  + 1/qbeta;      % posterior precision (policy)

% solve
%==========================================================================
for t = 1:T
    
    % processing time and reset
    %----------------------------------------------------------------------
    tstart = tic;
    for f = 1:Nf
        x{f} = spm_softmax(log(x{f})/2);
    end
    
    % Variational updates (hidden states) under sequential policies
    %======================================================================
    for i = 1:Ni
        px    = x;
        F     = zeros(Np,T);
        for f = 1:Nf
            for k = p
                for j = 1:T
                    
                    % evaluate free energy and gradients (v = dFdx)
                    %======================================================
                    
                    % entropy term
                    %------------------------------------------------------
                    qx     = log(x{f}(:,j,k));
                    v      = qx;
                    
                    ind    = 1:Nf;
                    ind(f) = [];
                    xq     = cell(1,(Nf - 1));
                    for  q = 1:numel(ind)
                        xq{q} = x{ind(q)}(:,j,k);
                    end

                    % likelihood
                    %------------------------------------------------------
                    if j <= t
                        for g = 1:Ng
                            Aq  = spm_dot(qA{g},xq,ind + 1);
                            v   = v - Aq(o(g,j),:)';
                        end
                    end
                    
                    % emprical prior
                    %------------------------------------------------------
                    if j == 1, v = v - qD{f};                                         end
                    if j >  1, v = v - log(sB{f}(:,:,V(j - 1,k,f))*x{f}(:,j - 1,k));  end
                    if j <  T, v = v - log(rB{f}(:,:,V(j    ,k,f))*x{f}(:,j + 1,k));  end
                    
                    % free energy and belief updating
                    %------------------------------------------------------
                    F(k,j)       = F(k,j) + x{f}(:,j,k)'*v;
                    px{f}(:,j,k) = spm_softmax(qx - v/8);
                    
                    % record neuronal activity
                    %------------------------------------------------------
                    xn{f}(i,:,j,t,k) = x{f}(:,j,k);
                    
                end
            end
        end
        
        % hidden state updates
        %----------------------------------------------------------
        x = px;

    end
    
    % (negative path integral of) free energy of policies (Q)
    %======================================================================
    Q     = zeros(Np,T);
    for k = p
        for j = 1:T
            for g = 1:Ng
                xq    = cell(1,Nf);
                ind   = 1:Nf;
                for f = 1:Nf
                    xq{f}  = x{f}(:,j,k);
                end
                qx     = spm_dot(A{g},xq,ind + 1);
                Q(k,j) = Q(k,j) + qx'*(Vo{g}(:,j) - log(qx));
                Q(k,j) = Q(k,j) + spm_dot(H{g},xq,ind);
            end
        end
    end
    
    
    % variational updates - policies and precision
    %======================================================================
    F     = sum(F,2);
    Q     = sum(Q,2);
    p     = p(softmax(-F(p)) > 1/32);
    for i = 1:Ni
        
        % policy (u)
        %------------------------------------------------------------------
        qu = spm_softmax(gu(t)*Q(p) - F(p));
        pu = spm_softmax(gu(t)*Q(p));
        
        % precision (gu) with free energy gradients (v = -dF/dw)
        %------------------------------------------------------------------
        if OPTIONS.gamma_u
            gu(t) = 1/beta;
        else
            v     = qbeta - beta + (qu - pu)'*Q(p);
            qbeta = qbeta - v/2;
            gu(t) = 1/qbeta;
        end
        
        % simulated dopamine responses (precision at each iteration)
        %------------------------------------------------------------------
        n       = (t - 1)*Ni + i;
        u(p,t)  = qu;
        wn(n,1) = gu(t);
        un(p,n) = qu;
        
    end
    
    
    % Bayesian model averaging of hidden states over policies
    %----------------------------------------------------------------------
    for f = 1:Nf
        for i = 1:T
            X{f}(:,i) = squeeze(x{f}(:,i,:))*u(:,t);
        end
    end
    
    % processing time
    %----------------------------------------------------------------------
    rt(t) = toc(tstart);
    
    
    % action selection and sampling of next state (outcome)
    %======================================================================
    if t < T
        
        % posterior expectations about (remaining) actions
        %==================================================================
        for g  = 1:Ng
            
            % predicted outcome
            %--------------------------------------------------------------
            ind   = 1:Nf;
            for f = ind
                xq{f} = X{f}(:,t + 1);
            end
            v     = log(spm_dot(A{g},xq,ind + 1));
            
            % outcome under unique actions
            %--------------------------------------------------------------
            up    = unique(shiftdim(V(t,p,:),1),'rows');
            for i = 1:size(up,1)
                for f = ind
                    xq{f} = B{f}(:,:,up(i,f))*X{f}(:,t);
                end
                qo     = spm_dot(A{g},xq,ind + 1);
                dP     = (v - log(qo))'*qo;
                if     numel(Nu) < 2
                    P(up(i,1),t) = P(up(i,1),t) + dP;
                elseif numel(Nu) < 3
                    P(up(i,1),up(i,2),t) = P(up(i,1),up(i,2),t) + dP;
                elseif numel(Nu) < 4
                    P(up(i,1),up(i,2),up(i,3),t) = P(up(i,1),up(i,2),up(i,3),t) + dP;
                end
            end
        end
        
        % action selection
        %------------------------------------------------------------------
        if     numel(Nu) < 2
            P(:,t) = spm_softmax(alpha*spm_vec(P(:,t)));
        elseif numel(Nu) < 3
            P(:,:,t) = spm_softmax(alpha*spm_vec(P(:,:,t)));
        elseif numel(Nu) < 4
            P(:,:,:,t) = spm_softmax(alpha*spm_vec(P(:,:,:,t)));
        end

        % next action
        %------------------------------------------------------------------
        try
            a(:,t) = MDP.u(:,t);
        catch
            if     numel(Nu) < 2
                ind = find(rand < cumsum(spm_vec(P(:,t))),1);
                [i1] = ind2sub(Nu,ind);
                a(:,t) = [i1];
            elseif numel(Nu) < 3
                ind = find(rand < cumsum(spm_vec(P(:,:,t))),1);
                [i1,i2] = ind2sub(Nu,ind);
                a(:,t) = [i1,i2];
            elseif numel(Nu) < 4
                ind = find(rand < cumsum(spm_vec(P(:,:,:,t))),1);
                [i1,i2,i3] = ind2sub(Nu,ind);
                a(:,t) = [i1,i2,i3];
            end
        end
        
        % next sampled state
        %------------------------------------------------------------------
        try
            s(:,t + 1) = MDP.s(t + 1);
        catch
            for f = 1:Nf
                s(f,t + 1) = find(rand < cumsum(B{f}(:,s(t),a(t))),1);
            end
        end
        
        % next observed state
        %------------------------------------------------------------------
        try
            o(:,t + 1) = MDP.o(:,t + 1);
        catch
            for g = 1:Ng
                if     Nf < 2
                    o(g,t + 1) = find(rand < cumsum(A{g}(:,s(1,t + 1))),1);
                elseif Nf < 3
                    o(g,t + 1) = find(rand < cumsum(A{g}(:,s(1,t + 1),s(2,t + 1))),1);
                elseif Nf < 4
                    o(g,t + 1) = find(rand < cumsum(A{g}(:,s(1,t + 1),s(2,t + 1),s(3,t + 1))),1);
                end
            end
        end
        
        % save expected precision
        %------------------------------------------------------------------
        gu(1,t + 1)   = gu(t);
        
    end
    
end

% learning
%==========================================================================
for t = 1:T
    
    % mapping from hidden states to outcomes: a
    %----------------------------------------------------------------------
    if isfield(MDP,'a')
        for g = 1:Ng
            da     = sparse(o(g,t),1,1,No(g),1);
            for  f = 1:Nf
                da = spm_cross(da,X{f}(:,t));
            end
            MDP.a{g} = MDP.a{g} + da.*(MDP.a{g} > 0);
        end
    end
    
    % mapping from hidden states to hidden states: b(u)
    %----------------------------------------------------------------------
    if isfield(MDP,'b') && t > 1
        for f = 1:Nf
            for k = 1:Np
                v   = V(t - 1,k,f);
                db  = u(k,t - 1)*x{f}(:,t,k)*x{f}(:,t - 1,k)';
                MDP.b{f}(:,:,v) = MDP.b{f}(:,:,v) + db.*(MDP.b{f}(:,:,v) > 0);
            end
        end
    end
    
end

% initial hidden states: d
%--------------------------------------------------------------------------
if isfield(MDP,'d')
    for f = 1:Nf
        i = MDP.d{f} > 0;
        MDP.d{f}(i) = MDP.d{f}(i) + X{f}(i,1) - (MDP.d{f}(i) - 1)/16;
    end
end

% simulated dopamine (or cholinergic) responses
%--------------------------------------------------------------------------
dn    = 8*gradient(wn) + wn/8;

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
MDP.o   = o;              % outcomes at 1,...,T
MDP.s   = s;              % states at 1,...,T
MDP.u   = a;              % action at 1,...,T
MDP.w   = gu;             % posterior expectations of precision (policy)
MDP.v   = gx;             % posterior expectations of precision (states)
MDP.C   = Vo;             % utility

MDP.un  = un;             % simulated neuronal encoding of policies
MDP.xn  = Xn;             % simulated neuronal encoding of hidden states
MDP.wn  = wn;             % simulated neuronal encoding of precision
MDP.dn  = dn;             % simulated dopamine responses (deconvolved)
MDP.rt  = rt;             % simulated reaction time

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
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
             A(:,i,j,k) = A(:,i,j,k)/sum(A(:,i,j,k),1);
        end
    end
end

function A = spm_psi(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
             A(:,i,j,k) = psi(A(:,i,j,k)) - psi(sum(A(:,i,j,k)));
        end
    end
end

function H = spm_ent(qA)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(qA,2)
    for j = 1:size(qA,3)
        for k = 1:size(qA,4)
             H(i,j,k) = spm_softmax(qA(:,i,j,k))'*qA(:,i,j,k);
        end
    end
end
