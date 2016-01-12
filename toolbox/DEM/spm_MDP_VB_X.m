function [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
% active inference and learning using variational Bayes (factorised)
% FORMAT [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
%
% MDP.V(T - 1,P,F)      - P allowable policies of T moves over F factors
% or
% MDP.U(1,P,F)          - P allowable actions at each move
% MDP.T                 - number of outcomes
%
% MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes given hidden states
% MDP.B{F}(NF,NF,MF)    - transitions among hidden under MF control states
% MDP.C{G}(O,T)         - prior preferences over O outsomes in modality G
% MDP.D{F}(NF,1)        - prior probabilities over initial states
%
% MDP.a{G}              - concentration parameters for A
% MDP.b{F}              - concentration parameters for B
% MDP.d{F}              - concentration parameters for D
%
% optional:
% MDP.s(F,T)            - vector of true states - for each hidden factor
% MDP.o(G,T)            - vector of outcome     - for each outcome modality
% MDP.u(F,T - 1)        - vector of action      - for each hidden factor
% MDP.w(1,T)            - vector of precisions
%
% MDP.alpha             - precision – action selection [16]
% MDP.beta              - precision over precision (Gamma hyperprior - [1])
%
% OPTIONS.plot          - switch to suppress graphics: (default: [0])
%
% produces:
%
% MDP.P(M1,...,MF,T)    - probability of emitting action over time
% MDP.Q{F}(NF,T,P)      - expected hidden states under each policy
% MDP.X{F}(NF,T)        - and Bayesian model averages over policies
% MDP.R(P,T)            - conditional expectations over policies
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
% hidden states.
%
% This implementation equips agents with the prior beliefs that they will
% maximise expected free energy: expected free energy is the free energy
% of future outcomes under the posterior predictive distribution. This can
% be interpreted in several ways – most intuitively as minimising the KL
% divergence between predicted and preferred outcomes (specified as prior
% beliefs) – while simultaneously minimising ambiguity.
%
% This particular scheme is designed for any allowable policies or control
% sequences specified in MDP.V. Constraints on allowable policies can limit
% the numerics or combinatorics considerably.furthermore, the outcome space
% and hidden states can be defined in terms of factors; corresponding to
% sensory modalities and (functionally) segregated representations,
% respectively. This means, for each factor or subset of hidden states
% there are corresponding control states that determine the transition
% probabilities.
%
% This specification simplifies the generative model, allowing a fairly
% exhaustive model of potential outcomes. In brief, the agent encodes
% beliefs about hidden states in the past and in the future conditioned
% on each policy. The conditional expectations determine the (path
% integral) of free energy that then determines the prior over policies.
% This prior is used to create a predictive distribution over outcomes,
% which specifies the next action.
%
% In addition to state estimation and policy selection, the scheme also
% updates model parameters; including the state transition matrices,
% mapping to outcomes and the initial state. This is useful for learning
% the context.
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
% $Id: spm_MDP_VB_X.m 6672 2016-01-12 12:28:31Z karl $


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
try
    V = MDP.U;                        % allowable actions (1,Np)
    T = MDP.T;                        % number of transitions
catch
    V = MDP.V;                        % allowable policies (T - 1,Np)
    T = size(MDP.V,1) + 1;            % number of transitions
end

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
Nf  = numel(MDP.B);                 % number of hidden state factors
Ng  = numel(MDP.A);                 % number of outcome factors
Np  = size(V,2);                    % number of allowable policies
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
        B{f}(:,:,j) = spm_norm(MDP.B{f}(:,:,j) + p0);
        
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
        qD{f} = log(spm_norm(MDP.D{f} + p0));
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
try, alpha = MDP.alpha; catch, alpha = 16; end
try, beta  = MDP.beta;  catch, beta  = 1;  end

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
    xn{f} = zeros(Ni,Ns(f),T,T,Np) + 1/Ns(f);
    x{f}  = zeros(Ns(f),T,Np)      + 1/Ns(f);
    X{f}  = zeros(Ns(f),T);
    for k = 1:Np
        x{f}(:,1,k) = spm_softmax(qD{f});
    end
    
end

% initialise posteriors over polices and action
%--------------------------------------------------------------------------
P  = zeros([Nu,(T - 1)]);
un = zeros(Np,T*Ni);
u  = zeros(Np,T - 1);
a  = zeros(Nf,T - 1);


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
gu    = zeros(1,T) + 1/qbeta;       % posterior precision (policy)

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
    S     = size(V,1) + 1;
    for i = 1:Ni
        px    = x;
        F     = zeros(Np,S);
        for f = 1:Nf
            for k = p
                for j = 1:S
                    
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
                    if j <  S, v = v - log(rB{f}(:,:,V(j    ,k,f))*x{f}(:,j + 1,k));  end
                    
                    % free energy and belief updating
                    %------------------------------------------------------
                    F(k,j)       = F(k,j) - x{f}(:,j,k)'*v;
                    px{f}(:,j,k) = spm_softmax(qx - v/4);
                    
                    % record neuronal activity
                    %------------------------------------------------------
                    xn{f}(i,:,j,t,k) = x{f}(:,j,k);
                    
                end
            end
        end
        
        % hidden state updates and convergence
        %------------------------------------------------------------------
        x = px;
        if i > 1
            if all(sum(F - G,2) < 1/128)
                for g = i:Ni
                    for f = 1:Nf
                        for k = p
                            for j = 1:S
                                xn{f}(g,:,j,t,k) = x{f}(:,j,k);
                            end
                        end
                    end
                end
                break
            end
        end
        G = F;
    end
    
    % (negative path integral of) free energy of policies (Q)
    %======================================================================
    Q     = zeros(Np,S);
    for k = p
        for j = 1:S
            xq    = cell(1,Nf);
            ind   = 1:Nf;
            for f = 1:Nf
                xq{f}  = x{f}(:,j,k);
            end
            for g = 1:Ng
                qo     = spm_dot(A{g},xq,ind + 1);
                Q(k,j) = Q(k,j) + qo'*(Vo{g}(:,j) - log(qo));
                Q(k,j) = Q(k,j) + spm_dot(H{g},xq,ind);
            end
        end
    end
    
    
    % variational updates - policies and precision
    %======================================================================
    F     = sum(F,2);
    Q     = sum(Q,2);
    if ~isfield(MDP,'U')
        p = p((F(p) - max(F(p))) > -3);
    end
    for i = 1:Ni
        
        % policy (u)
        %------------------------------------------------------------------
        qu = spm_softmax(gu(t)*Q(p) + F(p));
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
    
    
    % Bayesian model averaging of hidden states (over policies)
    %----------------------------------------------------------------------
    for f = 1:Nf
        for i = 1:S
            X{f}(:,i) = reshape(x{f}(:,i,:),Ns(f),Np)*u(:,t);
        end
    end
    
    % processing time
    %----------------------------------------------------------------------
    rt(t) = toc(tstart);
    
    
    % action selection and sampling of next state (outcome)
    %======================================================================
    if t < T
        
        % posterior potential for (allowable) actions (for each modality)
        %==================================================================
        
        % unique combinations of actions
        %------------------------------------------------------------------
        up    = unique(shiftdim(V(t,p,:),1),'rows');
        
        % predicted hidden statesat the next time step
        %------------------------------------------------------------------
        ind   = 1:Nf;
        for f = 1:Nf
            xp{f} = X{f}(:,t + 1);
        end
        
        % predicted hidden states under each action
        %------------------------------------------------------------------
        Pu    = zeros(Nu);
        for i = 1:size(up,1)
            
            for f = 1:Nf
                xq{f} = B{f}(:,:,up(i,f))*X{f}(:,t);
            end
            
            % accumulate action potential over outcomes
            %--------------------------------------------------------------
            for g = 1:Ng
                
                % predicted outcome
                %----------------------------------------------------------
                po = spm_dot(A{g},xp,ind + 1);
                qo = spm_dot(A{g},xq,ind + 1);
                dP = (log(po) - log(qo))'*qo;
                
                % augment action potential
                %----------------------------------------------------------
                sub        = num2cell(up(i,:));
                Pu(sub{:}) = Pu(sub{:}) + dP + 16;
                
            end
        end
        
        % action selection - a softmax function of action potential
        %------------------------------------------------------------------
        sub         = repmat({':'},1,Nf);
        Pu(:)       = spm_softmax(alpha*Pu(:));
        P(sub{:},t) = Pu;
        
        % next action - sampled from beliefs about control states
        %------------------------------------------------------------------
        try
            a(:,t)  = MDP.u(:,t);
        catch
            ind     = find(rand < cumsum(Pu(:)),1);
            a(:,t)  = spm_ind2sub(Nu,ind);
        end
        
        % next sampled state - based on the current action
        %------------------------------------------------------------------
        try
            s(:,t + 1) = MDP.s(:,t + 1);
        catch
            for f = 1:Nf
                s(f,t + 1) = find(rand < cumsum(B{f}(:,s(f,t),a(f,t))),1);
            end
        end
        
        % next observed state
        %------------------------------------------------------------------
        try
            o(:,t + 1) = MDP.o(:,t + 1);
        catch
            for g = 1:Ng
                ind        = num2cell(s(:,t + 1));
                o(g,t + 1) = find(rand < cumsum(A{g}(:,ind{:})),1);
            end
        end
        
        % next expected precision
        %------------------------------------------------------------------
        gu(1,t + 1)   = gu(t);
        
        % update policy if necessary
        %------------------------------------------------------------------
        if isfield(MDP,'U') && t < (T - 1)
            for f = 1:Nf
                V(t,:,f) = a(f,t);
            end
            for j = 1:size(MDP.U,1)
                if (t + j) < T
                    V(t + j,:,:) = MDP.U(j,:,:);
                end
            end
            
            % and reinitialise expectations about hidden states
            %--------------------------------------------------------------
            j     = 1:t;
            for f = 1:Nf
                for k = 1:Np
                    x{f}(:,j,k) = x{f}(:,j,a(f,t));
                end
            end
            
        end
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
for f = 1:Nf
    Xn{f}    = zeros(Ni,Ns(f),T,T);
    for i = 1:T
        for k = 1:Np
            Xn{f}(:,:,:,i) = Xn{f}(:,:,:,i) + xn{f}(:,:,:,i,k)*u(k,i);
        end
    end
end


% assemble results and place in NDP structure
%--------------------------------------------------------------------------
MDP.P   = P;              % probability of action at time 1,...,T - 1
MDP.Q   = x;              % conditional expectations over N hidden states
MDP.X   = X;              % Bayesian model averages
MDP.R   = u;              % conditional expectations over policies
MDP.V   = V;              % policies
MDP.o   = o;              % outcomes at 1,...,T
MDP.s   = s;              % states at 1,...,T
MDP.u   = a;              % action at 1,...,T
MDP.w   = gu;             % posterior expectations of precision (policy)
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
            for l = 1:size(A,5)
                A(:,i,j,k,l) = A(:,i,j,k,l)/sum(A(:,i,j,k,l),1);
            end
        end
    end
end

function A = spm_psi(A)
% normalisation of a probability transition rate matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                A(:,i,j,k,l) = psi(A(:,i,j,k,l)) - psi(sum(A(:,i,j,k,l)));
            end
        end
    end
end

function H = spm_ent(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                H(i,j,k,l) = spm_softmax(A(:,i,j,k,l))'*A(:,i,j,k,l);
            end
        end
    end
end

function sub = spm_ind2sub(siz,ndx)
% subscripts from linear index
%--------------------------------------------------------------------------
n = numel(siz);
k = [1 cumprod(siz(1:end-1))];
for i = n:-1:1,
    vi       = rem(ndx - 1,k(i)) + 1;
    vj       = (ndx - vi)/k(i) + 1;
    sub(i,1) = vj;
    ndx      = vi;
end

