function [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
% active inference and learning using variational Bayes (factorised)
% FORMAT [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
%
% MDP.V(T - 1,P,F)      - P allowable policies (T - 1 moves) over F factors
% or
% MDP.U(1,P,F)          - P allowable actions at each move
% MDP.T                 - number of outcomes
%
% MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes given hidden states
% MDP.B{F}(NF,NF,MF)    - transitions among states under MF control states
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
% MDP.u(F,T - 1)        - vector of actions     - for each hidden factor
%
% MDP.alpha             - precision – action selection [16]
% MDP.beta              - precision over precision (Gamma hyperprior - [1])
% MDP.tau               - time constant for gradient descent
% MDP.eta               - learning rate for a and b parameters
%
% OPTIONS.plot          - switch to suppress graphics:  (default: [0])
% OPTIONS.gamma         - switch to suppress precision: (default: [0])
%
% produces:
%
% MDP.P(M1,...,MF,T)    - probability of emitting action M1,.. over time
% MDP.Q{F}(NF,T,P)      - expected hidden states under each policy
% MDP.X{F}(NF,T)        - and Bayesian model averages over policies
% MDP.R(P,T)            - conditional expectations over policies
%
% MDP.un          - simulated neuronal encoding of hidden states
% MDP.vn          - simulated neuronal prediction error
% MDP.xn          - simulated neuronal encoding of policies
% MDP.wn          - simulated neuronal encoding of precision (tonic)
% MDP.dn          - simulated dopamine responses (phasic)
% MDP.rt          - simulated reaction times
%
% MDP.F           - (Np x T) free energies over time
% MDP.G           - (Np x T) expected free energies over time
%
% This routine provides solutions of active inference (minimisation of
% variational free energy) using a generative model based upon a Markov
% decision process. The model and inference scheme is formulated
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
% the numerics or combinatorics considerably. Further, the outcome space
% and hidden states can be defined in terms of factors; corresponding to
% sensory modalities and (functionally) segregated representations,
% respectively. This means, for each factor or subset of hidden states
% there are corresponding control states that determine the transition
% probabilities.
%
% This specification simplifies the generative model, allowing a fairly
% exhaustive model of potential outcomes. In brief, the agent encodes
% beliefs about hidden states in the past (and in the future) conditioned
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
% $Id: spm_MDP_VB_X.m 6902 2016-10-08 19:28:25Z karl $


% deal with a sequence of trials
%==========================================================================

% options
%--------------------------------------------------------------------------
try, OPTIONS.plot;  catch, OPTIONS.plot  = 0; end
try, OPTIONS.gamma; catch, OPTIONS.gamma = 0; end

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
        
        % Bayesian model reduction
        %------------------------------------------------------------------
        if isfield(OPTIONS,'BMR')
            OUT(i) = spm_MDP_VB_sleep(OUT(i),OPTIONS.BMR);
        end
        
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
    T = MDP.T;                      % number of updates
    V = MDP.U;                      % allowable actions (1,Np)
    
catch
    V = MDP.V;                      % allowable policies (T - 1,Np)
    T = size(MDP.V,1) + 1;          % number of transitions
end

% eensure ppolicy length iis less than the number of updates
%--------------------------------------------------------------------------
if size(V,1) > (T - 1)
    V = V(1:(T - 1),:,:);
end

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
Ng  = numel(MDP.A);                 % number of outcome factors
Nf  = numel(MDP.B);                 % number of hidden state factors
Np  = size(V,2);                    % number of allowable policies
for f = 1:Nf
    Ns(f) = size(MDP.B{f},1);       % number of hidden states
    Nu(f) = size(MDP.B{f},3);       % number of hidden controls
end
for g = 1:Ng
    No(g) = size(MDP.A{g},1);       % number of outcomes
end

% parameters of generative model and policies
%==========================================================================

% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
p0    = exp(-16);
for g = 1:Ng
    
    MDP.A{g}  = spm_norm(MDP.A{g});
    
    % parameters (concentration parameters): A
    %----------------------------------------------------------------------
    if isfield(MDP,'a')
        A{g}  = spm_norm(MDP.a{g});
        qA{g} = spm_psi(MDP.a{g} + 1/16);
        wA{g} = 1./spm_cum(MDP.a{g}) - 1./(MDP.a{g} + p0);
        wA{g} = wA{g}.*(MDP.a{g} > 0);
    else
        A{g}  = MDP.A{g};
    end
    
end

% transition probabilities (priors)
%--------------------------------------------------------------------------
for f = 1:Nf
    for j = 1:Nu(f)
        
        % controlable transition probabilities
        %------------------------------------------------------------------
        MDP.B{f}(:,:,j) = spm_norm(MDP.B{f}(:,:,j));
        
        % parameters (concentration parameters): B
        %------------------------------------------------------------------
        if isfield(MDP,'b')
            sB{f}(:,:,j) = spm_norm(MDP.b{f}(:,:,j)  + p0);
            rB{f}(:,:,j) = spm_norm(MDP.b{f}(:,:,j)' + p0);
        else
            sB{f}(:,:,j) = spm_norm(MDP.B{f}(:,:,j)  + p0);
            rB{f}(:,:,j) = spm_norm(MDP.B{f}(:,:,j)' + p0);
        end
        
    end
end


% priors over initial hidden states - concentration parameters
%--------------------------------------------------------------------------
for f = 1:Nf
    if isfield(MDP,'d')
        D{f} = spm_norm(MDP.d{f});
    elseif isfield(MDP,'D')
        D{f} = spm_norm(MDP.D{f});
    else
        D{f} = spm_norm(ones(Ns(f),1));
        MDP.D{f} = D{f};
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
    %----------------------------------------------------------------------
    if size(Vo{g},2) == 1
        Vo{g} = repmat(Vo{g},1,T);
    end
    Vo{g}     = spm_log(spm_softmax(Vo{g}));
end

% precision defaults
%--------------------------------------------------------------------------
try, alpha = MDP.alpha; catch, alpha = 16;   end
try, beta  = MDP.beta;  catch, beta  = 1;    end
try, eta   = MDP.eta;   catch, eta   = 1;    end
try, tau   = MDP.tau;   catch, tau   = 4;    end
try, chi   = MDP.chi;   catch, chi   = 1/64; end

% initialise posteriors over states
%--------------------------------------------------------------------------
Ni    = 16;                         % number of VB iterations
for f = 1:Nf
    xn{f} = zeros(Ni,Ns(f),1,1,Np) + 1/Ns(f);
    vn{f} = zeros(Ni,Ns(f),1,1,Np);
    x{f}  = zeros(Ns(f),T,Np)      + 1/Ns(f);
    X{f}  = repmat(D{f},1,1);
    for k = 1:Np
        x{f}(:,1,k) = D{f};
    end
end

% initialise posteriors over polices and action
%--------------------------------------------------------------------------
P  = zeros([Nu,1]);
un = zeros(Np,1);
u  = zeros(Np,1);
a  = zeros(Nf,1);


% expected rate parameter
%--------------------------------------------------------------------------
p     = 1:Np;                       % allowable policies
qbeta = beta;                       % initialise rate parameters
gu    = 1/qbeta;                    % posterior precision (policy)

% solve
%==========================================================================
for t = 1:T
    
    % generate true states and outcomes
    %======================================================================
    
    % sampled state - based on previous action
    %----------------------------------------------------------------------
    for f = 1:Nf
        try
            s(f,t) = MDP.s(f,t);
        catch
            if t > 1
                ps = MDP.B{f}(:,s(f,t - 1),a(f,t - 1));
            else
                ps = spm_norm(MDP.D{f});
            end
            s(f,t) = find(rand < cumsum(ps),1);
        end
    end
    
    % sample outcome from true state if not specified
    %----------------------------------------------------------------------
    ind   = num2cell(s(:,t));
    for g = 1:Ng
        try
            o(g,t) = MDP.o(g,t);
        catch
            po     = MDP.A{g}(:,ind{:});
            o(g,t) = find(rand < cumsum(po),1);
        end
    end
    
    % posterior predictive density (prior for suborinate level)
    %------------------------------------------------------------------
    for f = 1:Nf
        if t > 1
            xq{f} = sB{f}(:,:,a(f,t - 1))*X{f}(:,t - 1);
        else
            xq{f} = X{f}(:,t);
        end
    end
    for g = 1:Ng
        if isfield(MDP,'link') || isfield(MDP,'demi')
            O{g,t} = spm_dot(A{g},xq);
        else
            O{g,t} = sparse(o(g,t),1,1,No(g),1);
        end
    end
    
    % generate outcomes from a subordinate MDP
    %======================================================================
    if isfield(MDP,'link')
        
        % use previous inversions (if available) to reproduce outcomes
        %------------------------------------------------------------------
        try
            mdp = MDP.mdp(t);
        catch
            mdp = MDP.MDP;
        end
        link       = MDP.link;
        mdp.factor = find(any(link,2));

        % priors over states (of subordinate level)
        %------------------------------------------------------------------
        for f = 1:size(link,1)
            i = find(link(f,:));
            if numel(i)
                
                % empirical priors
                %----------------------------------------------------------
                mdp.D{f} = O{i,t};
                
                % true state for lower level is the true outcome
                %----------------------------------------------------------
                try
                    mdp.s(f,1) = mdp.s(f,1);
                catch
                    mdp.s(f,1) = o(i,t);
                end
                
            else
                
                % otherwise use subordinate priors over states
                %----------------------------------------------------------
                try
                    mdp.s(f,1) = mdp.s(f,1);
                catch
                    if isfield(mdp,'D')
                        ps = spm_norm(mdp.D{f});
                    else
                        ps = spm_norm(ones(Ns(f),1));
                    end
                    mdp.s(f,1) = find(rand < cumsum(ps),1);
                end  
            end
            
        end
        
        % infer hidden states at the lower level (outcomes at this level)
        %------------------------------------------------------------------
        MDP.mdp(t) = spm_MDP_VB_X(mdp);

        % get inferred outcomes from subordinate MDP
        %==================================================================
        for g = 1:Ng
            i = find(link(:,g));
            if numel(i)
                O{g,t} = MDP.mdp(t).X{i}(:,1);
            end
        end
    end
    
    % generate outcomes
    %======================================================================
    if isfield(MDP,'demi')
        
        % use previous inversions (if available)
        %------------------------------------------------------------------
        try
            MDP.dem(t) = spm_ADEM_update(MDP.dem(t - 1));
        catch
            MDP.dem(t) = MDP.DEM;
        end
              
        % get inferred outcome (from Bayesian filtering)
        %------------------------------------------------------------------
        MDP.dem(t) = spm_MDP_DEM(MDP.dem(t),MDP.demi,O(:,t),o(:,t));
        for g = 1:Ng
            O{g,t} = MDP.dem(t).X{g}(:,end);
        end
        
    end
    
    
    % Variational updates
    %======================================================================
    
    % processing time and reset
    %----------------------------------------------------------------------
    tstart = tic;
    for f = 1:Nf
        x{f} = spm_softmax(spm_log(x{f})/4);
    end
    
    % Variational updates (hidden states) under sequential policies
    %======================================================================
    S     = size(V,1) + 1;
    F     = zeros(Np,1);
    for k = p
        dF    = 1;
        for i = 1:Ni
            F(k)  = 0;
            for j = 1:S
                
                % marginal likelihood over outcome factors
                %----------------------------------------------------------
                if j <= t
                    for f = 1:Nf
                        xq{f} = x{f}(:,j,k);
                    end
                    for g = 1:Ng
                        Ao{g} = spm_dot(A{g},[O(g,j) xq],(1:Nf) + 1);
                    end
                end
                
                for f = 1:Nf
                    
                    % hidden states for this time and policy
                    %------------------------------------------------------
                    sx = x{f}(:,j,k);
                    v  = spm_zeros(sx);
                    
                    % evaluate free energy and gradients (v = dFdx)
                    %------------------------------------------------------
                    if dF > 0
                        
                        % marginal likelihood over outcome factors
                        %--------------------------------------------------
                        if j <= t
                            for g = 1:Ng
                                Aq = spm_dot(Ao{g},xq,f);
                                v  = v + spm_log(Aq(:));
                            end
                        end
                        
                        % entropy
                        %--------------------------------------------------
                        qx  = spm_log(sx);
                        
                        % emprical priors
                        %--------------------------------------------------
                        if j < 2, v = v - qx + spm_log(D{f});                                    end
                        if j > 1, v = v - qx + spm_log(sB{f}(:,:,V(j - 1,k,f))*x{f}(:,j - 1,k)); end
                        if j < S, v = v - qx + spm_log(rB{f}(:,:,V(j    ,k,f))*x{f}(:,j + 1,k)); end
                        
                        % (negative) expected free energy
                        %--------------------------------------------------
                        F(k) = F(k) + sx'*v/Nf;
                        
                        % update
                        %--------------------------------------------------
                        sx   = spm_softmax(qx + v/tau);
                        
                    else
                        F(k) = G(k);
                    end
                    
                    % store update neuronal activity
                    %------------------------------------------------------
                    x{f}(:,j,k)      = sx;
                    xn{f}(i,:,j,t,k) = sx;
                    vn{f}(i,:,j,t,k) = v - mean(v);
                    
                end
            end
            
            % convergence
            %--------------------------------------------------------------
            if i > 1
                dF = F(k) - G(k);
            end
            G = F;
            
        end
    end
    
    % accumulate expected free energy of policies (Q)
    %======================================================================
    Q     = zeros(Np,1);
    for k = p
        for j = 1:S
            
            % get expected states for this policy and time point
            %--------------------------------------------------------------
            for f = 1:Nf
                xq{f} = x{f}(:,j,k);
            end
            
            % (negative) expected free energy
            %==============================================================
            
            % Bayesian surprise about states
            %--------------------------------------------------------------
            Q(k) = Q(k) + spm_MDP_G(A,xq);
            
            for g = 1:Ng
                
                % prior preferences about outcomes
                %----------------------------------------------------------
                qo   = spm_dot(A{g},xq);
                Q(k) = Q(k) + qo'*(Vo{g}(:,j));
                
                % Bayesian surprise about parameters
                %----------------------------------------------------------
                if isfield(MDP,'a')
                    Q(k) = Q(k) - spm_dot(wA{g},[qo xq]);
                end
                
            end
            
        end
    end
    
    % eliminate unlikely policies
    %----------------------------------------------------------------------
    if ~isfield(MDP,'U')
        p = p((F(p) - max(F(p))) > -3);
    else
        OPTIONS.gamma = 1;
    end
    
    % variational updates - policies and precision
    %======================================================================
    
    % previous expected precision
    %----------------------------------------------------------------------
    if t > 1
        gu(t) = gu(t - 1);
    end
    for i = 1:Ni
        
        % posterior and prior beliefs about policies
        %------------------------------------------------------------------
        qu = spm_softmax(gu(t)*Q(p) + F(p));
        pu = spm_softmax(gu(t)*Q(p));
        
        % precision (gu) with free energy gradients (v = -dF/dw)
        %------------------------------------------------------------------
        if OPTIONS.gamma
            gu(t) = 1/beta;
        else
            eg    = (qu - pu)'*Q(p);
            dFdg  = qbeta - beta + eg;
            qbeta = qbeta - dFdg/2;
            gu(t) = 1/qbeta;
        end
        
        % simulated dopamine responses (precision at each iteration)
        %------------------------------------------------------------------
        n       = (t - 1)*Ni + i;
        wn(n,1) = gu(t);
        un(p,n) = qu;
        u(p,t)  = qu;
        
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
    rt(t)    = toc(tstart);
    
    % record (negative) free energies
    %----------------------------------------------------------------------
    MDP.F(:,t) = F;
    MDP.G(:,t) = Q;
    
    
    % check for residual uncertainty in hierarchical schemes
    %----------------------------------------------------------------------
    if isfield(MDP,'factor')
        
        for f = MDP.factor
            qx   = X{f}(:,1);
            H(f) = qx'*spm_log(qx);
        end
        
        % break if there is no further uncertainty to resolve
        %------------------------------------------------------------------
        if sum(H) > - chi
            T = t;
        end
    end
    
    
    % action selection and sampling of next state (outcome)
    %======================================================================
    if t < T
        
        % marginal posterior probability of action (for each modality)
        %------------------------------------------------------------------
        Pu    = zeros([Nu,1]);
        for i = 1:Np
            sub        = num2cell(V(t,i,:));
            Pu(sub{:}) = Pu(sub{:}) + u(i,t);
        end
        
        % action selection - a softmax function of action potential
        %------------------------------------------------------------------
        sub         = repmat({':'},1,Nf);
        Pu(:)       = spm_softmax(alpha*log(Pu(:)));
        P(sub{:},t) = Pu;
        
        % next action - sampled from marginal posterior
        %------------------------------------------------------------------
        try
            a(:,t)  = MDP.u(:,t);
        catch
            ind     = find(rand < cumsum(Pu(:)),1);
            a(:,t)  = spm_ind2sub(Nu,ind);
        end
        
        
        % update policy and states for moving policies
        %------------------------------------------------------------------
        if isfield(MDP,'U')
            
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
            for f = 1:Nf
                for k = 1:Np
                    x{f}(:,:,k) = 1/Ns(f);
                end
            end
            
        end
    elseif t == T
        break;
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
            da       = da.*(MDP.a{g} > 0);
            MDP.a{g} = MDP.a{g} + da*eta;
            MDP.Fa   = spm_vec(da)'*spm_vec(qA{g}) - sum(spm_vec(spm_betaln(MDP.a{g})));
        end
    end
    
    % mapping from hidden states to hidden states: b(u)
    %----------------------------------------------------------------------
    if isfield(MDP,'b') && t > 1
        for f = 1:Nf
            for k = 1:Np
                v   = V(t - 1,k,f);
                db  = u(k,t - 1)*x{f}(:,t,k)*x{f}(:,t - 1,k)';
                db  = db.*(MDP.b{f}(:,:,v) > 0);
                MDP.b{f}(:,:,v) = MDP.b{f}(:,:,v) + db*eta;
            end
        end
    end
    
end

% initial hidden states:
%--------------------------------------------------------------------------
if isfield(MDP,'d')
    for f = 1:Nf
        i = MDP.d{f} > 0;
        MDP.d{f}(i) = MDP.d{f}(i) + X{f}(i,1);
    end
end

% simulated dopamine (or cholinergic) responses
%--------------------------------------------------------------------------
dn    = 8*gradient(wn) + wn/8;

% Bayesian model averaging of expected hidden states over policies
%--------------------------------------------------------------------------
for f = 1:Nf
    Xn{f} = zeros(Ni,Ns(f),T,T);
    Vn{f} = zeros(Ni,Ns(f),T,T);
    for i = 1:T
        for k = 1:Np
            Xn{f}(:,:,:,i) = Xn{f}(:,:,:,i) + xn{f}(:,:,1:T,i,k)*u(k,i);
            Vn{f}(:,:,:,i) = Vn{f}(:,:,:,i) + vn{f}(:,:,1:T,i,k)*u(k,i);
        end
    end
end

% use penultimate beliefs about moving policies
%--------------------------------------------------------------------------
if isfield(MDP,'U')
    u(:,T)  = [];
    un(:,(end - Ni + 1):end) = [];
end

% assemble results and place in NDP structure
%--------------------------------------------------------------------------
MDP.T   = T;              % number of belief updates
MDP.P   = P;              % probability of action at time 1,...,T - 1
MDP.Q   = x;              % conditional expectations over N hidden states
MDP.X   = X;              % Bayesian model averages over T outcomes
MDP.R   = u;              % conditional expectations over policies
MDP.V   = V;              % policies
MDP.o   = o;              % outcomes at 1,...,T
MDP.s   = s;              % states   at 1,...,T
MDP.u   = a;              % action   at 1,...,T - 1
MDP.w   = gu;             % posterior expectations of precision (policy)
MDP.C   = Vo;             % utility

MDP.un  = un;             % simulated neuronal encoding of policies
MDP.vn  = Vn;             % simulated neuronal prediction error
MDP.xn  = Xn;             % simulated neuronal encoding of hidden states
MDP.wn  = wn;             % simulated neuronal encoding of precision
MDP.dn  = dn;             % simulated dopamine responses (deconvolved)
MDP.rt  = rt;             % simulated reaction time (seconds)


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

function A = spm_log(A)
% log of numeric array plus a small constant
%--------------------------------------------------------------------------
A  = log(A + 1e-16);


function A = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                S = sum(A(:,i,j,k,l),1);
                if S > 0
                    A(:,i,j,k,l) = A(:,i,j,k,l)/S;
                else
                    A(:,i,j,k,l) = 1/size(A,1);
                end
            end
        end
    end
end

function A = spm_cum(A)
% summation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                A(:,i,j,k,l) = sum(A(:,i,j,k,l),1);
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

function C = spm_joint(A,B)
% subscripts from linear index
%--------------------------------------------------------------------------
C = zeros(size(A,1),size(A,2),size(B,2));
for i = 1:size(A,1)
    C(i,:,:) = spm_cross(A(i,:),B(:,i));
end
C = spm_norm(C);


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



