function [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
% active inference and learning using variational message passing
% FORMAT [MDP] = spm_MDP_VB_X(MDP,OPTIONS)
%
% Input; MDP(m,n)       - structure array of m models over n epochs
%
% MDP.V(T - 1,P,F)      - P allowable policies (T - 1 moves) over F factors
% or
% MDP.U(1,P,F)          - P allowable actions at each move
% MDP.T                 - number of outcomes
%
% MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes given hidden states
% MDP.B{F}(NF,NF,MF)    - transitions among states under MF control states
% MDP.C{G}(O,T)         - (log) prior preferences for outcomes (modality G)
% MDP.D{F}(NF,1)        - prior probabilities over initial states
% MDP.E(P,1)            - prior probabilities over policies
%
% MDP.a{G}              - concentration parameters for A
% MDP.b{F}              - concentration parameters for B
% MDP.c{G}              - concentration parameters for C
% MDP.d{F}              - concentration parameters for D
% MDP.e                 - concentration parameters for E
%
% optional:
% MDP.s(F,T)            - matrix of true states - for each hidden factor
% MDP.o(G,T)            - matrix of outcomes    - for each outcome modality
% or .O{G}(O,T)         - likelihood matrix     - for each outcome modality
% MDP.u(F,T - 1)        - vector of actions     - for each hidden factor
%
% MDP.alpha             - precision - action selection [512]
% MDP.beta              - precision over precision (Gamma hyperprior - [1])
% MDP.tau               - time constant for gradient descent [4]
% MDP.eta               - learning rate for model parameters
% MDP.zeta              - Occam's window for polcies [3]
% MDP.erp               - resetting of initial states, to simulate ERPs [4]
%
% MDP.demi.C            - Mixed model: cell array of true causes (DEM.C)
% MDP.demi.U            - Bayesian model average (DEM.U) see: spm_MDP_DEM
% MDP.link              - link array to generate outcomes from
%                         subordinate MDP; for deep (hierarchical) models
%
% OPTIONS.plot          - switch to suppress graphics:  (default: [0])
% OPTIONS.gamma         - switch to suppress precision: (default: [0])
% OPTIONS.D             - switch to update initial states over epochs
% OPTIONS.BMR           - Bayesian model reduction for multiple trials
%                         see: spm_MDP_VB_sleep(MDP,BMR)
% Outputs:
%
% MDP.P(M1,...,MF,T)    - probability of emitting action M1,.. over time
% MDP.Q{F}(NF,T,P)      - expected hidden states under each policy
% MDP.X{F}(NF,T)        - and Bayesian model averages over policies
% MDP.R(P,T)            - response: conditional expectations over policies
%
% MDP.un          - simulated neuronal encoding of hidden states
% MDP.vn          - simulated neuronal prediction error
% MDP.xn          - simulated neuronal encoding of policies
% MDP.wn          - simulated neuronal encoding of precision (tonic)
% MDP.dn          - simulated dopamine responses (phasic)
% MDP.rt          - simulated reaction times
%
% MDP.F           - (P x T) (negative) free energies over time
% MDP.G           - (P x T) (negative) expected free energies over time
% MDP.H           - (1 x T) (negative) total free energy over time
% MDP.Fa          - (1 x 1) (negative) free energy of parameters (a)
% MDP.Fb          - ...
%
% This routine provides solutions of active inference (minimisation of
% variational free energy) using a generative model based upon a Markov
% decision process (or hidden Markov model, in the absence of action). The
% model and inference scheme is formulated in discrete space and time. This
% means that the generative model (and process) are  finite state machines
% or hidden Markov models whose dynamics are given by transition
% probabilities among states and the likelihood corresponds to a particular
% outcome conditioned upon hidden states.
%
% When supplied with outcomes, in terms of their likelihood (O) in the
% absence of any policy specification, this scheme will use variational
% message passing to optimise expectations about latent or hidden states
% (and likelihood (A) and prior (B) probabilities). In other words, it will
% invert a hidden Markov model. When  called with policies, it will
% generate outcomes that are used to infer optimal policies for active
% inference.
%
% This implementation equips agents with the prior beliefs that they will
% maximise expected free energy: expected free energy is the free energy of
% future outcomes under the posterior predictive distribution. This can be
% interpreted in several ways – most intuitively as minimising the KL
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
% the context. Likelihood and prior probabilities can be specified in terms
% of concentration parameters (of a Dirichlet distribution (a,b,c,..). If
% the corresponding (A,B,C,..) are supplied, they will be used to generate
% outcomes; unless called without policies (in hidden Markov model mode).
% In this case, the (A,B,C,..) are treated as posterior estimates.
%
% If supplied with a structure array, this routine will automatically step
% through the implicit sequence of epochs (implicit in the number of
% columns of the array). If the array has multiple rows, each row will be
% treated as a separate model or agent. This enables agents to communicate
% through acting upon a common set of hidden factors, or indeed sharing the
% same outcomes.
%
% See also: spm_MDP, which uses multiple future states and a mean field
% approximation for control states – but allows for different actions at
% all times (as in control problems).
%
% See also: spm_MDP_game_KL, which uses a very similar formulation but just
% maximises the KL divergence between the posterior predictive distribution
% over hidden states and those specified by preferences or prior beliefs.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_X.m 7398 2018-08-15 14:50:27Z thomas $


% deal with a sequence of trials
%==========================================================================

% options
%--------------------------------------------------------------------------
try, OPTIONS.plot;  catch, OPTIONS.plot  = 0; end
try, OPTIONS.gamma; catch, OPTIONS.gamma = 0; end
try, OPTIONS.D;     catch, OPTIONS.D     = 0; end

% check MDP specification
%--------------------------------------------------------------------------
MDP = spm_MDP_check(MDP);

% if there are multiple trials ensure that parameters are updated
%--------------------------------------------------------------------------
if size(MDP,2) > 1
    
    % plotting options
    %----------------------------------------------------------------------
    GRAPH        = OPTIONS.plot;
    OPTIONS.plot = 0;
    
    for i = 1:size(MDP,2)
        for m = 1:size(MDP,1)
            
            % update concentration parameters
            %--------------------------------------------------------------
            if i > 1
                try,  MDP(m,i).a = OUT(m,i - 1).a; end
                try,  MDP(m,i).b = OUT(m,i - 1).b; end
                try,  MDP(m,i).d = OUT(m,i - 1).d; end
                try,  MDP(m,i).c = OUT(m,i - 1).c; end
                try,  MDP(m,i).e = OUT(m,i - 1).e; end
                
                % update initial states (post-diction)
                %----------------------------------------------------------
                if OPTIONS.D
                    for f = 1:numel(MDP(m,i).D)
                        MDP(m,i).D{f} = OUT(m,i - 1).X{f}(:,1);
                    end
                end
            end
        end
        
        % solve this trial (for all models synchronously)
        %------------------------------------------------------------------
        OUT(:,i) = spm_MDP_VB_X(MDP(:,i),OPTIONS);
        
        % Bayesian model reduction
        %------------------------------------------------------------------
        if isfield(OPTIONS,'BMR')
            for m = 1:size(MDP,1)
                OUT(m,i) = spm_MDP_VB_sleep(OUT(m,i),OPTIONS.BMR);
            end
        end
        
    end
    
    % plot summary statistics - over trials
    %----------------------------------------------------------------------
    MDP = OUT;
    if GRAPH
        if ishandle(GRAPH)
            figure(GRAPH); clf
        else
            spm_figure('GetWin','MDP'); clf
        end
        spm_MDP_VB_game(MDP(1,:))
    end
    return
end


% set up and preliminaries
%==========================================================================

% defaults
%--------------------------------------------------------------------------
try, alpha = MDP(1).alpha; catch, alpha = 512;  end
try, beta  = MDP(1).beta;  catch, beta  = 1;    end
try, zeta  = MDP(1).zeta;  catch, zeta  = 3;    end
try, eta   = MDP(1).eta;   catch, eta   = 1;    end
try, tau   = MDP(1).tau;   catch, tau   = 4;    end
try, chi   = MDP(1).chi;   catch, chi   = 1/64; end
try, erp   = MDP(1).erp;   catch, erp   = 4;    end

% preclude precision updates for moving policies
%--------------------------------------------------------------------------
if isfield(MDP,'U'), OPTIONS.gamma = 1;         end


for m = 1:size(MDP,1)
    
    if isfield(MDP(m),'O') && size(MDP(m).U,2) < 2
        
        % no policies – assume hidden Markov model (HMM)
        %------------------------------------------------------------------
        T(m) = size(MDP(m).O{1},2);         % HMM mode
        V{m} = ones(T - 1,1);               % single 'policy'
        HMM  = 1;
        
    elseif isfield(MDP(m),'U')
        
        % called with repeatable actions (U,T)
        %------------------------------------------------------------------
        T(m) = MDP(m).T;                    % number of updates
        V{m} = MDP(m).U;                    % allowable actions (1,Np)
        HMM  = 0;
        
    elseif isfield(MDP(m),'V')
        
        % full sequential policies (V)
        %------------------------------------------------------------------
        V{m} = MDP(m).V;                    % allowable policies (T - 1,Np)
        T(m) = size(MDP(m).V,1) + 1;        % number of transitions
        HMM  = 0;
        
    else
        sprintf('Please specify MDP(%d).U, MDP(%i).V or MDP(%d).O',m), return
    end
    
end

% initialise model-specific variables
%--------------------------------------------------------------------------
T     = T(1);                              % number of time steps
Ni    = 16;                                % number of VB iterations
for m = 1:size(MDP,1)
    
    % ensure policy length is less than the number of updates
    %----------------------------------------------------------------------
    if size(V{m},1) > (T - 1)
        V{m} = V{m}(1:(T - 1),:,:);
    end
    
    % numbers of transitions, policies and states
    %----------------------------------------------------------------------
    Ng(m) = numel(MDP(m).A);               % number of outcome factors
    Nf(m) = numel(MDP(m).B);               % number of hidden state factors
    Np(m) = size(V{m},2);                  % number of allowable policies
    for f = 1:Nf(m)
        Ns(m,f) = size(MDP(m).B{f},1);     % number of hidden states
        Nu(m,f) = size(MDP(m).B{f},3);     % number of hidden controls
    end
    for g = 1:Ng(m)
        No(m,g) = size(MDP(m).A{g},1);     % number of outcomes
    end
    
    % parameters of generative model and policies
    %======================================================================
    
    % likelihood model (for a partially observed MDP)
    %----------------------------------------------------------------------
    for g = 1:Ng(m)
        
        % ensure probabilities are normalised  : A
        %------------------------------------------------------------------
        MDP(m).A{g} = spm_norm(MDP(m).A{g});
        
        % parameters (concentration parameters): a
        %------------------------------------------------------------------
        if isfield(MDP,'a')
            A{m,g}  = spm_norm(MDP(m).a{g});
        else
            A{m,g}  = spm_norm(MDP(m).A{g});
        end
        
        % prior concentration paramters for complexity (and novelty)
        %------------------------------------------------------------------
        if isfield(MDP,'a')
            pA{m,g} = MDP(m).a{g};
            wA{m,g} = spm_wnorm(MDP(m).a{g}).*(pA{m,g} > 0);
        end
        
    end
    
    % transition probabilities (priors)
    %----------------------------------------------------------------------
    for f = 1:Nf(m)
        for j = 1:Nu(m,f)
            
            % controlable transition probabilities : B
            %--------------------------------------------------------------
            MDP(m).B{f}(:,:,j) = spm_norm(MDP(m).B{f}(:,:,j));
            
            % parameters (concentration parameters): b
            %--------------------------------------------------------------
            if isfield(MDP,'b') && ~HMM
                sB{m,f}(:,:,j) = spm_norm(MDP(m).b{f}(:,:,j) );
                rB{m,f}(:,:,j) = spm_norm(MDP(m).b{f}(:,:,j)');
            else
                sB{m,f}(:,:,j) = spm_norm(MDP(m).B{f}(:,:,j) );
                rB{m,f}(:,:,j) = spm_norm(MDP(m).B{f}(:,:,j)');
            end
            
        end
        
        % prior concentration paramters for complexity
        %------------------------------------------------------------------
        if isfield(MDP,'b')
            pB{m,f} = MDP(m).b{f};
        end
        
    end
    
    % priors over initial hidden states - concentration parameters
    %----------------------------------------------------------------------
    for f = 1:Nf(m)
        if isfield(MDP,'d')
            D{m,f} = spm_norm(MDP(m).d{f});
        elseif isfield(MDP,'D')
            D{m,f} = spm_norm(MDP(m).D{f});
        else
            D{m,f} = spm_norm(ones(Ns(m,f),1));
            MDP(m).D{f} = D{m,f};
        end
        
        % prior concentration paramters for complexity
        %------------------------------------------------------------------
        if isfield(MDP,'d')
            pD{m,f} = MDP(m).d{f};
            wD{m,f} = spm_wnorm(MDP(m).d{f});
        end
    end
    
    % priors over policies - concentration parameters
    %----------------------------------------------------------------------
    if isfield(MDP,'e')
        E{m} = spm_norm(MDP(m).e);
    elseif isfield(MDP,'E')
        E{m} = spm_norm(MDP(m).E);
    else
        E{m} = spm_norm(ones(Np(m),1));
    end
    qE{m}    = spm_log(E{m});
    
    % prior concentration paramters for complexity
    %----------------------------------------------------------------------
    if isfield(MDP,'e')
        pE{m} = MDP(m).e;
    end
    
    % prior preferences (log probabilities) : C
    %----------------------------------------------------------------------
    for g = 1:Ng(m)
        if isfield(MDP,'c')
            C{m,g}  = spm_psi(MDP(m).c{g} + 1/32);
            pC{m,g} = MDP(m).c{g};
        elseif isfield(MDP,'C')
            C{m,g}  = MDP(m).C{g};
        else
            C{m,g}  = zeros(No(m,g),1);
        end
        
        % assume time-invariant preferences, if unspecified
        %------------------------------------------------------------------
        if size(C{m,g},2) == 1
            C{m,g} = repmat(C{m,g},1,T);
            if isfield(MDP,'c')
                MDP(m).c{g} = repmat(MDP(m).c{g},1,T);
                pC{m,g}     = repmat(pC{m,g},1,T);
            end
        end
        C{m,g} = spm_log(spm_softmax(C{m,g}));
    end
    
    % initialise  posterior expectations of hidden states
    %----------------------------------------------------------------------
    for f = 1:Nf(m)
        xn{m,f} = zeros(Ni,Ns(m,f),1,1,Np(m)) + 1/Ns(m,f);
        vn{m,f} = zeros(Ni,Ns(m,f),1,1,Np(m));
        x{m,f}  = zeros(Ns(m,f),T,Np(m))      + 1/Ns(m,f);
        X{m,f}  = repmat(D{m,f},1,1);
        for k = 1:Np(m)
            x{m,f}(:,1,k) = D{m,f};
        end
    end
    
    % initialise posteriors over polices and action
    %----------------------------------------------------------------------
    P{m}  = zeros([Nu(m,:),1]);
    un{m} = zeros(Np(m),1);
    u{m}  = zeros(Np(m),1);
    
    % if there is only one policy
    %----------------------------------------------------------------------
    if Np(m) == 1
        u{m} = ones(Np(m),T);
    end
    
    % if states have not been specified set to 0
    %----------------------------------------------------------------------
    s{m}  = zeros(Nf(m),T);
    try
        i = find(MDP(m).s);
        s{m}(i) = MDP(m).s(i);
    end
    MDP(m).s = s{m};
    
    % if outcomes have not been specified set to 0
    %----------------------------------------------------------------------
    o{m}  = zeros(Ng(m),T);
    try
        i = find(MDP(m).o);
        o{m}(i) = MDP(m).o(i);
    end
    MDP(m).o = o{m};
    
    % (indices of) plausible (allowable) policies
    %----------------------------------------------------------------------
    p{m}  = 1:Np(m);
    
    % expected rate parameter (precision of posterior over policies)
    %----------------------------------------------------------------------
    qb{m} = beta;                          % initialise rate parameters
    w{m}  = 1/qb{m};                       % posterior precision (policy)
    
end

% belief updating for shared states and outcomes (multiple models)
%==========================================================================

% ensure any outcome generating agent is updated first
%--------------------------------------------------------------------------
for m = 1:size(MDP,1)
    n      = -MDP(m).o;
    N(m,:) = mode(n.*(n > 0),1);
end
n     = mode(N,1);
for t = 1:T
    if n(t) > 0
        M(t,:) = circshift((1:size(MDP,1)),[0 (1 - n(t))]);
    else
        M(t,:) = 1;
    end
end


% belief updating over successive time points
%==========================================================================
for t = 1:T
    
    % generate hidden states and outcomes for each agent or model
    %======================================================================
    for m = M(t,:)
        
        if ~HMM % not required for HMM
            
            % sample state, if not specified
            %--------------------------------------------------------------
            for f = 1:Nf(m)
                
                % the next state is generated by action
                %----------------------------------------------------------
                if MDP(m).s(f,t) == 0
                    if t > 1
                        ps = MDP(m).B{f}(:,MDP(m).s(f,t - 1),MDP(m).u(f,t - 1));
                    else
                        ps = spm_norm(MDP(m).D{f});
                    end
                    MDP(m).s(f,t) = find(rand < cumsum(ps),1);
                end
                
            end
            
            % posterior predictive density
            %--------------------------------------------------------------
            for f = 1:Nf(m)
                
                % under selected action (xqq)
                %----------------------------------------------------------
                if t > 1
                    xqq{m,f} = sB{m,f}(:,:,MDP(m).u(f,t - 1))*X{m,f}(:,t - 1);
                else
                    xqq{m,f} = X{m,f}(:,t);
                end
                
                % Bayesian model average (xq)
                %----------------------------------------------------------
                xq{m,f} = X{m,f}(:,t);
                
            end
            
            % sample outcome, if not specified
            %--------------------------------------------------------------
            for g = 1:Ng(m)
                
                if MDP(m).o(g,t) < 0
                    
                    % outcome is generated by model n
                    %------------------------------------------------------
                    n = -MDP(m).o(g,t);
                    MDP(m).n(g,t) = n;
                    if n == m
                        
                        % outcome that minimises expected free energy
                        %--------------------------------------------------
                        po    = spm_dot(A{m,g},xqq(m,:));
                        px    = spm_vec(spm_cross(xqq(m,:)));
                        F     = zeros(No(m,g),1);
                        for i = 1:No(m,g)
                            xp   = MDP(m).A{g}(i,:);
                            xp   = spm_norm(spm_vec(xp));
                            F(i) = spm_vec(px)'*spm_log(xp) + spm_log(po(i));
                        end
                        po            = spm_softmax(F*512);
                        MDP(m).o(g,t) = find(rand < cumsum(po),1);
                        
                    else
                        
                        % outcome from model n
                        %--------------------------------------------------
                        MDP(m).o(g,t) = MDP(n).o(g,t);
                        
                    end
                    
                elseif MDP(m).o(g,t) == 0
                    
                    % sample outcome from the generative process
                    %------------------------------------------------------
                    ind           = num2cell(MDP(m).s(:,t));
                    po            = MDP(m).A{g}(:,ind{:});
                    MDP(m).o(g,t) = find(rand < cumsum(po),1);
                    
                end
            end
            
        end % HMM
        
        % get probabilistic outcomes from samples or subordinate level
        %==================================================================
        
        % get outcome likelihood (O)
        %------------------------------------------------------------------
        for g = 1:Ng(m)
            
            % specified as a likelihood or observation (HMM)
            %--------------------------------------------------------------
            if HMM
                O{m}{g,t} = MDP(m).O{g}(:,t);
            else
                O{m}{g,t} = sparse(MDP(m).o(g,t),1,1,No(m,g),1);
            end
        end
        
        % generate outcomes from a subordinate MDP
        %==================================================================
        if isfield(MDP,'link')

            % use previous inversions (if available) to reproduce outcomes
            %--------------------------------------------------------------
            try
                mdp = MDP(m).mdp(t);
            catch
                mdp = MDP(m).MDP;
            end

            % priors over states (of subordinate level)
            %--------------------------------------------------------------
            mdp.factor = [];
            for f = 1:size(MDP(m).link,1)
                for g = 1:size(MDP(m).link,2)
                    if ~isempty(MDP(m).link{f,g})
                        
                        % subiordinate state has hierarchical constraints
                        %--------------------------------------------------
                        mdp.factor(end + 1) = f;
                                                
                        % empirical priors
                        %--------------------------------------------------
                        O{m}{g,t} = spm_dot(A{m,g},xq(m,:));
                        mdp.D{f}  = MDP(m).link{f,g}*O{m}{g,t};
                        
                        % outcomes (i.e., states) are generated by model n
                        %--------------------------------------------------
                        if isfield(MDP(m),'n')
                            n    = MDP(m).n(g,t);
                            if m == n
                                ps         = MDP(m).link{f,g}(:,MDP(m).o(g,t));
                                mdp.s(f,1) = find(ps);
                            else
                                mdp.s(f,1) = MDP(n).mdp(t).s(f,1);
                            end
                        end
                        
                        % hidden state for lower level is the outcome
                        %--------------------------------------------------
                        try
                            mdp.s(f,1) = mdp.s(f,1);
                        catch
                            ps         = MDP(m).link{f,g}(:,MDP(m).o(g,t));
                            mdp.s(f,1) = find(ps);
                        end

                    end
                end
            end
            
            % infer hidden states at lower level (outcomes at this level)
            %==============================================================
            MDP(m).mdp(t) = spm_MDP_VB_X(mdp);
            
            % get inferred outcomes from subordinate MDP
            %--------------------------------------------------------------
            for f = 1:size(MDP(m).link,1)
                for g = 1:size(MDP(m).link,2)
                    if ~isempty(MDP(m).link{f,g})
                        O{m}{g,t} = MDP(m).link{f,g}'*MDP(m).mdp(t).X{f}(:,1);
                    end
                end
            end
            
        end % end of hierarchical mode
        
        
        % generate outcomes from a generalised Bayesian filter
        %==================================================================
        if isfield(MDP,'demi')
            
            % use previous inversions (if available)
            %--------------------------------------------------------------
            try
                MDP(m).dem(t) = spm_ADEM_update(MDP(m).dem(t - 1));
            catch
                MDP(m).dem(t) = MDP(m).DEM;
            end
            
            % get prior over outcomes
            %--------------------------------------------------------------
            for g = 1:Ng(m)
                O{m}{g,t} = spm_dot(A{m,g},xqq(m,:));
            end
            
            % get posterior outcome from Bayesian filtering
            %--------------------------------------------------------------
            MDP(m).dem(t) = spm_MDP_DEM(MDP(m).dem(t),MDP(m).demi,O{m}(:,t),MDP(m).o(:,t));
            for g = 1:Ng(m)
                O{m}{g,t} = MDP(m).dem(t).X{g}(:,end);
            end
        end
        
        % likelihood (for multiple modalities)
        %==================================================================
        L{m,t} = 1;
        for g = 1:Ng(m)
            L{m,t} = L{m,t}.*spm_dot(A{m,g},O{m}{g,t});
        end
        
        
        % Variational updates (skip to t = T in HMM mode)
        %==================================================================
        if ~HMM || T == t
            
            % eliminate unlikely policies
            %--------------------------------------------------------------
            if ~isfield(MDP,'U') && t > 1
                F    = log(u{m}(p{m},t - 1));
                p{m} = p{m}((F - max(F)) > -zeta);
            end
            
            % processing time and reset
            %--------------------------------------------------------------
            tstart = tic;
            for f = 1:Nf(m)
                x{m,f} = spm_softmax(spm_log(x{m,f})/erp);
            end
            
            % Variational updates (hidden states) under sequential policies
            %==============================================================
            
            
            % variational message passing (VMP)
            %--------------------------------------------------------------
            S     = size(V{m},1) + 1;   % horizon
            if isfield(MDP,'U')
                R = t;
            else
                R = S;
            end
            F     = zeros(Np(m),1);
            for k = p{m}                % loop over plausible policies
                dF    = 1;              % reset criterion for this policy
                for i = 1:Ni            % iterate belief updates
                    F(k)  = 0;          % reset free energy for this policy
                    for j = 1:S         % loop over future time points
                        
                        % curent posterior over outcome factors
                        %--------------------------------------------------
                        if j <= t
                            for f = 1:Nf(m)
                                xq{m,f} = full(x{m,f}(:,j,k));
                            end
                        end
                        
                        for f = 1:Nf(m)
                            
                            % hidden states for this time and policy
                            %----------------------------------------------
                            sx = full(x{m,f}(:,j,k));
                            qL = zeros(Ns(m,f),1);
                            v  = zeros(Ns(m,f),1);
                            
                            % evaluate free energy and gradients (v = dFdx)
                            %----------------------------------------------
                            if dF > exp(-8) || i > 4
                                
                                % marginal likelihood over outcome factors
                                %------------------------------------------
                                if j <= t
                                    qL = spm_dot(L{m,j},xq(m,:),f);
                                    qL = spm_log(qL(:));
                                end
                                
                                % entropy
                                %------------------------------------------
                                qx  = spm_log(sx);
                                
                                % emprical priors (forward messages)
                                %------------------------------------------
                                if j < 2
                                    px = spm_log(D{m,f});
                                    v  = v + px + qL - qx;
                                else
                                    px = spm_log(sB{m,f}(:,:,V{m}(j - 1,k,f))*x{m,f}(:,j - 1,k));
                                    v  = v + px + qL - qx;
                                end
                                
                                % emprical priors (backward messages)
                                %------------------------------------------
                                if j < R
                                    px = spm_log(rB{m,f}(:,:,V{m}(j    ,k,f))*x{m,f}(:,j + 1,k));
                                    v  = v + px + qL - qx;
                                end
                                
                                % (negative) expected free energy
                                %------------------------------------------
                                F(k) = F(k) + sx'*v/Nf(m);
                                
                                % update
                                %------------------------------------------
                                v    = v - mean(v);
                                sx   = spm_softmax(qx + v/tau);
                                
                            else
                                F(k) = G(k);
                            end
                            
                            % store update neuronal activity
                            %----------------------------------------------
                            x{m,f}(:,j,k)      = sx;
                            xq{m,f}            = sx;
                            xn{m,f}(i,:,j,t,k) = sx;
                            vn{m,f}(i,:,j,t,k) = v;
                            
                        end
                    end
                    
                    % convergence
                    %------------------------------------------------------
                    if i > 1
                        dF = F(k) - G(k);
                    end
                    G = F;
                    
                end
            end
            
            % accumulate expected free energy of policies (Q)
            %==============================================================
            pu  = 1;                               % empirical prior
            qu  = 1;                               % posterior
            Q   = zeros(Np(m),1);                  % expected free energy
            if Np(m) > 1
                for k = p{m}
                    
                    % Bayesian surprise about inital conditions
                    %------------------------------------------------------
                    if isfield(MDP,'d')
                        for f = 1:Nf(m)
                            Q(k) = Q(k) - spm_dot(wD{m,f},x{m,f}(:,1,k));
                        end
                    end
                    
                    for j = t:S
                        
                        % get expected states for this policy and time
                        %--------------------------------------------------
                        for f = 1:Nf(m)
                            xq{m,f} = x{m,f}(:,j,k);
                        end
                        
                        % (negative) expected free energy
                        %==================================================
                        
                        % Bayesian surprise about states
                        %--------------------------------------------------
                        Q(k) = Q(k) + spm_MDP_G(A(m,:),xq(m,:));
                        
                        for g = 1:Ng(m)
                            
                            % prior preferences about outcomes
                            %----------------------------------------------
                            qo   = spm_dot(A{m,g},xq(m,:));
                            Q(k) = Q(k) + qo'*(C{m,g}(:,j));
                            
                            % Bayesian surprise about parameters
                            %----------------------------------------------
                            if isfield(MDP,'a')
                                Q(k) = Q(k) - spm_dot(wA{m,g},{qo xq{m,:}});
                            end
                        end
                    end
                end
                
                
                % variational updates - policies and precision
                %==========================================================
                
                % previous expected precision
                %----------------------------------------------------------
                if t > 1
                    w{m}(t) = w{m}(t - 1);
                end
                for i = 1:Ni
                    
                    % posterior and prior beliefs about policies
                    %------------------------------------------------------
                    qu = spm_softmax(qE{m}(p{m}) + w{m}(t)*Q(p{m}) + F(p{m}));
                    pu = spm_softmax(qE{m}(p{m}) + w{m}(t)*Q(p{m}));
                    
                    % precision (w) with free energy gradients (v = -dF/dw)
                    %------------------------------------------------------
                    if OPTIONS.gamma
                        w{m}(t) = 1/beta;
                    else
                        eg      = (qu - pu)'*Q(p{m});
                        dFdg    = qb{m} - beta + eg;
                        qb{m}   = qb{m} - dFdg/2;
                        w{m}(t) = 1/qb{m};
                    end
                    
                    % simulated dopamine responses (expected precision)
                    %------------------------------------------------------
                    n             = (t - 1)*Ni + i;
                    wn{m}(n,1)    = w{m}(t);
                    un{m}(p{m},n) = qu;
                    u{m}(p{m},t)  = qu;
                    
                end               
            end % end of loop over multiple policies
            
            
            % Bayesian model averaging of hidden states (over policies)
            %--------------------------------------------------------------
            for f = 1:Nf(m)
                for i = 1:S
                    X{m,f}(:,i) = reshape(x{m,f}(:,i,:),Ns(m,f),Np(m))*u{m}(:,t);
                end
            end
            
            % processing (i.e., reaction) time
            %--------------------------------------------------------------
            rt{m}(t)      = toc(tstart);
            
            % record (negative) free energies
            %--------------------------------------------------------------
            MDP(m).F(:,t) = F;
            MDP(m).G(:,t) = Q;
            MDP(m).H(1,t) = qu'*MDP(m).F(p{m},t) - qu'*(log(qu) - log(pu));
            
            % check for residual uncertainty (in hierarchical schemes)
            %--------------------------------------------------------------
            if isfield(MDP,'factor')
                
                for f = MDP(m).factor(:)'
                    qx     = X{m,f}(:,1);
                    H(m,f) = qx'*spm_log(qx);
                end
                
                % break if there is no further uncertainty to resolve
                %----------------------------------------------------------
                if sum(H(:)) > - chi
                    T = t;
                end
            end
            
            
            % action selection
            %==============================================================
            if t < T
                
                % marginal posterior over action (for each modality)
                %----------------------------------------------------------
                Pu    = zeros([Nu(m,:),1]);
                for i = 1:Np(m)
                    sub        = num2cell(V{m}(t,i,:));
                    Pu(sub{:}) = Pu(sub{:}) + u{m}(i,t);
                end
                
                % action selection (softmax function of action potential)
                %----------------------------------------------------------
                sub            = repmat({':'},1,Nf(m));
                Pu(:)          = spm_softmax(alpha*log(Pu(:)));
                P{m}(sub{:},t) = Pu;
                
                % next action - sampled from marginal posterior
                %----------------------------------------------------------
                try
                    MDP(m).u(:,t) = MDP(m).u(:,t);
                catch
                    ind           = find(rand < cumsum(Pu(:)),1);
                    MDP(m).u(:,t) = spm_ind2sub(Nu(m,:),ind);
                end
                
                % update policy and states for moving policies
                %----------------------------------------------------------
                if isfield(MDP,'U')
                    
                    for f = 1:Nf(m)
                        V{m}(t,:,f) = MDP(m).u(f,t);
                    end
                    for j = 1:size(MDP(m).U,1)
                        if (t + j) < T
                            V{m}(t + j,:,:) = MDP(m).U(j,:,:);
                        end
                    end
                    
                    % and re-initialise expectations about hidden states
                    %------------------------------------------------------
                    for f = 1:Nf(m)
                        for k = 1:Np(m)
                            x{m,f}(:,:,k) = 1/Ns(m,f);
                        end
                    end
                end
                
            end % end of state and action selection
        end % end of variational updates over time
    end % end of loop over models (agents)
    
    % terminate evidence accumulation
    %----------------------------------------------------------------------
    if t == T
        if T == 1
            MDP(m).u = zeros(Nf(m),0);
        end
        if ~HMM
            MDP(m).o  = MDP(m).o(:,1:T);        % outcomes at 1,...,T
            MDP(m).s  = MDP(m).s(:,1:T);        % states   at 1,...,T
            MDP(m).u  = MDP(m).u(:,1:T - 1);    % actions  at 1,...,T - 1
        end
        break;
    end
    
end % end of loop over time

% learning – accumulate concentration parameters
%==========================================================================
for m = 1:size(MDP,1)
    
    for t = 1:T
        
        % mapping from hidden states to outcomes: a
        %------------------------------------------------------------------
        if isfield(MDP,'a')
            for g = 1:Ng(m)
                da     = O{m}(g,t);
                for  f = 1:Nf(m)
                    da = spm_cross(da,X{m,f}(:,t));
                end
                da     = da.*(MDP(m).a{g} > 0);
                MDP(m).a{g} = MDP(m).a{g} + da*eta;
            end
        end
        
        % mapping from hidden states to hidden states: b(u)
        %------------------------------------------------------------------
        if isfield(MDP,'b') && t > 1
            for f = 1:Nf(m)
                for k = 1:Np(m)
                    v   = V{m}(t - 1,k,f);
                    db  = u{m}(k,t)*x{m,f}(:,t,k)*x{m,f}(:,t - 1,k)';
                    db  = db.*(MDP(m).b{f}(:,:,v) > 0);
                    MDP(m).b{f}(:,:,v) = MDP(m).b{f}(:,:,v) + db*eta;
                end
            end
        end
        
        % accumulation of prior preferences: (c)
        %------------------------------------------------------------------
        if isfield(MDP,'c')
            for g = 1:Ng(m)
                dc = O{m}(g,t);
                if size(MDP(m).c{g},2) > 1
                    dc = dc.*(MDP(m).c{g}(:,t) > 0);
                    MDP(m).c{g}(:,t) = MDP(m).c{g}(:,t) + dc*eta;
                else
                    dc = dc.*(MDP(m).c{g}>0);
                    MDP(m).c{g} = MDP(m).c{g} + dc*eta;
                end
            end
        end
    end
    
    % initial hidden states:
    %----------------------------------------------------------------------
    if isfield(MDP,'d')
        for f = 1:Nf(m)
            i = MDP(m).d{f} > 0;
            MDP(m).d{f}(i) = MDP(m).d{f}(i) + X{m,f}(i,1);
        end
    end
    
    % policies
    %----------------------------------------------------------------------
    if isfield(MDP,'e')
        MDP(m).e = MDP(m).e + eta*u{m}(:,T);
    end
    
    % (negative) free energy of parameters (complexity): outcome specific
    %----------------------------------------------------------------------
    for g = 1:Ng(m)
        if isfield(MDP,'a')
            MDP(m).Fa(g) = - spm_KL_dir(MDP(m).a{g},pA{m,g});
        end
        if isfield(MDP,'c')
            MDP(m).Fc(f) = - spm_KL_dir(MDP(m).c{g},pC{g});
        end
    end
    
    % (negative) free energy of parameters: state specific
    %----------------------------------------------------------------------
    for f = 1:Nf(m)
        if isfield(MDP,'b')
            MDP(m).Fb(f) = - spm_KL_dir(MDP(m).b{f},pB{m,f});
        end
        if isfield(MDP,'d')
            MDP(m).Fd(f) = - spm_KL_dir(MDP(m).d{f},pD{m,f});
        end
    end
    
    % (negative) free energy of parameters: policy specific
    %----------------------------------------------------------------------
    if isfield(MDP,'e')
        MDP(m).Fe = - spm_KL_dir(MDP(m).e,pE{m});
    end
    
    % simulated dopamine (or cholinergic) responses
    %----------------------------------------------------------------------
    if Np(m) > 1
        dn{m} = 8*gradient(wn{m}) + wn{m}/8;
    else
        dn{m} = [];
        wn{m} = [];
    end
    
    % Bayesian model averaging of expected hidden states over policies
    %----------------------------------------------------------------------
    for f = 1:Nf(m)
        Xn{m,f} = zeros(Ni,Ns(m,f),T,T);
        Vn{m,f} = zeros(Ni,Ns(m,f),T,T);
        for i = 1:T
            for k = 1:Np(m)
                Xn{m,f}(:,:,:,i) = Xn{m,f}(:,:,:,i) + xn{m,f}(:,:,1:T,i,k)*u{m}(k,i);
                Vn{m,f}(:,:,:,i) = Vn{m,f}(:,:,:,i) + vn{m,f}(:,:,1:T,i,k)*u{m}(k,i);
            end
        end
    end
    
    % use penultimate beliefs about moving policies
    %----------------------------------------------------------------------
    if isfield(MDP,'U')
        u{m}(:,T)  = [];
        try un{m}(:,(end - Ni + 1):end) = []; catch, end
    end
    
    % assemble results and place in NDP structure
    %----------------------------------------------------------------------
    MDP(m).T  = T;            % number of belief updates
    MDP(m).V  = V{m};         % policies
    MDP(m).O  = O{m};         % policies
    MDP(m).P  = P{m};         % probability of action at time 1,...,T - 1
    MDP(m).R  = u{m};         % conditional expectations over policies
    MDP(m).Q  = x(m,:);       % conditional expectations over N states
    MDP(m).X  = X(m,:);       % Bayesian model averages over T outcomes
    MDP(m).C  = C(m,:);       % utility
    
    if HMM
        MDP(m).o  = zeros(Ng(m),0);      % outcomes at 1,...,T
        MDP(m).s  = zeros(Nf(m),0);      % states   at 1,...,T
        MDP(m).u  = zeros(Nf(m),0);      % actions  at 1,...,T - 1
        return
    end
    
    MDP(m).w  = w{m};         % posterior expectations of precision (policy)
    MDP(m).vn = Vn(m,:);      % simulated neuronal prediction error
    MDP(m).xn = Xn(m,:);      % simulated neuronal encoding of hidden states
    MDP(m).un = un{m};        % simulated neuronal encoding of policies
    MDP(m).wn = wn{m};        % simulated neuronal encoding of precision
    MDP(m).dn = dn{m};        % simulated dopamine responses (deconvolved)
    MDP(m).rt = rt{m};        % simulated reaction time (seconds)
    
end


% plot
%==========================================================================
if OPTIONS.plot
    if ishandle(OPTIONS.plot)
        figure(OPTIONS.plot); clf
    else
        spm_figure('GetWin','MDP'); clf
    end
    spm_MDP_VB_trial(MDP(1))
end


% auxillary functions
%==========================================================================

function A  = spm_log(A)
% log of numeric array plus a small constant
%--------------------------------------------------------------------------
A  = log(A + 1e-16);

function A  = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A           = bsxfun(@rdivide,A,sum(A,1));
A(isnan(A)) = 1/size(A,1);

function A  = spm_wnorm(A)
% summation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A   = A + 1e-16;
A   = bsxfun(@minus,1./sum(A,1),1./A);

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

% NOTES:

% generate least surprising outcome
%==========================================================================

% or least surprising outcome
%------------------------------------------------------
% j     = sub2ind(Ns(m,:),ind{:});
% F     = zeros(No(m,g),1);
% for i = 1:No(m,g)
%     po    = MDP(m).A{g}(i,:);
%     po    = spm_norm(spm_vec(po));
%     F(i)  = spm_log(po(j));
% end
% po        = spm_softmax(F*512);

