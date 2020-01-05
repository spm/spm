function [MDP] = spm_MDP_VB_XX(MDP,OPTIONS)
% active inference and learning using belief propagation
% FORMAT [MDP] = spm_MDP_VB_XX(MDP,OPTIONS)
%
% Input; MDP(m,n)       - structure array of m models over n epochs
% MDP.U(P,F)            - P allowable actions over F factors
% MDP.T                 - number of outcomes
%
% MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes given hidden states
% MDP.B{F}(NF,NF,PF)     - transitions among states under PF control states
% MDP.C{G}(O,T)         - (log) prior preferences for outcomes (modality G)
% MDP.D{F}(NF,1)        - prior probabilities over initial states
% MDP.E(P,1)            - prior probabilities over policies
%
% MDP.a{G}              - concentration parameters for A
% MDP.b{F}              - concentration parameters for B
% MDP.c{G}              - concentration parameters for C
% MDP.d{F}              - concentration parameters for D
% MDP.e{P}              - concentration parameters for E
%
% optional:
% MDP.s(F,T)            - matrix of true states - for each hidden factor
% MDP.o(G,T)            - matrix of outcomes    - for each outcome modality
% or .O{G}(O,T)         - likelihood matrix     - for each outcome modality
% MDP.u(F,T - 1)        - vector of actions     - for each hidden factor
%
% MDP.alpha             - precision - action selection [512]
% MDP.chi               - Occams window for deep updates
% MDP.eta               - learning rate for model parameters
% MDP.N                 - depth of deep policy search [N <= T]
%
% MDP.demi.C            - Mixed model: cell array of true causes (DEM.C)
% MDP.demi.U            - Bayesian model average (DEM.U) see: spm_MDP_DEM
% MDP.link              - link array to generate outcomes from
%                         subordinate MDP; for deep (hierarchical) models
%
% OPTIONS.plot          - switch to suppress graphics:  (default: [0])
% OPTIONS.D             - switch to update initial states over epochs
% OPTIONS.BMR           - Bayesian model reduction for multiple trials
%                         see: spm_MDP_VB_sleep(MDP,BMR)
% Outputs:
%
% MDP.P(N1,...,NF,T)    - action probability
% MDP.X{F}(NF,T)        - conditional expectation is over hidden states
% MDP.R(P,T)            - conditional expectations over policies
%
% MDP.un          - simulated neuronal encoding of hidden states
% MDP.xn          - simulated neuronal encoding of policies
% MDP.wn          - simulated neuronal encoding of precision (tonic)
% MDP.dn          - simulated dopamine responses (phasic)
%
% This routine provides solutions of active inference (minimisation of
% variational free energy) using a generative model based upon a Markov
% decision process. The model and inference scheme is formulated in
% discrete space and time. This means that the generative model (and
% process) are hidden Markov models whose dynamics are given by transition
% probabilities among states and the likelihood corresponds to a particular
% outcome conditioned upon hidden states.
%
% This implementation equips agents with the prior beliefs that they will
% maximise expected free energy. Variational free energy can be interpreted
% in several ways - most intuitively as minimising the KL divergence
% between predicted and preferred outcomes (specified as prior beliefs) -
% while simultaneously minimising ambiguity.
%
% This particular scheme is designed for any allowable policies or control
% variables specified in MDP.U. Constraints on allowable policies can limit
% the numerics or combinatorics considerably. Further, the outcome space
% and hidden states can be defined in terms of factors; corresponding to
% sensory modalities and (functionally) segregated representations,
% respectively. This means, for each factor or subset of hidden states
% there are corresponding control states that determine the transition
% probabilities. in this implementation, hidden factors are combined using
% a Kronecker intensive product to enable exact Bayesian inference using
% belief propagation (the Kronecker tensor form ensures that conditional
% dependencies among hidden factors are evaluated).
%
% In this belief propagation scheme, the next action is evaluated in terms
% of the free energy expected under all subsequent actions until some time
% horizon (specified by MDP.T). This expected free energy is accumulated
% along all allowable paths or policies (see the subroutine spm_forward);
% effectively, performing a deep tree search over future sequences of
% actions. Because actions are conditionally independent of previous
% actions, it is only necessary to update posterior beliefs over hidden
% states at the current time point (using a Bayesian belief updating) and
% then use the prior over actions (based upon expected free energy) to
% select the next action. Previous actions are realised variables and are
% used when evaluating the posterior beliefs over current states.
%
% In brief, the agent encodes beliefs about hidden states in the past
% conditioned on realised outcomes and actions. The resulting conditional
% expectations determine the (path integral) of free energy that then
% determines an empirical prior over the next action, from which the next
% realised action sampled
%
% In addition to state estimation and policy selection, the scheme also
% updates model parameters; including the state transition matrices,
% mapping to outcomes and the initial state. This is useful for learning
% the context. Likelihood and prior probabilities can be specified in terms
% of concentration parameters (of a Dirichlet distribution (a,b,c,..). If
% the corresponding (A,B,C,..) are supplied, they will be used to generate
% outcomes.
%
% If supplied with a structure array, this routine will automatically step
% through the implicit sequence of epochs (implicit in the number of
% columns of the array). If the array has multiple rows, each row will be
% treated as a separate model or agent. This enables agents to communicate
% through acting upon a common set of hidden factors, or indeed sharing the
% same outcomes.
%
% See also: spm_MDP, which uses multiple future states and a mean field
% approximation for control states - but allows for different actions at
% all times (as in control problems).
%
% See also: spm_MDP_VB_X,  which is the corresponding variational message
% passing scheme for fixed policies; i.e., ordered sequences of actions
% that are specified a priori.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_XX.m 7766 2020-01-05 21:37:39Z karl $


% deal with a sequence of trials
%==========================================================================

% options
%--------------------------------------------------------------------------
try, OPTIONS.plot;  catch, OPTIONS.plot  = 0; end
try, OPTIONS.D;     catch, OPTIONS.D     = 0; end
global COUNT

% check MDP specification
%--------------------------------------------------------------------------
MDP = spm_MDP_check(MDP);

% handle multiple trials, ensuring parameters (and posteriors) are updated
%==========================================================================
if size(MDP,2) > 1
    
    % plotting options
    %----------------------------------------------------------------------
    GRAPH        = OPTIONS.plot;
    OPTIONS.plot = 0;
    
    for i = 1:size(MDP,2)                  % number of MDPs
        for m = 1:size(MDP,1)              % number of trials
            if i > 1                       % if previous inversions
                
                % update concentration parameters
                %----------------------------------------------------------
                MDP(m,i)  = spm_MDP_update(MDP(m,i),OUT(m,i - 1));
                
                % update initial states (post-diction)
                %----------------------------------------------------------
                if any(OPTIONS.D)
                    nD = numel(MDP(m,i).D);
                    if numel(OPTIONS.D) ~= nD
                        OPTIONS.D = ones(nD,1);
                    end
                    for f = 1:nD
                        if OPTIONS.D(f)
                            MDP(m,i).D{f} = OUT(m,i - 1).X{f}(:,end);
                        end
                    end
                end
            end
        end
        
        % solve this trial (for all models synchronously)
        %------------------------------------------------------------------
        OUT(:,i) = spm_MDP_VB_XX(MDP(:,i),OPTIONS);
        
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

% number of outcomes T and policies U (and V for consistency with plotting
% schemes)
%--------------------------------------------------------------------------
[T,U,V] = spm_MDP_get_T(MDP);

% defaults
%--------------------------------------------------------------------------
try, alpha = MDP(1).alpha; catch, alpha = 512;  end % action precision
try, eta   = MDP(1).eta;   catch, eta   = 1;    end % learning rate
try, chi   = MDP(1).chi;   catch, chi   = 1/64; end % Occam window updates
try, N     = MDP(1).N;     catch, N     = T;    end % depth of policy search
N          = min(N,T);

% initialise model-specific parameters
%==========================================================================
for m = 1:size(MDP,1)
    
    % numbers of transitions, policies and states
    %----------------------------------------------------------------------
    Ng(m) = numel(MDP(m).A);               % number of outcome factors
    Nf(m) = numel(MDP(m).B);               % number of hidden state factors
    Np(m) = size(U{m},1);                  % number of allowable actions
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
        
        % prior concentration parameters and novelty (W)
        %------------------------------------------------------------------
        if isfield(MDP,'a')
            pA{m,g} = MDP(m).a{g};
            W{m,g}  = spm_wnorm(pA{m,g}).*(pA{m,g} > 0);
            W{m,g}  = W{m,g}(:,:);
        else
            W{m,g}  = spm_zeros(A{m,g});
            W{m,g}  = W{m,g}(:,:);
        end
        
        % and ambiguity (H) (for computation of expected free energy: G)
        %------------------------------------------------------------------
        H{m,g}      = sum(A{m,g}.*spm_log(A{m,g}),1);
        H{m,g}      = spm_vec(shiftdim(H{m,g},1));
 
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
            if isfield(MDP,'b')
                fB{m,f}(:,:,j) = spm_norm(MDP(m).b{f}(:,:,j) );
            else
                fB{m,f}(:,:,j) = spm_norm(MDP(m).B{f}(:,:,j) );
            end
            
        end
        
        % prior concentration paramters for novelty
        %------------------------------------------------------------------
        if isfield(MDP,'b')
            pB{m,f} = MDP(m).b{f};
        end
        
    end
    
    % Assemble Kronecker form of policies (for vectorised states)
    %----------------------------------------------------------------------
    for k = 1:Np(m)
        
        % belief propagation (B): tensor product over hidden factors
        %------------------------------------------------------------------
        B{m,k} = 1;
        for f = 1:Nf(m)
            B{m,k} = spm_kron(fB{m,f}(:,:,U{m}(k,f)),B{m,k});
        end
    end
    
    % priors over initial hidden states: concentration parameters
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
        
        % prior concentration paramters for novelty
        %------------------------------------------------------------------
        if isfield(MDP,'d')
            pD{m,f} = MDP(m).d{f};
        end
    end
    
    % priors over policies: concentration parameters
    %----------------------------------------------------------------------
    if isfield(MDP,'e')
        E{m} = spm_norm(MDP(m).e);
    elseif isfield(MDP,'E')
        E{m} = spm_norm(MDP(m).E);
    else
        E{m} = spm_norm(ones(Np(m),1));
    end
    E{m}     = spm_log(E{m});
    
    % prior concentration paramters for novellty
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
    
    
    % initialise posterior expectations (X) of hidden states at the current
    % time and over time (S)
    %======================================================================
    for f = 1:Nf(m)
        S{m,f} = zeros(Ns(m,f),T,T) + 1/Ns(m,f);
        X{m,f} = repmat(D{m,f},1,T);
    end
    for t = 1:T
        Q{m,t} = spm_kron(D(m,:));
    end
    
    % initialise posteriors over polices and action
    %----------------------------------------------------------------------
    P{m}  = zeros([Nu(m,:),1]);
    
    % if states have not been specified, set to 0
    %----------------------------------------------------------------------
    s{m}  = zeros(Nf(m),T);
    try
        i       = find(MDP(m).s);
        s{m}(i) = MDP(m).s(i);
    end
    MDP(m).s    = s{m};
    
    % if outcomes have not been specified set to 0
    %----------------------------------------------------------------------
    o{m}  = zeros(Ng(m),T);
    try
        i       = find(MDP(m).o);
        o{m}(i) = MDP(m).o(i);
    end
    MDP(m).o    = o{m};
    
end

% ensure any outcome generating agent is updated first
%--------------------------------------------------------------------------
[M,MDP] = spm_MDP_get_M(MDP,T,Ng);


% belief updating over successive time points
%==========================================================================
for t = 1:T
    
    % generate hidden states and outcomes for each agent or model
    %======================================================================
    for m = M(t,:)
        
        % sample state, if not specified
        %------------------------------------------------------------------
        for f = 1:Nf(m)
            
            % the next state is generated by action on external states
            %--------------------------------------------------------------
            if MDP(m).s(f,t) == 0
                if t > 1
                    ps = MDP(m).B{f}(:,MDP(m).s(f,t - 1),MDP(m).u(f,t - 1));
                else
                    ps = spm_norm(MDP(m).D{f});
                end
                MDP(m).s(f,t) = find(rand < cumsum(ps),1);
            end
            
        end
        
        % posterior predictive density over hidden (external) states
        %------------------------------------------------------------------
        for f = 1:Nf(m)
            
            % under realised action (xq)
            %--------------------------------------------------------------
            if t > 1
                xq{m,f} = fB{m,f}(:,:,MDP(m).u(f,t - 1))*X{m,f}(:,t - 1);
            else
                xq{m,f} = X{m,f}(:,t);
            end
            
        end
        
        % sample outcome, if not specified
        %------------------------------------------------------------------
        for g = 1:Ng(m)
            
            % if outcome is not specified
            %--------------------------------------------------------------
            if ~MDP(m).o(g,t)
                
                % outcome is generated by model n
                %----------------------------------------------------------
                if MDP(m).n(g,t)
                    
                    n    = MDP(m).n(g,t);
                    if n == m
                        
                        % outcome that minimises free energy (i.e.,
                        % maximises accuracy)
                        %----------------------------------------------
                        F             = spm_dot(spm_log(A{m,g}),xq(m,:));
                        po            = spm_softmax(F*512);
                        MDP(m).o(g,t) = find(rand < cumsum(po),1);

                    else
                        
                        % outcome from model n
                        %--------------------------------------------------
                        MDP(m).o(g,t) = MDP(n).o(g,t);
                        
                    end
                    
                else
                    
                    % or sample from likelihood given hidden state
                    %------------------------------------------------------
                    ind           = num2cell(MDP(m).s(:,t));
                    po            = MDP(m).A{g}(:,ind{:});
                    MDP(m).o(g,t) = find(rand < cumsum(po),1);
                    
                end
            end
        end
        
        % get probabilistic outcomes from samples or subordinate level
        %==================================================================
        
        % get outcome probability (O{m})
        %------------------------------------------------------------------
        for g = 1:Ng(m)
            O{m}{g,t} = sparse(MDP(m).o(g,t),1,1,No(m,g),1);
        end
        
        % or generate outcomes from a subordinate MDP
        %==================================================================
        if isfield(MDP,'link')
            
            % use previous inversions (if available) to reproduce outcomes
            %--------------------------------------------------------------
            try
                mdp = MDP(m).mdp(t);
            catch
                try
                    mdp     = spm_MDP_update(MDP(m).MDP(t),MDP(m).mdp(t - 1));
                catch
                    try
                        mdp = spm_MDP_update(MDP(m).MDP(1),MDP(m).mdp(t - 1));
                    catch
                        mdp = MDP(m).MDP(1);
                    end
                end
            end
            
            % priors over states (of subordinate level)
            %--------------------------------------------------------------
            mdp.factor = [];
            for f = 1:size(MDP(m).link,1)
                for g = 1:size(MDP(m).link,2)
                    if ~isempty(MDP(m).link{f,g})
                        
                        % subordinate state has hierarchical constraints
                        %--------------------------------------------------
                        mdp.factor(end + 1) = f;
                        
                        % empirical priors over initial states
                        %--------------------------------------------------
                        O{m}{g,t} = spm_dot(A{m,g},xq(m,:));
                        mdp.D{f}  = MDP(m).link{f,g}*O{m}{g,t};
                        
                        % outcomes (i.e., states) are generated by model n
                        %--------------------------------------------------
                        if MDP(m).n(g,t)
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
            
            
            % empirical prior preferences
            %--------------------------------------------------------------
            if isfield(MDP,'linkC')
                for f = 1:size(MDP(m).linkC,1)
                    for g = 1:size(MDP(m).linkC,2)
                        if ~isempty(MDP(m).linkC{f,g})
                            O{m}{g,t} = spm_dot(A{m,g},xq(m,:));
                            mdp.C{f}  = spm_log(MDP(m).linkC{f,g}*O{m}{g,t});
                        end
                    end
                end
            end
            
            % empirical priors over policies
            %--------------------------------------------------------------
            if isfield(MDP,'linkE')
                mdp.factorE = [];
                for g = 1:size(MDP(m).linkE,2)
                    if ~isempty(MDP(m).linkE{g})
                        O{m}{g,t} = spm_dot(A{m,g},xq(m,:));
                        mdp.E     = MDP(m).linkE{g}*O{m}{g,t};
                    end
                end
            end
            
            
            
            % infer hidden states at lower level (outcomes at this level)
            %==============================================================
            MDP(m).mdp(t) = spm_MDP_VB_XX(mdp);
            
            
            % get inferred outcomes from subordinate MDP
            %==============================================================
            for f = 1:size(MDP(m).link,1)
                for g = 1:size(MDP(m).link,2)
                    if ~isempty(MDP(m).link{f,g})
                        O{m}{g,t} = MDP(m).link{f,g}'*MDP(m).mdp(t).X{f}(:,1);
                    end
                end
            end
            
            % if hierarchical preferences, these contribute to outcomes ...
            %--------------------------------------------------------------
            if isfield(MDP,'linkC')
                for f = 1:size(MDP(m).linkC,1)
                    for g = 1:size(MDP(m).linkC,2)
                        if ~isempty(MDP(m).linkC{f,g})
                            indC      = sparse(MDP(m).mdp(t).o(f,:)',1:length(MDP(m).mdp(t).o(f,:)),ones(length(MDP(m).mdp(t).o(f,:)),1),size(MDP(m).mdp(t).C{f},1),size(MDP(m).mdp(t).C{f},2));
                            O{m}{g,t} = spm_softmax(spm_log(O{m}{g,t}) + MDP(m).linkC{f,g}'*sum((indC.*(MDP(m).mdp(t).C{f})),2));
                        end
                    end
                end
            end
            
            % ... and the same for policies
            %--------------------------------------------------------------
            if isfield(MDP,'linkE')
                for g = 1:size(MDP(m).linkE,2)
                    if ~isempty(MDP(m).linkE{g})
                        O{m}{g,t} = spm_softmax(spm_log(O{m}{g,t}) + spm_log(MDP(m).linkE{g}'*MDP(m).mdp(t).R(:,end)));
                    end
                end
            end
            
            % Ensure DEM starts with final states from previous inversion
            %--------------------------------------------------------------
            if isfield(MDP(m).MDP,'demi')
                MDP(m).MDP.DEM.G(1).x = MDP(m).mdp(t).dem(end).pU.x{1}(:,end);
                MDP(m).MDP.DEM.M(1).x = MDP(m).mdp(t).dem(end).qU.x{1}(:,end);
            end
            
        end % end of hierarchical mode (link)
        
        
        % or generate outcome likelihoods from a variational filter
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
                O{m}{g,t} = spm_dot(A{m,g},xq(m,:));
            end
            
            % get posterior outcome from Bayesian filtering
            %--------------------------------------------------------------
            MDP(m).dem(t) = spm_MDP_DEM(MDP(m).dem(t),...
                MDP(m).demi,O{m}(:,t),MDP(m).o(:,t));
            
            for g = 1:Ng(m)
                O{m}{g,t} = MDP(m).dem(t).X{g}(:,end);
            end
            
        end % end outcomes from Bayesian filter
        
        
        % or generate outcome likelihoods from voice recognition
        %==================================================================
        if isfield(MDP,'VOX')
            
            % get predictive prior over outcomes if MDP.VOX = 2
            %--------------------------------------------------------------
            if MDP(m).VOX == 2
                
                % current outcome
                %----------------------------------------------------------
                for g = 1:Ng(m)
                    O{m}{g,t} = spm_dot(A{m,g},xq(m,:));
                end
                
                % and next outcome if available
                %----------------------------------------------------------
                try
                    for f = 1:Nf(m)
                        pq{f} = fB{m,f}(:,:)*xq{m,f};
                    end
                    for g = 1:Ng(m)
                        O{m}{g,t + 1} = spm_dot(A{m,g},pq);
                    end
                end
                
            end
            
            % get likelihood over outcomes - or articulate phrase
            %--------------------------------------------------------------
            O{m}  = spm_MDP_VB_VOX(MDP(m),O{m},t);
            
            % update outcomes
            %--------------------------------------------------------------
            for g = 1:Ng(m)
                po            = spm_softmax(O{m}{g,t}*512);
                MDP(m).o(g,t) = find(rand < cumsum(po),1);
            end
            
        end % end outcomes from voice recognition
        
        
        % Bayesian belief updating hidden states (Q) and policies (G)
        %==================================================================
        
        % empirical prior over hidden states at this time
        %------------------------------------------------------------------
        if t > 1
            
            % use transition probabilities (using current action)
            %--------------------------------------------------------------
            Q{m,t} = B{m,K(m,t - 1)}*Q{m,t - 1};
            
        else
            % use empirical priors over initial states
            %--------------------------------------------------------------
            Q{m,t} = spm_kron(D(m,:));
            
        end
        
        
        %  posterior over hidden states (Q) and expected free energy (G)
        %==================================================================
        COUNT = 0;
        [G,Q(m,:)] = spm_forwards(...
            O{m}(:,t),Q(m,:),A(m,:),B(m,:),C(m,:),E{m},H(m,:),W(m,:),t,T,min(T,t + N));
        
        % save marginal posteriors over hidden states
        %------------------------------------------------------------------
        for i = 1:T
            qx    = reshape(Q{m,i},[Ns,1]);
            for f = 1:Nf(m)
                S{m,f}(:,i,t) = spm_margin(qx,f);
            end
        end
        for f = 1:Nf(m)
            X{m,f}(:,t) = S{m,f}(:,t,t);
        end
        
        % posterior beliefs about policies (u) and precision (w)
        %------------------------------------------------------------------
        u{m,t}    = spm_softmax(G);
        w{m}(t)   = u{m,t}'*spm_log(u{m,t});
        
        % end policy search
        %==================================================================
        % disp(COUNT)
        
        
        % check for residual uncertainty (in hierarchical schemes)
        %------------------------------------------------------------------
        if isfield(MDP,'factor')
            
            for f = MDP(m).factor(:)'
                qx     = X{m,f}(:,1);
                S(m,f) = qx'*spm_log(qx);
            end
            
            % break if there is no further uncertainty to resolve
            %--------------------------------------------------------------
            if sum(S(:)) > - chi && ~isfield(MDP,'VOX')
                T = t;
            end
        end
        
        % check for end of sentence (' ') if in VOX mode
        %------------------------------------------------------------------
        XVOX = 0;
        if isfield(MDP,'VOX') && t > 1
            if strcmp(MDP.label.outcome{1}{MDP(m).o(1,t)},' ')
                T    = t;
                XVOX = 1;
            end
        end
        
        % action selection
        %==================================================================
        if t < T
            
            % record emprical prior over policies (R) and sample a policy (K)
            %--------------------------------------------------------------
            R{m}(:,t)      = u{m,t};
            Ru             = spm_softmax(alpha*log(u{m,t}));
            K(m,t)         = find(rand < cumsum(Ru),1);
            
            % marginal posterior (P) over action (for each factor)
            %--------------------------------------------------------------
            Pu    = zeros([Nu(m,:),1]);
            for k = 1:Np(m)
                sub        = num2cell(U{m}(k,:));
                Pu(sub{:}) = Pu(sub{:}) + Ru(k);
            end
            sub            = repmat({':'},1,Nf(m));
            P{m}(sub{:},t) = Pu;
            
            % realised action
            %--------------------------------------------------------------
            MDP(m).v(:,t)     = K(m,t);
            try
                MDP(m).u(:,t) = MDP(m).u(:,t);
            catch
                MDP(m).u(:,t) = U{m}(K(m,t),:);
            end
            
        end % end of state and action selection
    end % end of loop over models (agents)
    
    
    % Evaluate reporting function specified
    %======================================================================
    if isfield(MDP,'FCN')
        try
            MDP.FCN(MDP,X);
        end
    end
    
    % terminate evidence accumulation
    %----------------------------------------------------------------------
    if t == T
        MDP(m).o  = MDP(m).o(:,1:T);        % outcomes at 1,...,T
        MDP(m).s  = MDP(m).s(:,1:T);        % states   at 1,...,T
        MDP(m).u  = MDP(m).u(:,1:T - 1);    % actions  at 1,...,T - 1
        break;
    end
    
end % end of loop over time

% loop over models to accumulate Dirichlet parameters and prepare outputs
%==========================================================================
for m = 1:size(MDP,1)
    
    
    % learning - accumulate concentration parameters
    %======================================================================
    for t = 1:T
        
        % mapping from hidden states to outcomes: a
        %------------------------------------------------------------------
        if isfield(MDP,'a')
            for g = 1:Ng(m)
                da = spm_cross(O{m}(g,t),Q{m,t});
                da = reshape(da,[No(g),Ns]);
                da = da.*(MDP(m).a{g} > 0);
                MDP(m).a{g} = MDP(m).a{g} + da*eta;
            end
        end
        
        % mapping from hidden states to hidden states: b(u)
        %------------------------------------------------------------------
        if isfield(MDP,'b') && t < T
            for f = 1:Nf(m)
                j   = U{m}(K(m,t),f);
                db  = spm_cross(X{m,f}(:,t + 1),X{m,f}(:,t)');
                db  = db.*(MDP(m).b{f}(:,:,j) > 0);
                MDP(m).b{f}(:,:,j) = MDP(m).b{f}(:,:,j) + db*eta;
            end
        end
        
        % accumulation of prior preferences: (c)
        %------------------------------------------------------------------
        if isfield(MDP,'c') && t < T
            for g = 1:Ng(m)
                dc = O{m}{g,t + 1};
                if size(MDP(m).c{g},2) > 1
                    dc = dc.*(MDP(m).c{g}(:,t) > 0);
                    MDP(m).c{g}(:,t) = MDP(m).c{g}(:,t) + dc*eta;
                else
                    dc = dc.*(MDP(m).c{g} > 0);
                    MDP(m).c{g} = MDP(m).c{g} + dc*eta;
                end
            end
        end
    end
    
    % initial hidden states:
    %----------------------------------------------------------------------
    if isfield(MDP,'d')
        
        % posterior over initial states
        %------------------------------------------------------------------
        L     = spm_backwards(O{m},Q(m,:),A(m,:),B(m,:),K(m,:),T);
        L     = reshape(L,Ns);
        for f = 1:Nf(m)
            x{f} = spm_margin(L,f);
        end
        
        %  accumulate Dirichlet counts
        %------------------------------------------------------------------
        for f = 1:Nf(m)
            i = MDP(m).d{f} > 0;
            MDP(m).d{f}(i) = MDP(m).d{f}(i) + x{f}(i);
        end
    end
    
    % policies
    %----------------------------------------------------------------------
    if isfield(MDP,'e')
        de   = 0;
        for t = 1:(T - 1)
            de = de + u{m,t};
        end
        MDP(m).e = MDP(m).e + de*eta;
    end
    
    % (negative) free energy of parameters (complexity): outcome specific
    %======================================================================
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
    
    
    % simulated electrophysiological responses
    %======================================================================
    
    % simulated dopamine (or cholinergic) responses: assuming a
    % monoexponential kernel
    %----------------------------------------------------------------------
    n     = 16;
    h     = exp(-(0:(n - 1))/2);
    h     = h/sum(h);
    wn{m} = kron(w{m},ones(1,n));
    wn{m} = conv(wn{m},[spm_zeros(h) h],'same');
    dn{m} = gradient(wn{m}(:));
    
    
    % Belief updating about hidden states: assuming a kernel or impulse
    % response function with a cumulative gamma distribution
    %----------------------------------------------------------------------
    h     = spm_Gcdf(0:(n - 1),n/4,1);
    for f = 1:Nf(m)
        for i = 1:Ns(m,f)
            for j = 1:T
                for k = 1:T
                    if k == 1
                        h0 = 1/Ns(m,f);
                    else
                        h0 = S{m,f}(i,j,k - 1);
                    end
                    ht     = S{m,f}(i,j,k);
                    xn{m,f}(:,i,j,k) = h*(ht - h0) + h0;
                end
            end
        end
    end
    
    % sum to one contraint
    %----------------------------------------------------------------------
    for i = 1:n
        for j = 1:T
            for k = 1:T
                xn{m,f}(i,:,j,k) = xn{m,f}(i,:,j,k)/sum(xn{m,f}(i,:,j,k));
            end
        end
    end
    
    % belief updating about policies
    %----------------------------------------------------------------------
    u0    = spm_softmax(E{m});
    for k = 1:Np(m)
        for t = 1:(T - 1)
            if t == 1
                h0 = u0(k);
            else
                h0 = u{m,t - 1}(k);
            end
            ht         = u{m,t}(k);
            j          = (1:n) + (t - 1)*n;
            un{m}(k,j) = (h*(ht - h0) + h0);
        end
    end
    
    
    % assemble results and place in NDP structure
    %======================================================================
    MDP(m).T  = T;            % number of outcomes
    MDP(m).V  = V{m};         % policies
    MDP(m).O  = O{m};         % outcomes
    MDP(m).P  = P{m};         % probability of action over factors and time
    MDP(m).R  = R{m};         % conditional expectations over policies
    MDP(m).X  = X(m,:);       % conditional expectations over states
    MDP(m).C  = C(m,:);       % utility
    
    MDP(m).w  = w{m};         % precision of beliefs about policies
    MDP(m).xn = xn(m,:);      % simulated neuronal encoding of states
    MDP(m).un = un{m};        % simulated neuronal encoding of policies
    MDP(m).wn = wn{m};        % simulated neuronal encoding of precision
    MDP(m).dn = dn{m};        % simulated dopamine responses (phasic)
    
    
end % end loop over models (m)


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
function [G,P] = spm_forwards(O,P,A,B,C,E,H,W,t,T,N)
% deep tree search over policies or paths
%--------------------------------------------------------------------------
% FORMAT [G,Q] = spm_G(O,P,A,B,C,E,H,W,t,T);
% O{g}   - cell array of outcome probabilities for modality g
% P      - empirical prior over (vectorised) states
% A{g}   - likelihood mappings from hidden states
% B{k}   - belief propagators (action dependent probability transitions)
% C{g}   - cost: log priors over future outcomes)
% E      - empirical prior over actions
% H{g}   - state dependent ambiguity
% W{g}   - state and outcome dependent novelty
% t      - current time point
% T      - time horizon
% N      - policy horizon
%
% Q      - posterior over (vectorised) states for t
%
%  This subroutine performs a deep tree search over sequences of actions to
%  evaluate the expected free energy over policies or paths. Crucially, it
%  only searches likely policies under likely hidden states in the future.
%  This search is sophisticated; in the sense that posterior beliefs are
%  updated on the basis of future outcomes to evaluate the free energy
%  under each outcome. The resulting  average is then accumulated to
%  furnish a path integral of expected free energy for the next action.
%  This routine operates recursively by updating predictive posteriors over
%  hidden states and their most likely outcomes.
%__________________________________________________________________________

global COUNT
COUNT = COUNT + 1;

% Posterior over hidden states based on likelihood (L) and priors (P)
%==========================================================================
G     = E;                                   % log priors over actions
L     = 1;
for g = 1:numel(A)
    L = L.*spm_dot(A{g},O{g});
end
P{t}  = spm_norm(L(:).*P{t});

% terminate search at time horizon
%--------------------------------------------------------------------------
if t == T, return, end

% Expected free energy of subsequent action
%==========================================================================
for k = 1:size(B,2)                       % search over actions
    
    % (negative) expected free energy
    %----------------------------------------------------------------------
    Q{k}   = B{k}*P{t};                   % predictive posterior Q{k}
    for g  = 1:numel(A)
        
        % predictive posterior and prior over outcomes
        %--------------------------------------------------------------
        qo   = A{g}(:,:)*Q{k};
        po   = C{g}(:,t);
        
        if 0
            
            % Bayesian risk only
            %--------------------------------------------------------------
            G(k) = G(k) + qo'*po;
            
        else
            
            % G(k)      = ambiguity  + risk
            %--------------------------------------------------------------
            G(k) = G(k) + Q{k}'*H{g} - qo'*(spm_log(qo) - po);
            
            % Bayesian surprise about parameters (i.e., novelty)
            %--------------------------------------------------------------
            G(k) = G(k) - qo'*W{g}*Q{k};
            
        end
    end
end

% deep (recursive) search over action sequences ( i.e., paths)
%==========================================================================
if t < N
    
    % probability over action (terminating search at a suitable threshold)
    %----------------------------------------------------------------------
    u     = spm_softmax(G);
    k     = u <= 1/16;
    u(k)  = 0;
    G(k)  = - 64;
    for k = 1:size(B,2)                      % search over actions
        if u(k) > 1/16                       % evaluating plausible paths
            
            %  evaluate  expected free energy for plausible hidden states
            %--------------------------------------------------------------
            j     = find(Q{k} > 1/16);
            if isempty(j)
                j = find(Q{k} > 1/numel(Q{k}));
            end
            for i = j(:)'
                
                % outcome probabilities under hidden state (i)
                %----------------------------------------------------------
                for g = 1:numel(A)
                    O{g} = A{g}(:,i);
                end
                
                % prior over subsequent action under this hidden state
                %----------------------------------------------------------
                P{t + 1} = Q{k};
                F        = spm_forwards(O,P,A,B,C,E,H,W,t + 1,T,N);
                
                % expected free energy marginalised over subsequent action
                %----------------------------------------------------------
                K(i)     = spm_softmax(F)'*F;
                
            end
            
            % accumulate expected free energy marginalised over states
            %--------------------------------------------------------------
            G(k) = G(k) + K(j)*Q{k}(j);
            
        end % plausible paths
    end % search over actions
    
    
    % Predictive posterior over hidden states
    %----------------------------------------------------------------------
    u     = spm_softmax(G);
    R     = 0;
    for k = 1:size(B,2)
        R = R + u(k)*Q{k};
    end
    P{t + 1} = R;
    
end


function [L] = spm_backwards(O,Q,A,B,u,T)
% Backwards smoothing to evaluate posterior over initial states
%--------------------------------------------------------------------------
% O{g}   - cell array of outcome probabilities for modality g
% Q{t}   - posterior expectations over vectorised hidden states
% A{g}   - likelihood mappings from hidden states
% B{k}   - belief propagators (action dependent probability transitions)
% u{t}   - posterior expectations over actions
% T      - time horizon
%
% L      - posterior over initial states
%
%  This subroutine evaluate the posterior over initial states using a
%  backwards algorithm; namely, by evaluating the likelihood of each
%  initial state, given subsequent outcomes, under posterior expectations
%  about state transitions.


% initialise to posterior and accumulate likelihoods for each initial state
%--------------------------------------------------------------------------
L     = Q{1};
q     = spm_zeros(L);
for i = 1:numel(q)
    qi    = q;
    qi(i) = 1;
    for t = 2:T
        
        % next distribution over hidden states
        %------------------------------------------------------------------
        qi    = B{u(t - 1)}*qi;
        
        % and accumulate likelihood
        %------------------------------------------------------------------
        for g = 1:numel(A)
            L(i) = L(i).*(O{g,t}'*A{g}(:,:)*qi);
        end
    end
end

% marginal distribution over initial states
%--------------------------------------------------------------------------
L     = spm_norm(L(:));



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
A   = bsxfun(@minus,1./sum(A,1),1./A)/2;


function A  = spm_margin(A,i)
% marginalise a joint distribution
%--------------------------------------------------------------------------
d    = 1:ndims(A);
d(i) = [];
A    = sum(A,d);
A    = A(:);


return


function [T,U,V] = spm_MDP_get_T(MDP)
% FORMAT [T,U,V] = spm_MDP_get_T(MDP)
% returns number of policies, policy cell array and HMM flag
% MDP(m) - structure array of m MPDs
% T      - number of trials or updates
% U(m)   - indices of actions for m-th MDP
% HMM    - flag indicating a hidden Markov model
%
% This subroutine returns the policy matrix as a cell array (for each
% model) and the maximum number of updates. If outcomes are specified
% probabilistically in the field MDP(m).O, and there is only one policy,
% the partially observed MDP reduces to a hidden Markov model.
%__________________________________________________________________________

for m = 1:size(MDP,1)
    
    if isfield(MDP(m),'U')
        
        % called with repeatable actions (U,T)
        %------------------------------------------------------------------
        T(m) = MDP(m).T;                    % number of updates
        U{m} = MDP(m).U;                    % allowable actions (Np,Nf)
        V{m}(1,:,:) = MDP(m).U;             % allowable actions (1,Np,Nf)
        
    elseif isfield(MDP(m),'V')
        
        % full sequential policies (V)
        %------------------------------------------------------------------
        T(m) = MDP(m).T;                    % number of updates
        U{m} = MDP(m).V(1,:,:);             % allowable actions (1,Np,Nf)
        V{m} = MDP(m).V;                    % allowable actions (Np,Nf)
    end
    
end

% number of time steps
%--------------------------------------------------------------------------
T = max(T);

return


function [M,MDP] = spm_MDP_get_M(MDP,T,Ng)
% FORMAT [M,MDP] = spm_MDP_get_M(MDP,T,Ng)
% returns an update matrix for multiple models
% MDP(m) - structure array of m MPDs
% T      - number of trials or updates
% Ng(m)  - number of output modalities for m-th MDP
%
% M      - update matrix for multiple models
% MDP(m) - structure array of m MPDs
%
% In some applications, the outcomes are generated by a particular model
% (to maximise free energy, based upon the posterior predictive density).
% The generating model is specified in the matrix MDP(m).n, with a row for
% each outcome modality, such that each row lists the index of the model
% responsible for generating outcomes.
%__________________________________________________________________________

% check for VOX and ensure the agent generates outcomes when speaking
%--------------------------------------------------------------------------
if numel(MDP) == 1
    if isfield(MDP,'MDP')
        if isfield(MDP.MDP,'VOX')
            MDP.n = [MDP.MDP.VOX] == 1;
        end
    end
end

for m = 1:size(MDP,1)
    
    % check size of outcome generating agent, as specified by MDP(m).n
    %----------------------------------------------------------------------
    if ~isfield(MDP(m),'n')
        MDP(m).n = zeros(Ng(m),T);
    end
    if size(MDP(m).n,1) < Ng(m)
        MDP(m).n = repmat(MDP(m).n(1,:),Ng(m),1);
    end
    if size(MDP(m).n,1) < T
        MDP(m).n = repmat(MDP(m).n(:,1),1,T);
    end
    
    % mode of generating model (most frequent over outcome modalities)
    %----------------------------------------------------------------------
    n(m,:) = mode(MDP(m).n.*(MDP(m).n > 0),1);
    
end

% reorder list of model indices for each update
%--------------------------------------------------------------------------
n     = mode(n,1);
for t = 1:T
    if n(t) > 0
        M(t,:) = circshift((1:size(MDP,1)),[0 (1 - n(t))]);
    else
        M(t,:) = 1;
    end
end


return

function MDP = spm_MDP_update(MDP,OUT)
% FORMAT MDP = spm_MDP_update(MDP,OUT)
% moves Dirichlet parameters from OUT to MDP
% MDP - structure array (new)
% OUT - structure array (old)
%__________________________________________________________________________

% check for concentration parameters at this level
%--------------------------------------------------------------------------
try,  MDP.a = OUT.a; end
try,  MDP.b = OUT.b; end
try,  MDP.c = OUT.c; end
try,  MDP.d = OUT.d; end
try,  MDP.e = OUT.e; end

% check for concentration parameters at nested levels
%--------------------------------------------------------------------------
try,  MDP.MDP(1).a = OUT.mdp(end).a; end
try,  MDP.MDP(1).b = OUT.mdp(end).b; end
try,  MDP.MDP(1).c = OUT.mdp(end).c; end
try,  MDP.MDP(1).d = OUT.mdp(end).d; end
try,  MDP.MDP(1).e = OUT.mdp(end).e; end

return


function L = spm_MDP_VB_VOX(MDP,L,t)
% FORMAT L = spm_MDP_VB_VOX(MDP,L,t)
% returns likelihoods from voice recognition (and articulates responses)
% MDP - structure array
% L   - predictive prior over outcomes
% t   - current trial
%
% L   - likelihood of lexical and prosody outcomes
%
% this subroutine determines who is currently generating auditory output
% and produces synthetic speech - or uses the current audio recorder object
% to evaluate the likelihood of the next word
%__________________________________________________________________________


% check for VOX structure
%--------------------------------------------------------------------------
global VOX
global TRAIN
if ~isstruct(VOX), load VOX; VOX.RAND = 0; end
if isempty(TRAIN), TRAIN = 0;              end
if t == 1,         pause(1);               end


if ~isfield(VOX,'msg')
    
    % prepare useful fields in (global) VOX structure
    %----------------------------------------------------------------------
    Data    = imread('recording','png');
    VOX.msg = msgbox('Recording','','custom',Data);
    set(VOX.msg,'Visible','off'), drawnow
    
    % indices of words in lexicon and inferred prosody states
    %----------------------------------------------------------------------
    VOX.io  = spm_voice_i(MDP.label.outcome{1});
    VOX.ip  = find(ismember({VOX.PRO.str},{'amp','dur','Tf','p0','p1','p2'}));
    
    % check for audio recorder
    %----------------------------------------------------------------------
    if ~isfield(VOX,'audio')
        VOX.audio  = audiorecorder(22050,16,1);
    end
    
end

if MDP.VOX == 0 || MDP.VOX == 1
    
    % Agent: computer
    %----------------------------------------------------------------------
    str = MDP.label.outcome{1}(MDP.o(1,1:t));
    fprintf('%i: %s\n',MDP.VOX, str{t});
    i   = ismember(str,' ');
    str = str(~i);
    eof = sum(i);
    
    % if this is the end of a sentence
    %----------------------------------------------------------------------
    if eof == 1 || (~eof && t == MDP.T)
        
        % get lexical and prosody
        %------------------------------------------------------------------
        lexical = spm_voice_i(str);
        prosody = VOX.prosody(:,lexical);
        if MDP.VOX == 0
            speaker = [12;12];
        else
            speaker = [3; 3 ];
        end
        
        % add prosody and articulate
        %------------------------------------------------------------------
        prosody(VOX.ip,:) = MDP.o(2:end,1:numel(str));
        prosody(1,:)      = min(8,prosody(1,:) + 2);
        spm_voice_speak(lexical,prosody,speaker);
        
        
        % TRAIN: prompt for prosody
        %------------------------------------------------------------------
        if TRAIN
            VOX.mute  = 0;
            VOX.depth = 1;
            [i,P]     = spm_voice_i(str);
            prosody   = [];
            while size(prosody,2) ~= numel(str)
                clc
                disp('Please repeat:'), disp(str)
                [SEG,W,prosody] = spm_voice_read(VOX.audio,P);
            end
            clc, disp('Thank you')
            
            % prompt for prosody
            %--------------------------------------------------------------
            for i = 1:numel(VOX.ip)
                for j = 1:size(P,2)
                    L{i + 1,j} = sparse(prosody(VOX.ip(i),j),1,1,8,1);
                end
            end
            
            % uniform priors for spce  (' ')
            %--------------------------------------------------------------
            L{i + 1,t} = ones(8,1)/8;
            
        end
        
        
    end
    
    
elseif MDP.VOX == 2
    
    % user
    %----------------------------------------------------------------------
    if t == 1
        VOX.IT = 1;
        stop(VOX.audio)
        record(VOX.audio,8);
        set(VOX.msg,'Visible','on')
        pause(1);
        set(VOX.msg,'Visible','off')
        
        % toggle to see spectral envelope
        %------------------------------------------------------------------
        VOX.onsets = 0;
        
    end
    
    % get prior over outcomes and synchronise with Lexicon
    %----------------------------------------------------------------------
    io  = VOX.io;                            % indices words in lexicon
    ip  = VOX.ip;                            % indices of prosody
    no  = numel(io);                         % number of outcomes
    nw  = numel(VOX.LEX);                    % number of words in lexicon
    nk  = size(L,2) - t + 1;                 % number of predictions
    P   = zeros(nw,nk);                      % prior over lexicon
    for k = 1:nk
        for i = 1:no
            j = io(i);
            if j
                P(j,k) = L{1,t + k - 1}(i);
            end
        end
    end
    
    % deep segmentation: check for last word (P(:,2) = 0)
    %----------------------------------------------------------------------
    VOX.LL = -128;
    VOX.LW = 0;
    if size(P,2) > 1
        if any(P(:,2))
            if any(P(:,2) < (1 - 1/8))
                VOX.LL = 4;                  % there may be no next word
                VOX.LW = 0;
            else
                VOX.LL = -128;               % there is a subsequent word
                VOX.LW = 0;                  % and this is not the last
            end
        else
            P      = P(:,1);
            VOX.LW = 1;                      % this is the last word
        end
    end
    
    % or direct segmentation (comment out to suppress)
    %----------------------------------------------------------------------
    P  = P(:,1);
    
    % get likelihood of discernible words
    %----------------------------------------------------------------------
    P      = P > 1/128;
    if any(P(:,1))
        
        % log likelihoods
        %------------------------------------------------------------------
        O  = spm_voice_get_word(VOX.audio,bsxfun(@rdivide,P,sum(P)));
        
        % check for end of sentence
        %------------------------------------------------------------------
        try
            A = spm_softmax(O{2}(:,1));    % P(lowest amplitude)
        catch
            A = 1;
        end
        if isempty(O) || A(1) > 1/2
            
            % end of sentence or indiscernible word
            %--------------------------------------------------------------
            L{1,t}( ~io) = 1;
            L{1,t}(~~io) = 0;
            L{1,t}       = L{1,t}/sum(L{1,t});
            
            % prosody likelihoods
            %--------------------------------------------------------------
            for g = 2:numel(L)
                L{g,t} = spm_softmax(spm_zeros(L{g,t}));
            end
            
        else
            
            % lexical likelihoods
            %--------------------------------------------------------------
            LL    = zeros(no,1);
            for i = 1:no
                j = io(i);
                if j
                    LL(i) = O{1}(j);
                end
            end
            L{1,t}  = spm_softmax(LL);
            
            % prosody likelihoods
            %--------------------------------------------------------------
            for g = 2:numel(L)
                L{g,t} = spm_softmax(O{2}(:,ip(g - 1)));
            end
            
        end
        
    end % discernible words
    
    % display word
    %----------------------------------------------------------------------
    [d,w]  = max(L{1,t});
    fprintf('%i: %s\n',MDP.VOX, MDP.label.outcome{1}{w});
    
    
end


