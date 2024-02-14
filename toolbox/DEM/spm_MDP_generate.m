function [MDP] = spm_MDP_generate(MDP)
% active inference and learning using belief propagation (factorised)
% FORMAT [MDP] = spm_MDP_generate(MDP)
%
% Input; MDP(m,n)       - structure array of m models over n epochs
% MDP.U(1, Nf)          - controllable factors
% MDP.V(Np,Nf)          - combinations of latent factors (policies)
% MDP.T                 - number of time steps
%
% MDP.A{Ng}(No(g),Ns(1),...,Ns(Nf)) - likelihood of outcomes, given hidden states
% MDP.B{Nf}(Ns(f),Ns(f),Nu(f))      - transitions among states under control states
% MDP.C{Ng}(No(g),1)    - prior probabilities over outcomes       (Dirichlet counts)
% MDP.D{Nf}(Ns(f),1)    - prior probabilities over initial states (Dirichlet counts)
% MDP.E{Nf}(Nu(f),1)    - prior probabilities over paths          (Dirichlet counts)
% MDP.H{Nf}(Ns(f),1)    - prior probabilities over final states   (Dirichlet counts)
%
% OPTIONS.A             - switch to evaluate explicit action
%
% Outputs:
% MDP.s(Nf,T)           - states   - for each hidden factor
% MDP.o(Ng,T)           - outcomes - for each outcome modality
% MDP.O{Ng,T}           - likelihoods   - for each outcome modality
% MDP.u(Nu,T)           - controls - for each hidden factor
%
% This routine provides solutions of active inference (minimisation of
% variational free energy) using a generative model based upon a Markov
% decision process. The model and inference scheme is formulated in
% discrete space and time. This means that the generative model (and
% process) are hidden Markov models whose dynamics are given by transition
% probabilities among states and the likelihood corresponds to a particular
% outcome conditioned upon hidden states.
%
%
% See also: spm_MDP_VB_XXX, which is the corresponding variational message
% passing scheme.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_XXX.m 8460 2023-12-07 14:01:28Z karl $


% set up and preliminaries
%==========================================================================

% check MDP specification
%--------------------------------------------------------------------------
MDP   = spm_MDP_checkX(MDP);

% initialise model-specific parameters
%==========================================================================
T     = MDP(1).T;                              % number of updates
Nm    = numel(MDP);
for m = 1:Nm

    % number of outcomes and latent states
    %----------------------------------------------------------------------
    Ng(m) = numel(MDP(m).A);                   % number of modalities
    Nf(m) = numel(MDP(m).B);                   % number of factors
    for g = 1:Ng(m)
        No(m,g) = size(MDP(m).A{g},1);         % number of outcomes
    end
    for f = 1:Nf(m)
        Ns(m,f) = size(MDP(m).B{f},1);         % number of hidden states
        Nu(m,f) = size(MDP(m).B{f},3);         % number of hidden paths
    end

end

% parameters of generative model and policies
%==========================================================================
O     = cell(Nm,max(Ng),T);                   % outcomes
for m = 1:Nm

    % likelihood model (for a partially observed MDP)
    %----------------------------------------------------------------------
    for g = 1:Ng(m)

        % normalised likelihood
        %------------------------------------------------------------------
        A{m,g} = spm_norm(MDP(m).A{g});

    end

    % transition probabilities (priors)
    %----------------------------------------------------------------------
    for f = 1:Nf(m)

        % normalised transition probabilities
        %------------------------------------------------------------------
        B{m,f} = spm_norm(MDP(m).B{f});

        % priors over initial hidden states
        %------------------------------------------------------------------
        D{m,f} = spm_norm(MDP(m).D{f});

        % priors over paths (control states)
        %------------------------------------------------------------------
        E{m,f} = spm_norm(MDP(m).E{f});

    end


    % controllable factors
    %======================================================================

    % policies (V) - latent factors
    %----------------------------------------------------------------------
    U{m}       = any(MDP(m).U,1);
    k          = find(U{m});
    u          = spm_combinations(Nu(m,k));
    V{m}       = zeros(size(u,1),Nf(m));
    V{m}(:,k)  = u;
    Np(m)      = size(V{m},1);                 % number of latent policies


    % if states have not been specified, set to 0
    %----------------------------------------------------------------------
    k        = zeros(Nf(m),T);
    try
        i    = find(MDP(m).s);
        k(i) = MDP(m).s(i);
    end
    MDP(m).s = k;

    % if paths have not been specified, set to 0
    %----------------------------------------------------------------------
    k        = zeros(Nf(m),T);
    try
        i    = find(MDP(m).u);
        k(i) = MDP(m).u(i);
    end
    MDP(m).u = k;

    % set outcomes to 0
    %======================================================================
    k        = zeros(Ng(m),T);
    MDP(m).o = k;

    % domains (id)
    %----------------------------------------------------------------------
    id{m}     = MDP(m).id;

end % end model (m)

% inital policy and transitions
%--------------------------------------------------------------------------
try
    tau = MDP(m).tau;
catch
    tau = 2;
end
K     = 1;
PK    = (1 - 1/tau)*eye(Np(m),Np(m)) + (1/tau)/(Np(m));
PK    = spm_norm(PK);

% updating over successive time points
%==========================================================================
for t = 1:T

    % generate hidden states and controls for each agent or model
    %======================================================================
    for m = 1:Nm

        % initialise and propagate control state (path)
        %==================================================================
        for f = 1:Nf(m)
            if ~MDP(m).u(f,t)
                if t > 1

                    % previous path
                    %------------------------------------------------------
                    MDP(m).u(f,t) = MDP(m).u(f,t - 1);

                else

                    % otherwise sample a path
                    %------------------------------------------------------
                    pu            = spm_norm(E{m,f});
                    MDP(m).u(f,t) = spm_sample(pu);

                end
            end
        end

        % action generating outcomes (for controllable factors)
        %==================================================================
        if t > 1

            % implicit action
            %--------------------------------------------------------------
            if isfield(id{m},'hid')

                % prior over policy
                %----------------------------------------------------------
                K   = spm_sample(spm_softmax(G));

            else

                % smooth policy transitions
                %----------------------------------------------------------
                K   = spm_sample(PK(:,K));

            end

            for f = 1:Nf(m)
                if U{m}(f)

                    % selected action
                    %------------------------------------------------------
                    MDP(m).u(f,t - 1) = V{m}(K,f);

                end
            end

        end % end generating paths or actions


        % sample state if not specified
        %==================================================================
        for f = 1:Nf(m)
            if ~MDP(m).s(f,t)

                if t > 1

                    % the next state is generated by state transititions
                    %----------------------------------------------------------
                    ps = B{m,f}(:,MDP(m).s(f,t - 1),MDP(m).u(f,t - 1));

                else

                    % unless it is the initial state
                    %------------------------------------------------------
                    ps = spm_norm(D{m,f});

                end

                MDP(m).s(f,t) = spm_sample(ps);

            end

        end % end generating states

    end % end generating over models


    % generate outcomes O{m,g,t} for each agent or model
    %======================================================================
    for m = 1:Nm

        % generate outcomes
        %------------------------------------------------------------------
        for g = 1:Ng(m)

            % domain and codomain of A{g}
            %--------------------------------------------------------------
            [j,i] = spm_get_edges(id{m},g,MDP(m).s(:,t));
            for o = i

                % sample from likelihood, given hidden state
                %==========================================================
                ind           = num2cell(MDP(m).s(j,t));
                O{m,o,t}      = MDP(m).A{g}(:,ind{:});
                MDP(m).o(o,t) = spm_sample(O{m,o,t});

            end
        end
    end


    % Updating of controls (via G)
    %======================================================================
    for m = 1:Nm

        % prior over current path or control state (Pu)
        %==================================================================
        G  = zeros(Np(m),1);
        if isfield(id{m},'hid')
            for f = 1:Nf(m)
                for k = 1:Np(m)

                    if V{m}(k,f)

                        % transitions for this policy
                        %--------------------------------------------------
                        BP{m,f,k} = B{m,f}(:,:,V{m}(k,f));

                    else

                        % transitions and novelty for current path
                        %--------------------------------------------------
                        BP{m,f,k} = B{m,f}(:,:,MDP(m).u(f,t));

                    end
                end
            end

            % hidden states (Q) and induction
            %==============================================================
            for f = 1:Nf(m)
                Q{m,f} = sparse(MDP(m).s(f,t),1,1,Ns(m,f),1);
            end
            [R,r] = spm_induction(BP(m,:,:),Q,(T - t),id{m});

            % Expected free energy of subsequent action
            %==============================================================
            if any(R,'all')
                for k = 1:Np(m)

                    % predictive posterior
                    %------------------------------------------------------
                    for f = r
                        P{f,k} = BP{m,f,k}(:,MDP(m).s(f,t));
                    end

                    % inductive constraints over states
                    %------------------------------------------------------
                    G(k) = R*P{r,k};
                end
            end
        end

    end % end of loop over models (agents)

end % end of loop over time


% loop over models and prepare outputs
%==========================================================================
for m = 1:size(MDP,1)

    % assemble results and place in NDP structure
    %======================================================================
    MDP(m).T  = T;            % number of outcomes
    MDP(m).O  = O(m,:,:);     % outcomes
    MDP(m).O  = shiftdim(MDP(m).O,1);

end % end loop over models (m)

return


% auxillary functions
%==========================================================================

function i  = spm_sample(P)
% log of numeric array plus a small constant
%--------------------------------------------------------------------------
if islogical(P)
    i = find(P,1);
else
    i = find(rand < cumsum(P),1);
end


function A  = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
if isnumeric(A)
    A           = rdivide(A,sum(A,1));
    A(isnan(A)) = 1/size(A,1);
end


function [R,hif] = spm_induction(B,Q,N,id)
% Inductive inference about next state
% FORMAT [R,hif] = spm_induction(B,Q,N,id)
%--------------------------------------------------------------------------
% A{1,g}   - likelihood mappings from hidden states
% B{1,f,k} - belief propagators (policy-dependent probability transitions)
% Q{1,f}   - cell array of posteriors over states
%
% N      - induction depth
% id     - domain structure
%
%   id.hid(Nf,Ni)  - indices of Ni intended  states
%   id.cid(Nf,Ni)  - indices of Ni suprising states
%
% R      - tensor encoding unconstrained states over hif factors
% hif    - factors of tensor
%
% This subroutine returns constraints on the next state based upon
% backwards induction of a simple sort; i.e., using backwards propagators
% to identify paths of least action using logical operators.
%
% In addition, constraints can be specified in a small number of latent
% factors by supplying a matrix of constraints, where each column
% corresponds to a distinct constraint and the number of rows corresponds
% to the number of hidden factors. This constraint matrix contains the
% index of the costly state for each factor. If an index is zero, the
% constraint is taken to be independent of the corresponding factor;
% otherwise, this conditional factor has to have a high posterior over the
% next state, before the constraint is implemented.
%__________________________________________________________________________

% Preliminary checks (for no priors over end states)
%==========================================================================

% Convert marginals to indices and find factors
%--------------------------------------------------------------------------
if isfield(id,'hid')
    hid   = id.hid;                           % intended states
    hif   = find(any(hid,2))';                % in hif factors
else
    hid   = [];
    hif   = [];
end

% Deal with constraints
%--------------------------------------------------------------------------
if isfield(id,'cid')
    cid   = id.cid;                           % contrained factors
    nid   = cid;                              % conditioning factors
    hif   = find(all(cid,2))';                % in hif factors
    nid(hif,:) = 0;

    % size of contrained factors
    %----------------------------------------------------------------------
    for f = hif
        Ns(f) = size(B{f},1);
    end
    Ns    = [Ns,1];

    % constraint tensor over hid factors
    %----------------------------------------------------------------------
    D     = true(Ns);                        % unconstrained states
    for i = 1:size(cid,2)

        % posterior of constraint violation
        %------------------------------------------------------------------
        q     = 1;
        for f = find(nid(:,i))'
            q = q*Q{f}(cid(f,i));
        end
        if q > (1 - 1/8)
            ind  = num2cell(cid(hif,i));
            j    = sub2ind(Ns,ind{:});
            D(j) = false;
        end
    end
else
    D = true;
end

% Return if there are no intended states or constraints
%--------------------------------------------------------------------------
if isempty(hif), R = false; return, end
if isempty(hid), R = 32*D;  return, end

% Threshold transition probabilities
%--------------------------------------------------------------------------
u     = 1/16;                            % probability threshold
b     = cell(1,numel(hif));
for f = hif
    b{f}  = false;
    for k = 1:size(B,3)
        b{f} = b{f} | gt(B{1,f,k}, max(B{1,f,k})*u);
    end
end

% Kronecker tensor products (sparse)
%--------------------------------------------------------------------------
Bf    = 1;
Qf    = 1;
for f = hif
    Ns(f) = size(B{f},1);                % numer of states for hif
    Bf    = spm_kron(b{f},Bf);           % unconstrained transitions
    Qf    = spm_kron(Q{f},Qf);           % posterior over states
end

% Expected cost in latent state space
%==========================================================================
Bf    = and(Bf,D(:));                    % constrained transitions


% Backwards induction: from end states
%==========================================================================

% hid are indices (of multiple endpoints)
%--------------------------------------------------------------------------
for i = 1:size(hid,2)
    for f = hif
        h{f} = false(Ns(f),1);
        h{f}(hid(f,i)) = true;
    end
    I     = true;
    for f = hif
        I = spm_kron(h{f},I);
    end
    Pf(:,i) = logical(I);
end


% Backwards induction: paths of least action
%==========================================================================
for i = 1:size(Pf,2)

    % for this end state
    %----------------------------------------------------------------------
    I = logical(Pf(:,i));

    % backwards protocol (for paths with a well-defined end state)
    %----------------------------------------------------------------------
    for n = 1:min(N,16)

        % any preceding states % & that have not been previously occupied
        %------------------------------------------------------------------
        I(:,n + 1) = any(Bf(I(:,n),:),1)'; % & ~any(I,2);
    end

    % Find most likely point on paths of least action
    %----------------------------------------------------------------------
    G(:,i) = I'*Qf;
    P{i}   = I;

end

% precise log prior over next state
%==========================================================================
G(1,:) = 0;                        % preclude current states
[d,n]  = max(G,[],1);              % next intended state
i      = d > u;                    % provided it exists
if any(i)

    % eliminate inaccessible end states
    %----------------------------------------------------------------------
    P     = P(i);
    n     = n(i);
    [n,i] = min(n);                % to the i-th end state

    % precise log prior over next state
    %----------------------------------------------------------------------
    P     = P{i}(:,max(n - 1,1));
    R     = single(reshape(full(P),[Ns 1]));
    R     = shiftdim(32*R,-1);
else
    R = false;
end

return




