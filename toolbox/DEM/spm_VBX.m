function [P,F] = spm_VBX(O,P,A,id)
% variational Bayes estimate of categorical posterior over factors
% FORMAT [Q,F] = spm_VBX(O,P,A,id)
%
% O{g}    -  outcome probabilities over each of G modalities
% P{f}    -  (empirical) prior over each of F factors
% A{g}    -  likelihood tensor for modality g
%
% Q{f}    -  variational posterior for each of F factors
% F       -  (-ve)  variational free energy or ELBO
%
% This routine is a simple implementation of variational Bayes for discrete
% state space models under a mean field approximation, in which latent
% states are partitioned into factors (and the distribution over outcomes
% is also assumed to be conditionally independent). It takes cell arrays of
% outcome probabilities, prior probabilities over factors and a likelihood
% tensor parameterising the likelihood of an outcome for any combination of
% latent states. This version uses belief propagation, read as:
%
% Propagation:  a non-iterative heuristic but numerically accurate scheme
% that replaces the variational density over hidden factors with the
% marginal over the exact posterior
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if isfield(id,'i')                    % use a subset of outcomes
    ig = id.g(id.i);
else
    ig = id.g;                        % use complete partition
end
Nf  = numel(P);                       % number of factors
F   = 0;                              % ELBO


% belief propagation with marginals of exact posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intial update of domain factors Q{ff} with children
%==========================================================================
if isfield(id,'ff')
    for p = 1:numel(ig)               % for each partition of modalities
        for g = ig{p}                 % for each modality in partition

            % Get parents and children of this modality
            %--------------------------------------------------------------
            [j,i] = spm_parents(id,g,P);

            % update parents iff (sic) they are domain factors
            %--------------------------------------------------------------
            if any(ismember(j,id.ff))
                [P,f] = spm_VBX_update(P,A{g},O,i,j);
                F     = F + f;
            end
        end
    end
end


% Update domain factors Q{ff} without children
%==========================================================================
if isfield(id,'ff')

    % reduced prior R{f} with plausible (sparse) states s{f} 
    %----------------------------------------------------------------------
    R     = cell(1,Nf);
    s     = cell(1,Nf);
    for f = 1:Nf
        s{f}  = find(P{f} > exp(-8)); % reduced states
        R{f}  = P{f}(s{f});           % reduced prior
    end

    % cycle over partition subsets
    %----------------------------------------------------------------------
    for p = 1:numel(ig)

        % domain factors for this partition
        %------------------------------------------------------------------
        if iscell(id.ff)
            ff = id.ff{p};
        else
            ff = id.ff;
        end

        % size of domain states
        %------------------------------------------------------------------
        Nff   = numel(ff);
        Ns    = ones(1,Nff);
        for f = 1:Nff
            Ns(f) = numel(s{ff(f)});
        end

        % if there is uncertainity over combinations
        %------------------------------------------------------------------
        if prod(Ns) > 1

            % accumulate log-likelihood of domain factors
            %--------------------------------------------------------------
            L     = 0;
            for g = ig{p}
                L = L + spm_VBX_update_L(A,O,P,R,s,g,ff,id);
            end

            % exact posterior and marginals
            %--------------------------------------------------------------
            U     = spm_times(exp(L),R{ff});
            Z     = sum(U,'all');
            if Z
                Q = spm_margin(U/Z);
            else
                Q = spm_margin(U + 1/numel(U));
            end

            % update marginals
            %--------------------------------------------------------------
            for i = 1:numel(Q)
                R{ff(i)}   = Q{i};
            end

            % update full posterior
            %--------------------------------------------------------------
            for f = ff
                P{f}(:)    = 0;
                P{f}(s{f}) = R{f};
            end

        end % numel(L)

    end % partition of modalities

end % state-dependent domains and co-domains


% accumulate likelihoods (L) over (partition) of modalities
%==========================================================================

% number of plausible combinations (Nq) of domain states
%--------------------------------------------------------------------------
Nq    = numel(spm_edges(id,1,P));
F     = repmat(F,Nq,1);
S     = repmat(P,Nq,1);
for p = 1:numel(ig)
    for g = ig{p}

        % Get parents and children of this modality
        %------------------------------------------------------------------
        [jq,iq] = spm_edges(id,g,P);

        % for each combination of domain states
        %------------------------------------------------------------------
        for q = 1:Nq
            i      = iq{q};
            j      = jq{q};
            [Q,Fq] = spm_VBX_update(S(q,:),A{g},O,i,j);
            S(q,j) = Q(j);
            F(q)   = F(q) + Fq;
        end

    end % modalities

end % partition of modalities

% Bayesian model average over (Nq) combinations
%==========================================================================

% full posterior and average free energy
%--------------------------------------------------------------------------
p     = spm_softmax(F(:));
F     = F'*p;
for f = 1:Nf
    P{f}(:) = 0;
    for q = 1:Nq
        P{f} = P{f} + S{q,f}*p(q);
    end
end

return

% subfunctions
%==========================================================================

function [P,F] = spm_VBX_update(P,A,O,i,j)
% Update (sparse) posteriors
% FORMAT [P,F] = spm_VBX_update(P,A,O,i,j)
% P  - marginal probabilities
% A  - likelihood tensors or function handles
% O  - outcome probabilities
% i  - children of likelihood mapping (modalities)
% j  - parents  of likelihood mapping (factors)
%
% This subroutine performs exact Bayesian inference or belief updating,
% given current priors and outcomes. It uses sparse tensor operations by
% removing near zero probabilities from P{f} and recording that the nonzero
% probabilities in R{f} and the accompanying indices, s{f}. The ELBO
% (variational free energy) is computed explicitly as a partition function.
%
% This routine will test for likelihood mappings specified with function
% handles. If the likelihood has a functional form, there is a test for
% joint likelihoods (numeric arrays) or marginal likelihoods (cell arrays).
% This allows one to specify the marginal likelihood given an outcome, as
% opposed to a tensor over joint states. This is appropriate when the
% implicit likelihood tensor can be marginalised. A special case of this is
% when a sub- tensor can be marginalised. These sub-tensor are specified by
% a restricted number of dimensions, returned in the vector r.
%__________________________________________________________________________

% sparse form
%--------------------------------------------------------------------------
Nf    = numel(j);
[R,s] = spm_VBX_sparse(P(j));

% likelihoods
%--------------------------------------------------------------------------
M     = num2cell(zeros(1,Nf));                % marginal log-likelihoods
L     = 0;                                    % joint log-likelihood
for o = i
    if isa(A,'function_handle')

        % Functional likelihood
        %--------------------------------------------------------------
        if nargout(A) > 1
            [U,r] = A(O{o},s,P(j));
        else
            U     = A(O{o},s);
            r     = 1:numel(j);
        end

        if iscell(U)

            % marginal log-likelihoods (cell array)
            %--------------------------------------------------------------
            for f = 1:numel(U)
                M{f} = M{f} + spm_log(U{f});
            end

        else

            % joint log-likelihood (numeric tensor)
            %--------------------------------------------------------------
            L = L + spm_log(U);
        end

    else

        % joint log-likelihood from tensor
        %------------------------------------------------------------------
        U = spm_dot(A(:,s{:}),O{o});
        L = L + spm_log(U);

    end
end

% Posterior marginals
%==========================================================================
if iscell(U)

    % loop over marginals
    %----------------------------------------------------------------------
    F     = 0;
    for f = 1:numel(U)
        U{f}  = exp(M{f});                 % exponentiate
        U{f}  = U{f}.*R{r(f)};             % posterior unnormalised
        Z     = sum(U{f},'all');           % partition coefficient
        if Z
            F    = F + spm_log(Z);         % negative free energy
            Q{f} = U{f}/Z;                 % marginal  posteriors
        else
            F    = -32;
            Q{f} = U{f} + 1/numel(U{f});
        end
    end

    % update posteriors
    %----------------------------------------------------------------------
    j     = j(r);
    for f = 1:numel(j)
        P{j(f)}(:)       = 0;
        P{j(f)}(s{r(f)}) = Q{f};
    end

else

    % if there is any posterior uncertainty
    %----------------------------------------------------------------------
    if numel(L) > 1

        % Posterior marginals
        %------------------------------------------------------------------
        L     = exp(L);                    % exponentiate
        U     = spm_times(L,R{:});         % posterior unnormalised
        Z     = sum(U,'all');              % partition coefficient
        if Z
            F = spm_log(Z);                % negative free energy
            Q = spm_margin(U/Z);           % marginal  posteriors
        else
            F = -32;
            Q = spm_margin(U + 1/numel(U));
        end

        % update marginals
        %------------------------------------------------------------------
        for f = 1:min(numel(j),numel(Q))
            P{j(f)}(:)    = 0;
            P{j(f)}(s{f}) = Q{f};
        end

    elseif numel(L)        
        F = L;
    else
        F = -32;
    end

end

return

function  L = spm_VBX_update_L(A,O,P,R,s,g,ff,id)
% FORMAT: L = spm_VBX_update_L(A,O,P,R,s,g,ff,id)
% A  - likelihood tensor or function handles
% O  - cell array of outcomes
% P  - tensor of joint posteriors over states
% R  - reduced or sparse joint posterior
% s  - sparsity indices
% g  - modality to accumulate
% ff - the main factors
% id - structure of indices or identifiers
%
% This subroutine evaluates the log likelihood (L) of all combinations of
% domain states (ff) for modality g. The domain states can either be
% precomputed and stored in a cell array (with the same size as the domain
% states), were each cell contains the indices of the parents associated
% with the domain state. Alternatively, the domain states can be returned
% by a function of the modality (G) and the indices of the domain states
% (ind). See spm_parents.m
%__________________________________________________________________________

% size of domain states
%--------------------------------------------------------------------------
Nff   = numel(ff);
Ns    = ones(1,Nff);
for f = 1:Nff
    Ns(f) = numel(s{ff(f)});
end

% for all combinations of domain states
%==========================================================================
L     = zeros(Ns);
c     = spm_combinations(Ns);
for i = 1:size(c,1)

    % parents (f) and children (j), given states (ind)
    %----------------------------------------------------------------------
    ind   = cell(1,Nff);
    for j = 1:Nff
        ind{j} = s{ff(j)}(c(i,j));
    end

    % parents
    %----------------------------------------------------------------------
    if isfield(id,'fg')
        if iscell(id.fg)
            f = id.fg{g}{ind{:}};
        else
            f = id.fg(g,[ind{:}]);
        end
    else
        f = id.A{g};
    end

    % children
    %----------------------------------------------------------------------
    if isfield(id,'gg')
        if iscell(id.gg)
            j = id.gg{g}{ind{:}};
        else
            j = id.gg(g,[ind{:}]);
        end
    else
        j = g;
    end

    % predicted outcomes
    %----------------------------------------------------------------------
    if isa(A{g},'function_handle')
        qo = A{g}(P(f));
    else
        qo = spm_dot(A{g}(:,s{f}),[O{g},R(f)],1);
    end

    % likelihood for these domains and codomains
    %----------------------------------------------------------------------
    for o = j
        L(i) = L(i) + spm_log(spm_dot(qo,O{o}));
    end

end

return

function [R,s] = spm_VBX_sparse(P)
% subsamples non zero marginals
%__________________________________________________________________________
Nf    = numel(P);
R     = cell(1,Nf);
s     = cell(1,Nf);
for f = 1:Nf
    s{f}  = P{f} > exp(-16);           % reduced states
    R{f}  = P{f}(s{f});                % reduced prior
end

return

function A  = spm_log(A)
% log of numeric array plus a small constant
%--------------------------------------------------------------------------
if islogical(A)
    A = -32*(~A);
else
    A = max(log(A),-32);
end

function [Y] = spm_margin(X)
% Marginal densities over a multidimensional array of probabilities
% FORMAT [Y] = spm_margin(X)
% X  - numeric array of probabilities
% Y  - cell array of marginals
%
% See also: spm_dot, spm_marginal (with vectorization)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging

% evaluate marginals
%--------------------------------------------------------------------------
n     = ndims(X);
Y     = cell(n,1);
for i = 1:n
    j    = 1:n;
    j(i) = [];
    Y{i} = reshape(spm_sum(X,j),[],1);
end


return

function [X] = spm_times(X,varargin)
% Outer product
% FORMAT [X] = spm_times(X,x,...)
% X  - numeric array of probabilities
% x  - cell array of marginals
%
% X  = X.*spm_cross(x,...)
%
% See also: spm_dot, spm_cross
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging

% evaluate outer product
%--------------------------------------------------------------------------
x     = varargin;
n     = numel(x);
for i = 1:n
    j    = ones(1,n);
    j(i) = numel(x{i});
    X    = times(X,reshape(full(x{i}),[j,1]));
end


