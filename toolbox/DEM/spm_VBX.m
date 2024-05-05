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
%
% see: spm_MDP_VB_XXX.m (NOTES)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if isfield(id,'ig')                          % use a subset of outcomes
    i  = id.ig(end);
    ig = id.g(i);
else
    ig = id.g;                               % use complete partition
end

% belief propagation with marginals of exact posterior
%==========================================================================

% prior : s{f} plausible states of factor f
%--------------------------------------------------------------------------
Nf    = numel(P);
R     = cell(1,Nf);
for f = 1:Nf
    s{f}  = find(P{f} > exp(-8));     % reduced states
    Ns(f) = numel(s{f});              % number of reduced states
    R{f}  = P{f}(s{f});               % reduced prior
end
F     = 0;                            % ELBO

% intial updpate for domain factors with children
%==========================================================================
if isfield(id,'ff')
    for p = 1:numel(ig)               % for each partition of modalities
        for g = ig{p}                 % for each modality in partition

            % Get the latent factors (and outcomes) for this modality
            %--------------------------------------------------------------
            [j,i] = spm_get_edges(id,g,P);

            % if the parent of this modality is a domain factor
            %--------------------------------------------------------------
            if any(ismember(j,id.ff))

                % likelihoods
                %----------------------------------------------------------
                L = 0;
                for o = i
                    L = L + spm_log(spm_dot(A{g}(:,s{j}),O{o}));
                end
                L = exp(L);

                % Posterior marginals
                %==========================================================

                % eliminate factors with precise posteriors
                %----------------------------------------------------------
                i = size(L);
                r = find(i > 1);
                L = reshape(L,[i(r) 1 1]);

                % if L is not a vector
                %----------------------------------------------------------
                if numel(j) > 1
                    j = j(r);
                end

                if numel(j)

                    % Posterior marginals
                    %------------------------------------------------------
                    U     = spm_times(L,R{j});   % posterior unnormalised
                    Z     = sum(U,'all');        % partition coefficient
                    if Z
                        F = F + spm_log(Z);      % negative free energy
                        Q = spm_margin(U/Z);     % marginal  posteriors
                    else
                        F = F - 32;
                        Q = spm_margin(U + 1/numel(U));
                    end

                    % Posterior
                    %------------------------------------------------------
                    for i = 1:numel(j)
                        R{j(i)} = Q{i};
                    end

                else
                    F = F + spm_log(L);
                end

            end

        end % modalities

    end % partition of modalities

    % full posterior and average free energy
    %----------------------------------------------------------------------
    for f = 1:Nf
        P{f}(s{f}) = R{f};                      % posterior
        s{f}       = find(P{f} > exp(-8));      % reduced states
        Ns(f)      = numel(s{f});               % number of reduced states
        R{f}       = P{f}(s{f});                % reduced prior
    end
end

% Update domain factors Q{ff} using exact inference
%==========================================================================
if isfield(id,'ff')

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
        Nff    = numel(ff);

        % likelihood of domain factors 
        %------------------------------------------------------------------
        L     = zeros(Ns(ff));
        ind   = cell(1,Nff);
        for g = ig{p}
            for i = 1:numel(L)

                % parents (f) and children (j), given states (ind)
                %----------------------------------------------------------
                sf    = spm_index(Ns(ff),i);
                for k = 1:Nff
                   ind{k} = s{ff(k)}(sf(k));
                end

                % parents
                %----------------------------------------------------------
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
                %----------------------------------------------------------
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
                %----------------------------------------------------------
                qo    = spm_dot(A{g}(:,s{f}),[O(g) R(f)],1);

                % likelihood for these domains and codomains
                %----------------------------------------------------------
                for o = j
                    L(i) = L(i) + spm_log(spm_dot(qo,O{o}));
                end
            end
        end

        % exact posterior and marginals
        %------------------------------------------------------------------
        U     = spm_times(exp(L),R{ff});
        Z     = sum(U,'all');
        if Z
            Q = spm_margin(U/Z);
        else
            Q = spm_margin(U + 1/numel(U));
        end

        % update marginals
        %------------------------------------------------------------------
        for i = 1:numel(Q)
            R{ff(i)} = Q{i};
        end

    end % partition of modalities

    % update full posterior
    %----------------------------------------------------------------------
    for f = 1:Nf
        P{f}(:)    = 0;
        P{f}(s{f}) = R{f};
    end

end % end state-dependent domains and co-domains



% accumulate likelihoods (L) over partition of modalities
%==========================================================================

% number of plausible domain states
%--------------------------------------------------------------------------
Nq = numel(spm_edges(id,1,P));
F  = repmat(F,1,Nq);
R  = repmat(R,Nq,1);

for p = 1:numel(ig)
    for g = ig{p}

        % Get the latent factors (and outcomes) for this modality
        %------------------------------------------------------------------
        [jq,iq] = spm_edges(id,g,P);

        for q = 1:Nq

            % indices for this domain state
            %--------------------------------------------------------------
            j = jq{q};
            i = iq{q};

            % likelihoods
            %--------------------------------------------------------------
            L = 0;
            for o = i
                L = L + spm_log(spm_dot(A{g}(:,s{j}),O{o}));
            end
            L = exp(L);

            % Posterior marginals
            %==============================================================

            % factors to update (eliminate factors with precise posteriors)
            %--------------------------------------------------------------
            i = size(L);
            r = find(i > 1);
            L = reshape(L,[i(r) 1 1]);

            % if L is not a vector
            %--------------------------------------------------------------
            if numel(j) > 1
                j = j(r);
            end

            if numel(j)

                % Posterior marginals
                %----------------------------------------------------------
                U     = spm_times(L,R{q,j});       % posterior unnormalised
                Z     = sum(U,'all');              % partition coefficient
                if Z
                    F(q) = F(q) + spm_log(Z);      % negative free energy
                    Q    = spm_margin(U/Z);        % marginal  posteriors
                else
                    F(q) = F(q) - 32;
                    Q    = spm_margin(U + 1/numel(U));
                end

                % Posterior
                %----------------------------------------------------------
                for i = 1:numel(j)
                    R{q,j(i)} = Q{i};
                end

            else
                F(q) = F(q) + spm_log(L);
            end

        end % domain models

    end % modalities

end % partition of modalities

% Bayesian model average
%==========================================================================

% full posterior and average free energy
%--------------------------------------------------------------------------
p     = spm_softmax(F(:));
F     = F*p;
for f = 1:Nf
    P{f}(:) = 0;
    for q = 1:Nq
        P{f}(s{f}) = P{f}(s{f}) + R{q,f}*p(q);
    end
end

return


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


