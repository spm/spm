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
% tensor parameterising the likelihood of an outcome for any combination
% of latent states. METHOD switches among number of approximate schemes:
%
% 'propagation' :  a non-iterative heuristic but numerically accurate scheme
% that replaces the variational density over hidden factors with the
% marginal over the exact posterior [default]
%
% 'variational' :  a vanilla variational scheme that uses a coordinate descent
% over a small number (8) iterations (fixed point iteration). This is
% computationally efficient, using only nontrivial posterior probabilities
% and the domain of the likelihood mapping in tensor operations.
%
% 'sparse'  :  as for the exact scheme but suitable for sparse tensors
%
% see: spm_MDP_VB_XXX.m (NOTES)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
METHOD = 'propagation';                      % belief propagation scheme
if isfield(id,'ig')                          % use a subset of outcomes
    i  = id.ig(end);
    ig = id.g(i);
else
    ig = id.g;                               % use complete partition 
end

switch METHOD

    case 'propagation'

        % belief propagation with marginals of exact posterior
        %==================================================================

        % log prior : s{f} plausible states of factor f
        %------------------------------------------------------------------
        Nf    = numel(P);
        for f = 1:Nf
            s{f}  = find(P{f} > exp(-8));
            Ns(f) = numel(s{f});
            R{f}  = P{f}(s{f});
        end

        % Update domain factors Q{ff} using exact inference
        %==================================================================
        if isfield(id,'ff')

            % cycle over partition subsets
            %--------------------------------------------------------------
            for p = 1:numel(ig)

                % domain factors for this partition
                %----------------------------------------------------------
                if iscell(id.ff)
                    ff = id.ff{p};
                else
                    ff = id.ff;
                end

                % plausible domains id.fg{g}
                %----------------------------------------------------------
                for g = ig{p}
                    if isfield(id,'fg')
                        fg{g} = id.fg{g}(s{ff});
                    else
                        fg{g} = id.A{g};
                    end
                    if isfield(id,'gg')
                        gg{g} = id.gg{g}(s{ff});
                    else
                        gg{g} = g;
                    end
                end

                % likelihood of domain factors
                %----------------------------------------------------------
                L     = ones(Ns(ff));
                for g = ig{p}
                    for i = 1:numel(L)
                        if iscell(fg{g})
                            f = fg{g}{i};
                        else
                            f = fg{g};
                        end
                        if iscell(gg{g})
                            j = gg{g}{i};
                        else
                            j = gg{g};
                        end

                        % predicted outcomes
                        %--------------------------------------------------
                        qo    = spm_dot(A{g}(:,s{f}),[O(g) R(f)],1);

                        % likelihood for these domains and codomains
                        %--------------------------------------------------
                        for o = j
                            L(i) = L(i)*spm_dot(qo,O{o});
                        end
                    end
                end

                % exact posterior and marginals
                %----------------------------------------------------------
                U     = spm_times(L,R{ff});
                Z     = sum(U,'all');
                if Z
                    Q = spm_margin(U/Z);
                else
                    Q = spm_margin(U + 1/numel(U));
                end

                % update marginals
                %----------------------------------------------------------
                for i = 1:numel(Q)
                    f          = ff(i);
                    P{f}(:)    = 0;
                    P{f}(s{f}) = Q{i};
                end

            end % partition of modalities

        end % end state-dependent domains and co-domains


        % accumulate likelihoods (L) over partition of modalities
        %==================================================================
        L     = 1;
        for p = 1:numel(ig)
            for g = ig{p}

                % Get the latent factors (and outcomes) for this modality
                %----------------------------------------------------------
                [j,i] = spm_get_edges(id,g,P);

                % deal with degenerate mappings
                %----------------------------------------------------------
                l     = (1:numel(j))';
                [j,k] = unique(j,'stable');
                l(k)  = [];
                k     = [0;k;l] + 1;

                % likelihoods
                %----------------------------------------------------------
                LL    = 1;
                for o = i
                    Ag = permute(A{g},k);
                    LL = LL.*spm_dot(Ag(:,s{j}),O{o});
                end

                % accumulate
                %----------------------------------------------------------
                if numel(j) > 1
                    [j,i] = sort(j);
                    LL    = permute(LL,i);
                end
                k     = ones(1,Nf + 1); k(j) = size(LL,1:numel(j));
                L     = times(L,reshape(LL,k));
            end

        end % partition of modalities

        % factors to update (eliminate factors with precise posteriors)
        %------------------------------------------------------------------
        i = size(L);
        r = find(i > 1);
        L = reshape(L,[i(r) 1 1]);

        if isempty(r), F = spm_log(L); return; end

        % Posterior marginals
        %------------------------------------------------------------------
        U     = spm_times(L,R{r});                 % posterior unnormalised
        Z     = sum(U,'all');                      % partition coefficient
        if Z
            F = spm_log(Z);                        % negative free energy
            Q = spm_margin(U/Z);                   % marginal  posteriors
        else
            F = -32;                               % negative free energy
            Q = spm_margin(U + 1/numel(U));        % marginal  posteriors
        end

        % Posterior
        %------------------------------------------------------------------
        for i = 1:numel(r)
            f          = r(i);
            P{f}(:)    = 0;
            P{f}(s{f}) = Q{i};
        end


    case 'variational'

        % (iterative) variational scheme (Not suitable for partitions)
        %==================================================================

        % log prior : i{f} plausible states of factor f
        %------------------------------------------------------------------
        Nf    = numel(P);
        for f = 1:Nf
            s{f}  = find(P{f} > exp(-8));
            Q{f}  = P{f}(s{f});
            R{f}  = Q{f};
            Ns(f) = numel(s{f});
            LP{f} = spm_vec(spm_log(Q{f}));
        end


        % Update domain factors Q{ff} using exact inference
        %==================================================================
        if isfield(id,'ff')

            % plausible domains id.fg{g}
            %--------------------------------------------------------------
            Ns    = Ns(id.ff);
            for g = ig
                if isfield(id,'fg')
                    fg{g} = id.fg{g}(s{id.ff});
                else
                    fg{g} = id.A{g};
                end
                if isfield(id,'gg')
                    gg{g} = id.gg{g}(s{id.ff});
                else
                    gg{g} = g;
                end
            end
            

            % likelihood of domain factors
            %--------------------------------------------------------------
            L     = ones(Ns);
            for g = ig
                for i = 1:numel(L)
                    if iscell(fg{g})
                        f = fg{g}{i};
                    else
                        f = fg{g};
                    end
                    if iscell(gg{g})
                        j = gg{g}{i};
                    else
                        j = gg{g};
                    end
                    for o = j
                        L(i) = L(i)*spm_dot(A{g}(:,s{f}),{O{o},Q{f}});
                    end
                end
            end

            % exact posterior and marginals
            %--------------------------------------------------------------
            U     = spm_times(L,R{id.ff});
            Z     = sum(U,'all');
            L     = spm_margin(U/Z);
            for i = 1:numel(L)
                f          = id.ff(i);
                Q{f}       = L{i};
                P{f}(s{f}) = Q{f};
            end

        end

        % accumulate log likelihoods over modalities
        %==================================================================
        L     = 0;
        for g = ig
            [j,i] = spm_get_edges(id,g,P);
            j     = unique(j,'stable');
            LL    = 0;
            for o = i
                LL = LL + spm_log(spm_dot(A{g}(:,s{j}),O{o}));
            end
            if numel(j) > 1
                [j,i] = sort(j);
                LL    = permute(LL,i);
            end
            k     = ones(1,Nf + 1); k(j) = size(LL,1:numel(j));
            L     = plus(L,reshape(LL,k));
        end

        % factors to update (eliminate factors with precise posteriors)
        %------------------------------------------------------------------
        i = size(L);
        r = find(i > 1);
        L = reshape(L,[i(r) 1 1]);

        if isempty(r), F = spm_log(L); return; end

        % variational iterations
        %------------------------------------------------------------------
        Z     = -Inf;
        for v = 1:8
            F     = 0;
            for i = 1:numel(r)

                % marginal log likelihood
                %----------------------------------------------------------
                f    = r(i);
                LL   = spm_vec(spm_dot(L,Q(r),i));

                % posterior
                %----------------------------------------------------------
                Q{f} = spm_softmax(LL + LP{f});

                % (-ve) free energy (partition coefficient)
                %----------------------------------------------------------
                F    = F + Q{f}'*(LL + LP{f} - spm_log(Q{f}));

            end

            % convergence
            %--------------------------------------------------------------
            dF = F - Z;
            if dF < 1/128
                break
            else
                Z = F;
            end
        end

        % Posterior
        %------------------------------------------------------------------
        for f = r
            P{f}(:)    = 0;
            P{f}(s{f}) = Q{f};
        end


    case 'sparse'

        % approximation with marginals suitable for sparse tensors
        %==================================================================
        Nf    = size(A{1});
        L     = 1;
        for g = ig
            L = L.*(O{g}'*A{g}(:,:));              % likelihood over modalities
        end
        U     = spm_vec(L).*spm_vec(spm_cross(P)); % posterior unnormalised
        Z     = sum(U,'all');                      % partition coefficient
        F     = spm_log(Z);                        % negative free energy 
        U     = reshape(U/Z,[Nf(2:end),1]);        % joint posterior
        Q     = spm_margin(U);                     % marginal  posteriors

    otherwise
        disp('unknown method')

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


