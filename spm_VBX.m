function [P,F] = spm_VBX(O,P,A,id)
% vvariational Bayes estimate of categorical posterior over factors
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
% of latent states. The optional argument METHOD [default: full] switches
% among number of approximate schemes:
%
% 'full'    :  a vanilla variational scheme that uses a coordinate descent
% over a small number (8) iterations (fixed point iteration). This is
% computationally efficient, using only nontrivial posterior probabilities
% and the domain of the likelihood mapping in tensor operations.
% 
%
% 'exact'   :  a non-iterative heuristic but numerically accurate scheme
% that replaces the variational density over hidden factors with the
% marginal over the exact posterior
%
% 'sparse'  :  as for the exact scheme but suitable for sparse tensors
%
% 'marginal':  a heuristic scheme  that uses the log of the marginalised
% likelihood and log prior to estimate the log posterior
%
% see: spm_MDP_VB_XXX.m (NOTES)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
METHOD = 'full';                             % belief propagation scheme
ig     = id.g{id.ig(end)}(:)';               % modalities to sample

switch METHOD

    case 'full'

        %  (iterative) variational scheme
        %==================================================================

        % log prior : i{f} plausible states of factor f
        %------------------------------------------------------------------
        Nf    = numel(P);
        for f = 1:Nf
            s{f}  = find(P{f} > exp(-8));
            Q{f}  = P{f}(s{f});
            LP{f} = spm_vec(spm_log(Q{f}));
        end


        % Update domain factors Q{f}, if specified
        %==================================================================
        if isfield(id,'ff')

            % plausible domains id.fg{g}
            %--------------------------------------------------------------
            Nd    = numel(id.ff);
            for g = 1:numel(id.fg)
                fg{g} = id.fg{g}(s{id.ff});
            end
            Ns    = size(fg{1});

            % likelihood of domain factors
            %--------------------------------------------------------------
            L     = ones(Ns);
            for g = 1:numel(fg)
                for i = 1:numel(fg{g})
                    f    = fg{g}{i};
                    L(i) = L(i)*spm_dot(A{g}(:,s{f}),{O{g},Q{f}});
                end
            end

            % prior
            %--------------------------------------------------------------
            R     = cell(1,Nd);
            for i = 1:Nd
                f    = id.ff(i);
                R{i} = P{f}(s{f});
            end

            % posterior
            %--------------------------------------------------------------
            L     = L.*reshape(spm_cross(R{:}),Ns);
            Z     = sum(L,'all');
            L     = spm_marginal(L/Z);
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
            j     = spm_get_edges(id,g,P);
            j     = unique(j,'stable');
            LL    = spm_log(spm_dot(A{g}(:,s{j}),O{g}));
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


    case 'exact'

        % belief propagation with marginals of exact posterior
        %==================================================================

        % prior
        %------------------------------------------------------------------
        Nf    = numel(P);
        for f = 1:Nf
            s{f}  = find(P{f} > exp(-8));
            R{f}  = P{f}(s{f});
            Ns(f) = numel(P{f});
        end

        % accumulate likelihoods over modalities
        %------------------------------------------------------------------
        L     = 0;
        for g = ig
            j  = unique(id.A{g},'stable');
            LL = spm_dot(A{g}(:,s{j}),O{g});
            if numel(j) > 1
                [j,i] = sort(j);
                LL    = permute(LL,i);
            end
            k  = ones(1,Nf + 1); k(j) = size(LL,1:numel(j));
            L  = times(L,reshape(LL,k));
        end

        % marginal posteriors and free energy (partition function)
        %------------------------------------------------------------------
        U     = L.*spm_cross(R);                   % posterior unnormalised
        Z     = sum(U,'all');                      % partition coefficient
        if Z
            F = spm_log(Z);                        % negative free energy
            Q = spm_marginal(U/Z);                 % marginal  posteriors
        else
            F = -32;                               % negative free energy
            Q = spm_marginal(U + 1/numel(U));      % marginal  posteriors
        end

        % (-ve) free energy (partition coefficient)
        %------------------------------------------------------------------
        F     = 0;
        for f = 1:numel(P)
           LL = spm_vec(spm_dot(spm_log(L),Q,f));
           LP = spm_vec(spm_log(R{f}));
           F  = F + Q{f}'*(LL + LP - spm_log(Q{f}));
        end

        % Posterior
        %------------------------------------------------------------------
        for f = 1:numel(P)
            P{f}(:)    = 0;
            if numel(s{f}) > 1
                P{f}(s{f}) = Q{f};
            else
                P{f}(s{f}) = 1;
            end
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
        Q     = spm_marginal(U);                   % marginal  posteriors

    otherwise
        disp('unknown method')

end

return