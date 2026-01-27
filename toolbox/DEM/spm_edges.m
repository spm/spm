function [j,i,q] = spm_edges(id,g,Q)
% Returns parents and children of MDP likelihood mapping
% FORMAT [j,i,q] = spm_edges(id,g,Q)
% 
% id - identifier or index structure
%  id.A{g} State-independent domain
%
%  id.ff - List of domain factors
%  id.fg - List of parents  for A{g} under each combination of domains 
%  id.gg - List of children for A{g} under each combination of domains
% 
% g  - index of likelihood mapping A{g}
% Q  - posterior over domain factors (id.ff)
%
% j{}  -  parents of A{g}: hidden factors
% i{}  - children of A{g}: outcome modalities
% q    - posterior over parents and children
%
% Returns the domain [codomains] of factors [modalities] for this
% likelihood mapping. These [co] domains may or may not be state-dependent.
% If they are state-dependent, the posteriors are returned.
%__________________________________________________________________________
%           
% id.A{g} - Indices of parents of A{g}
% id.D{f} - Indices of parents of initial states of factor B{f}
% id.E{f} - Indices of parents of paths of factor B{f} 
%
% id.ff   - List of domain factors
% id.fg   - List of parents  of A{g} under each combination of domains 
% id.gg   - List of children of A{g} under each combination of domains
%
% id.g{k} - Partition of selected (i.e., attended) children of A{g}
% id.ge   - Indices of outcomes subtending expected free energy
%
% id.fu   - List of factors with contolled paths
% id.fp   - List of factors with inferred  paths 
%
% see: spm_parents.m
%__________________________________________________________________________


% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging

% If there are state-dependent domains
%--------------------------------------------------------------------------
if isfield(id,'ff')

    % Most likely states
    %----------------------------------------------------------------------
    Nff   = numel(id.ff);
    for f = 1:Nff
        ff    = id.ff(f);                      % domain factor
        r{f}  = find(Q{ff} > max(Q{ff})/16);   % likely states
        R{f}  = Q{ff}(r{f});                   % reduced states
        Nr(f) = numel(r{f});                   % number of reduced states
    end
    q     = spm_cross(R);
    q     = q(:);
    q     = q/sum(q);
    iq    = find(q > max(q)/16);
    q     = q(iq);

    % posterior over domain states
    %----------------------------------------------------------------------
    s     = cell(1,Nff);
    for k = 1:numel(q)

        % combination of domain states
        %------------------------------------------------------------------
        ind    = spm_index(Nr,iq(k));
        for ff = 1:Nff
            s{ff} = r{ff}(ind(ff));
        end

        % MAP domain (parents)
        %------------------------------------------------------------------
        if isfield(id,'fg')
            if iscell(id.fg)
                j{k} = id.fg{g}{s{:}};
            else
                j{k} = id.fg(g,[s{:}]);
            end
        else
            j{k} = id.A{g};
        end

        % MAP codomain (children)
        %------------------------------------------------------------------
        if isfield(id,'gg')
            if iscell(id.gg)
                i{k} = id.gg{g}{s{:}};
            else
                i{k} = id.gg(g,[s{:}]);
            end
        else
            i{k} = g;
        end
    end

else

    % State-independent domain
    %----------------------------------------------------------------------
    j = id.A(g);
    i = {g};
    q = 1;

end

