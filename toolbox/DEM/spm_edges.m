function [j,i,q] = spm_edges(id,g,Q)
% Returns parents and children of MDP likelihood mapping
% FORMAT [j,i,q] = spm_get_edges(id,g,Q)
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
%
% see: spm_MDP_VB_XXX.m; spm_get_edges
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging

% If there are state-dependent domains
%--------------------------------------------------------------------------
if isfield(id,'ff')

    % Most likely states
    %----------------------------------------------------------------------
    q   = spm_cross(Q{id.ff});
    Ns  = size(q);
    q   = q(:);
    iq  = find(q > max(q)/16);

    % posterior over domain states
    %----------------------------------------------------------------------
    q   = q(iq);
    q   = q/sum(q);

    for k = 1:numel(q)

        % combination of domain states
        %------------------------------------------------------------------
        s   = spm_index(Ns,iq(k));
        s   = num2cell(s);

        % MAP domain (parents)
        %------------------------------------------------------------------
        if isfield(id,'fg')
            j{k} = id.fg{g}{s{:}};
        else
            j{k} = id.A{g};
        end

        % MAP codomain (children)
        %------------------------------------------------------------------
        if isfield(id,'gg')
            i{k} = id.gg{g}{s{:}};
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

