function [j,i] = spm_get_edges(id,g,Q)
% Returns parents and children of MDP likelihood mapping
% FORMAT [j,i] = spm_get_edges(id,g,Q)
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
% j  -   domain of A{g}: hidden factors     (i.e., parents)
% i  - codomain of A{g}: outcome modalities (i.e., children)
%
% Returns the domain [codomains] of factors [modalities] for this
% likelihood mapping. These [co] domains may or may not be state-dependent.
% If they are state-dependent, the maximum a posteriori [co] domain is
% returned.
%
% see: spm_MDP_VB_XXX.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging

% If there are state-dependent domains
%--------------------------------------------------------------------------
if isfield(id,'ff')

    % Most likely states
    %----------------------------------------------------------------------
    if iscell(Q)
        Ns    = numel(id.ff);
        s     = cell(Ns,1);
        for f = 1:Ns
            [~,m] = max(Q{id.ff(f)});
            s{f}  = m;
        end
    else
        s = num2cell(Q(id.ff));
    end

    % MAP domain (parents)
    %----------------------------------------------------------------------
    if isfield(id,'fg')
        j = id.fg{g}{s{:}};
    else
        j = id.A{g};
    end

    % MAP codomain (children)
    %----------------------------------------------------------------------
    if isfield(id,'gg')
        i = id.gg{g}{s{:}};
    else
        i = g;
    end

else

    % State-independent domain
    %----------------------------------------------------------------------
    j = id.A{g};
    i = g;

end

