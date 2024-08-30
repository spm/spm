function MDP = spm_merge_structure_learning(O,MDP)
% RG structure learning of a hierarchical POMDP
% FORMAT MDP = spm_merge_structure_learning(O,MDP)
% O{N,T} - probabilitic exemplars of paths (cell aray)
% MDP    - cell array of (renormalising) models
%
% This routine uses fast structure learning to append or merge new
% trajectories to a recursive or renormalising generative model specified
% as a cell array of Markov decision processes. It augments past
% (generalised) states and transitions with new states and transitions
% under the assumption that the first states have been seen before.
%
% See: spm_fast_structure_learning.m
%      spm_faster_structure_teaching.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% options for model inversion (and evaluation)
%==========================================================================
O         = {O};                               % outcomes

% learn the dynamics in the form of a hierarchical MDP
%==========================================================================
for n = 1:numel(MDP)


    % Grouping into a partition of outcomes
    %----------------------------------------------------------------------
    G     = MDP{n}.G;

    % likelihood and transitions for each group
    %======================================================================
    Ng    = numel(G);
    Nt    = size(O{n},2);
    X     = cell(Ng,Nt - 1);
    P     = cell(Ng,Nt - 1);  
    for g = 1:Ng

        % structure_learning from unique exemplars
        %------------------------------------------------------------------
        mdp  = spm_merge_fast(O{n}(G{g},:),MDP{n}.A(G{g}),MDP{n}.B(g));

        % likelihoods and priors for this group
        %------------------------------------------------------------------
        MDP{n}.A(G{g}) = mdp.a;
        MDP{n}.B(g)    = mdp.b;

        % intial states and paths
        %------------------------------------------------------------------
        X(g,:)  = mdp.X;
        P(g,:)  = mdp.P;
        
    end

    % outcomes at next time scale
    %----------------------------------------------------------------------
    t        = 1:MDP{n}.T:(Nt - 1);
    O{n + 1} = [X(:,t); P(:,t)];

    % break if no new (generalised) states
    %----------------------------------------------------------------------
    if size(O{n + 1},2) < 1, break, end

end

return


function mdp = spm_merge_fast(O,A,B)
% A fast form of structure learning
% FORMAT mdp = spm_structure(O)
% O   - Cell array of (cells of) a sequence of probabilistic outcomes
% mdp - likelihood (a) and transition (b) tensors for this sequence
%
% This routine emulates structure learning from deterministic sequences of
% observations or outcomes. Effectively, it identifies all unique
% combinations of outcomes and transitions among those combinations. The
% likelihood matrix then maps from all unique combinations (i.e., latent
% states) to outcomes, and the transition tensor encodes observed
% transitions among latent states.
%__________________________________________________________________________


% Sizes
%--------------------------------------------------------------------------
B     = logical(B{1});             % old transitions
Ng    = size(A,1);                 % number of modalities
Ns    = size(B,1);                 % number of old states
Nu    = size(B,3);                 % number of old paths
for g = 1:Ng
    Na     = size(A{g,1},1);
    No     = size(O{g,1},1);
    if isnumeric(O{g,1})
        a      = zeros(No,Ns);
    else
        a      = false(No,Ns);
    end
    i      = 1:Na;
    a(i,:) = A{g,1};
    A{g,1} = a;
end

% Unique outputs
%--------------------------------------------------------------------------
o       = spm_cat([A O])';
o       = fix(o*4);
[~,i,j] = unique(o,'rows','stable');

k     = 1:Ns;                      % old states
n     = i(~ismember(i,k));         % new states
j(k)  = [];                        % and sequence

% Likelihood tensors
%--------------------------------------------------------------------------
a     = cell(Ng,1);
for g = 1:Ng
    a{g} = full(spm_cat([A{g} O(g,n - Ns)]));
end


% Transition tensors
%--------------------------------------------------------------------------
Nn    = size(a{1},2);              % number of new states
b     = false(Nn,Nn,Nu); b(k,k,:) = B;
Nt    = numel(j) - 1;
for t = 1:Nt

    % find paths
    %----------------------------------------------------------------------
    if ~ any(b(j(t + 1),j(t),:),'all')
        u  = find(~any(b(:,j(t),:),1),1,'first');
        if isempty(u)
            b(j(t + 1),j(t),end + 1) = true;
        else
            b(j(t + 1),j(t),u) = true;
        end
    end
end

% place likelihood and prior tensors in structure
%--------------------------------------------------------------------------
mdp.a    = a;
mdp.b{1} = b;

% add probabilities over inital states and paths
%==========================================================================
Ns    = size(b,2);
Nu    = size(b,3);
X     = false(Ns,1);
mdp.X = cell(1,Nt);
mdp.P = cell(1,Nt);

% states
%--------------------------------------------------------------------------
for t = 1:(Nt + 1)
    s        = X;
    s(j(t))  = true;
    mdp.X{t} = s;
end

% paths
%--------------------------------------------------------------------------
mdp.P = cell(1,Nt);
for t = 1:Nt
    if Nu > 1
        mdp.P{t} = squeeze(b(j(t + 1),j(t),:));
    else
        mdp.P{t} = true;
    end
end
mdp.X = mdp.X(1:Nt);

return



