function MDP = spm_merge_structure_learning(O,MDP)
% RG structure learning of a hierarchical POMDP
% FORMAT MDP = spm_merge_structure_learning(O,MDP)
% O{N,T} - cell array of (probabilisitc) outcomes
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


% learn the dynamics in the form of a hierarchical MDP
%==========================================================================
for n = 1:numel(MDP)

    % for each sector or segregated stream
    %----------------------------------------------------------------------
    N     = {};                             % outcomes for next level
    Ns    = 0;                              % number of groups
    No    = 0;                              % number of outcomes
    for s = 1:numel(MDP{n}.G)

        % Grouping into a partition of outcomes
        %----------------------------------------------------------------------
        G     = MDP{n}.G{s};

        % likelihood and transitions for each group
        %======================================================================
        Ng    = numel(G);
        Nt    = size(O,2);
        X     = cell(Ng,Nt - 1);
        P     = cell(Ng,Nt - 1);
        for g = 1:Ng

            % structure_learning from unique exemplars
            %------------------------------------------------------------------
            gg   = No + G{g};
            fg   = Ns + g;
            mdp  = spm_merge_fast(O(gg,:),MDP{n}.A(gg),MDP{n}.B(fg));

            % likelihoods and priors for this group
            %------------------------------------------------------------------
            MDP{n}.A(gg) = mdp.a;
            MDP{n}.B(fg) = mdp.b;

            % intial states and paths
            %------------------------------------------------------------------
            X(g,:)  = mdp.X;
            P(g,:)  = mdp.P;

        end

        % outcomes at next time scale
        %------------------------------------------------------------------
        t  = 1:MDP{n}.T:(Nt - 1);
        N  = [N; [X(:,t); P(:,t)]];
        Ns     = Ns + Ng;
        No     = No + G{end}(end);

    end

    % break if no new (generalised) states
    %----------------------------------------------------------------------
    if size(N,2) < 1, break, end

    % outcomes at next level
    %----------------------------------------------------------------------
    O  = N;

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


% sizes
%--------------------------------------------------------------------------
B     = logical(B{1});             % old transitions
Ng    = size(A,1);                 % number of modalities
Ns    = size(B,1);                 % number of old states
Nu    = size(B,3);                 % number of old paths

% append new states to A
%--------------------------------------------------------------------------
for g = 1:Ng
    Na     = size(A{g,1},1);
    No     = size(O{g,1},1);
    if isnumeric(O{g,1})
        a  = zeros(No,Ns);
    else
        a  = false(No,Ns);
    end
    i      = 1:Na;
    a(i,:) = A{g,1};
    A{g,1} = a;
end

% Unique outputs
%==========================================================================

% distance matrix (i.e., normalised vectors on a hypersphere)
%--------------------------------------------------------------------------
[D,C]   = spm_information_distance([A O]);

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[~,i,j] = unique(D < 1,'rows','stable');

k     = 1:Ns;                      % old states
n     = i(~ismember(i,k));         % new states
k     = j;                         % old sequence
j     = j((Ns + 1):end);           % new sequence


% accumulate likelihood tensors
%--------------------------------------------------------------------------
a     = cell(Ng,1);
for g = 1:Ng

    % use first mapping
    %----------------------------------------------------------------------
    a{g} = full(spm_cat({A{g} O(g,n - Ns)}));

%     % or average
%     %----------------------------------------------------------------------
%     a{g}  = full(A{g});
%     for s = 1:numel(n)
%         a{g}(:,Ns + s) = mean(spm_cat(O(g,ismember(j,k(n(s))))),2);
%     end

end

% Transition tensors
%--------------------------------------------------------------------------
k     = 1:Ns;                      % old states
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



