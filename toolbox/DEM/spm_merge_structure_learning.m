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
    sg    = {};                             % squences
    N     = {};                             % outcomes for next level
    Ns    = 0;                              % number of groups
    No    = 0;                              % number of outcomes
    for s = 1:numel(MDP{n}.G)

        % Group into a partition of outcomes
        %------------------------------------------------------------------
        G     = MDP{n}.G{s};

        % likelihood and transitions for each group
        %==================================================================
        Ng    = numel(G);
        Nt    = size(O,2);
        X     = cell(Ng,Nt - 1);
        P     = cell(Ng,Nt - 1);
        for g = 1:Ng

            % structure_learning from unique exemplars
            %--------------------------------------------------------------
            gg         = No + G{g};
            fg         = Ns + g;
            [mdp,j]    = spm_merge_fast(O(gg,:),MDP{n}.a(gg),MDP{n}.b(fg));
            sg{s}(g,:) = j;

            % likelihoods and priors for this group
            %--------------------------------------------------------------
            MDP{n}.a(gg,1) = mdp.a;
            MDP{n}.b(fg,1) = mdp.b;

            % intial states and paths
            %--------------------------------------------------------------
            X(g,:)  = mdp.X;
            P(g,:)  = mdp.P;

        end

        % outcomes at next time scale
        %------------------------------------------------------------------
        t  = 1:MDP{n}.T:(Nt - 1);
        N  = [N; [X(:,t); P(:,t)]];
        Ns = Ns + Ng;
        No = No + G{end}(end);

    end

    % Link principal (leading) stream to trading streamss
    %======================================================================
    if n > 1
        gi = sum(MDP{n}.sA == MDP{n}.sC) + 1;
        si = 1;                                     % source stream
        for sj = 2:numel(MDP{n}.G)                  % target stream

            % size of likelihood mappings
            %--------------------------------------------------------------
            fi = find(ismember(MDP{n}.sB,si));      % source factors
            fj = find(ismember(MDP{n - 1}.sB,sj));  % target factors
            t  = numel(sg{1}(1,:));                 % number of outcomes

            for i = 1:numel(fi)
                for j = 1:numel(fj)

                    % sizes
                    %------------------------------------------------------
                    Ni = size(MDP{n}.b{fi(i)},    1);
                    Nj = size(MDP{n - 1}.b{fj(j)},2);
                    Nu = size(MDP{n - 1}.b{fj(j)},3);

                    % prediction of initial states (D)
                    %------------------------------------------------------
                    g  = MDP{n - 1}.id.D{fj(j)}(1);
                    A  = spm_dir_norm(MDP{n}.a{g});
                    a  = zeros(Nj,Ni);

                    % augment likelihood mappings
                    %------------------------------------------------------
                    [x,y] = size(MDP{n}.a{gi});
                    a(1:x,1:y) = MDP{n}.a{gi};
                    for f = 1:t
                        ii = sg{si}(i,f);
                        ij = sg{sj}(j,f);
                        a(:,ii) = a(:,ii) + A(:,ij);
                    end
                    MDP{n}.a{gi} = a;
                    gi = gi + 1;


                    % prediction of initial paths (E)
                    %------------------------------------------------------
                    g  = MDP{n - 1}.id.E{fj(j)}(1);
                    A  = spm_dir_norm(MDP{n}.a{g});
                    a  = zeros(Nu,Ni);

                    % augment likelihood mappings
                    %------------------------------------------------------
                    [x,y] = size(MDP{n}.a{gi});
                    a(1:x,1:y) = MDP{n}.a{gi};
                    for f = 1:t
                        ii = sg{si}(i,f);
                        ij = sg{sj}(j,f);
                        a(:,ii) = a(:,ii) + A(:,ij);
                    end
                    MDP{n}.a{gi} = a;
                    gi = gi + 1;

                end
            end
        end
    end

    % break if no new (generalised) states
    %----------------------------------------------------------------------
    if size(N,2) < 1, break, end

    % outcomes at next level
    %----------------------------------------------------------------------
    O  = N;

end

return


function [mdp,j] = spm_merge_fast(O,A,B)
% A fast form of structure learning
% FORMAT [mdp,j] = spm_merge_fast(O,A,B)
% O   - Cell array of (cells of) a sequence of probabilistic outcomes
% A   - Previouly encountered outcomes:    i.e., likelihood
% B   - Previouly encountered transitions: i.e., priors
%
% mdp - likelihood (a) and transition (b) tensors for this sequence
% j   - sequence of unique identifiers for new outcomes
%
% This routine emulates structure learning from deterministic sequences of
% observations or outcomes. Effectively, it identifies all unique
% combinations of outcomes and transitions among those combinations. The
% likelihood matrix then maps from all unique combinations (i.e., latent
% states) to outcomes, and the transition tensor encodes observed
% transitions among latent states. This version appends some now outcomes
% to previously encoutered outcomes encoded in A and B.
%__________________________________________________________________________


% sizes
%--------------------------------------------------------------------------
B  = B{1};                      % old transitions
Ng = size(A,1);                 % number of modalities
Ns = size(B,1);                 % number of old states
Nu = size(B,3);                 % number of old paths

% append new states to each a{g}
%--------------------------------------------------------------------------
for g = 1:Ng
    Na     = size(A{g,1},1);
    No     = size(O{g,1},1);
    a      = zeros(No,Ns);
    i      = 1:Na;
    a(i,:) = A{g,1};
    A{g,1} = a;
end

% Unique outputs
%==========================================================================

% distance matrix (i.e., normalised vectors on a hypersphere)
%--------------------------------------------------------------------------
D       = spm_information_distance([A O]);

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[~,i,j] = unique(D < sqrt(2),'rows','stable');


% accumulate likelihood tensors
%--------------------------------------------------------------------------
Nn    = numel(i);               % number of new states
j     = j((Ns + 1):end);
a     = cell(Ng,1);
for g = 1:Ng

    % previously encountered outcomes
    %----------------------------------------------------------------------
    No   = numel(O{g});
    a{g} = zeros(No,Nn); a{g}(1:No,1:Ns) = A{g};

    % add new outcomes
    %----------------------------------------------------------------------
    for t = 1:numel(j)
        a{g}(:,j(t)) = a{g}(:,j(t)) + O{g,t};
    end

end

% Transition tensors
%--------------------------------------------------------------------------
Nt    = numel(j) - 1;
b     = zeros(Nn,Nn,Nu); b(1:Ns,1:Ns,:) = B;
for t = 1:Nt

    % Is this an existing transition?
    %----------------------------------------------------------------------
    u  = find(b(j(t + 1),j(t),:),1,'first');
    if numel(u)

        % accumulate Drichlet counts for this transition
        %------------------------------------------------------------------
        b(j(t + 1),j(t),u) = b(j(t + 1),j(t),u) + 1;

    else

        % find the first path that does not have a successor of j(t)
        %------------------------------------------------------------------
        u  = find(~any(b(:,j(t),:),1),1,'first');
        if numel(u)

            % equip this path with a successor
            %--------------------------------------------------------------
            b(j(t + 1),j(t),u) = 1;

        else

            % otherwise create a new path
            %--------------------------------------------------------------
            b(j(t + 1),j(t),end + 1) = 1;
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
        mdp.P{t} = logical(squeeze(b(j(t + 1),j(t),:)));
    else
        mdp.P{t} = true;
    end
end
mdp.X = mdp.X(1:Nt);

return



