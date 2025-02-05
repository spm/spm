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
    sg    = {};                                     % squences
    N     = {};                                     % next level outcomes
    for s = 1:numel(MDP{n}.G)

        % temporal down-sampling
        %------------------------------------------------------------------
        t     = 1:MDP{n}.T:(size(O,2) - 1);

        % likelihood and transitions for each group
        %==================================================================
        G     = MDP{n}.G{s};                        % groups in stream
        Ng(s) = numel(G); 
        for g = 1:Ng(s)

            % structure_learning from unique exemplars
            %--------------------------------------------------------------
            gg         = G{g};
            fg         = MDP{n}.id.A{min(gg)};
            [mdp,j]    = spm_merge_fast(O(gg,:),MDP{n}.a(gg),MDP{n}.b(fg));
            sg{s}(g,:) = j;

            % likelihoods and priors for this group
            %--------------------------------------------------------------
            MDP{n}.a(gg,1) = mdp.a;
            MDP{n}.b(fg,1) = mdp.b;

            % initial states and paths : outcomes at next time scale
            %--------------------------------------------------------------
            iD = min(MDP{n}.id.D{fg});                % parents of state
            iE = min(MDP{n}.id.E{fg});                % parents of paths
            if numel(iD)
                N(iD,:) = mdp.X(:,t);
                N(iE,:) = mdp.P(:,t);
            end
            
        end
    end

    % Link principal (leading) stream to trading streamss
    %======================================================================
    if n > 1
        si = 1;                                      % source stream
        st = find(Ng);                               % target stream
        for sj = st(2:end)

            % size of likelihood mappings
            %--------------------------------------------------------------
            fsi = find(ismember(MDP{n}.sB,si));      % source factors
            fsj = find(ismember(MDP{n - 1}.sB,sj));  % target factors
            t   = numel(sg{1}(1,:));                 % number of outcomes

            for i = 1:numel(fsi)
                for j = 1:numel(fsj)

                    % sizes
                    %------------------------------------------------------
                    fi = fsi(i);
                    fj = fsj(j);
                    Ni = size(MDP{n}.b{fi},    1);
                    Nj = size(MDP{n - 1}.b{fj},2);
                    Nu = size(MDP{n - 1}.b{fj},3);

                    % prediction of initial states (D)
                    %------------------------------------------------------
                    gi = MDP{n - 1}.ss.D{si,sj}(fi,fj);
                    gj = MDP{n - 1}.ss.D{sj,sj}(fj,fj);

                    % augment likelihood mappings
                    %--------------------------------------------------
                    A  = spm_dir_norm(MDP{n}.a{gj});
                    a  = zeros(Nj,Ni);
                    [x,y] = size(MDP{n}.a{gi});
                    a(1:x,1:y) = MDP{n}.a{gi};
                    for f = 1:t
                        ii = sg{si}(i,f);
                        ij = sg{sj}(j,f);
                        a(:,ii) = a(:,ii) + A(:,ij);
                    end
                    MDP{n}.a{gi} = a;


                    % prediction of initial paths (E)
                    %------------------------------------------------------
                    gi = MDP{n - 1}.ss.E{si,sj}(fi,fj);
                    gj = MDP{n - 1}.ss.E{sj,sj}(fj,fj);

                    % augment likelihood mappings
                    %--------------------------------------------------
                    A  = spm_dir_norm(MDP{n}.a{gj});
                    a  = zeros(Nu,Ni);
                    [x,y] = size(MDP{n}.a{gi});
                    a(1:x,1:y) = MDP{n}.a{gi};
                    for f = 1:t
                        ii = sg{si}(i,f);
                        ij = sg{sj}(j,f);
                        a(:,ii) = a(:,ii) + A(:,ij);
                    end
                    MDP{n}.a{gi} = a;
                    
                end
            end
        end
    end

    % compress tensors with no parents
    %----------------------------------------------------------------------
    for f = 1:numel(MDP{n}.b)
        if isempty(MDP{n}.id.D{f})

            % compress prior transitions
            %--------------------------------------------------------------
            MDP{n}.b{f} = sum(MDP{n}.b{f},'all');

            % and children
            %--------------------------------------------------------------
            for g = find(ismember([MDP{n}.id.A{:}],f))
                MDP{n}.a{g} = sum(MDP{n}.a{g},2);
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

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[i,j] = spm_unique([A,O]);

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



