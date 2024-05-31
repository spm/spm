function [MDP,RG,S] = spm_fast_structure_learning(O,S)
% RG structure learning of a hierarchical POMDP
% FORMAT [MDP,RG,S] = spm_fast_structure_learning(O,S)
% O{N,T} - probabilitic exemplars of paths (cell aray)
% S      - size of pixel array (2-vector)
%
% This routine returns a hierarchy of generative models (MDP{n}), given a
% sequence of outcomes using RG (renormalization group) operators.
%
% See: spm_MDP_structure_learning.m
%      spm_MDP_structure_teaching.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________

% options for model inversion (and evaluation)
%==========================================================================
dt    = 2;                              % time scaling
S     = {S};                            % size of image
O     = {O};                            % outcomes

% learn the dynamics in the form of a hierarchical MDP
%==========================================================================
MDP   = {};
for n = 1:8

    % Outcomes per tile or group
    %----------------------------------------------------------------------
    m     = numel(O{n}(:,1))/prod(S{n});

    % Grouping into a partition of outcomes
    %----------------------------------------------------------------------
    G     = spm_tile(S{n}(1),S{n}(2),m);
    T     = spm_time(size(O{n},2),dt);
    for g = 1:numel(G)

        % structure_learning from unique exemplars
        %------------------------------------------------------------------
        mdp  = spm_structure_fast(O{n}(G{g},:));

        % place in hierarchical structure
        %------------------------------------------------------------------
        MDP{n}.A(G{g},1) = mdp.a;
        MDP{n}.B(g,1)    = mdp.b;
        for j = 1:numel(G{g})
            MDP{n}.id.A{G{g}(j)} = uint16(g);
        end
    end

    % solve at next time scale
    %----------------------------------------------------------------------
    for t = 1:numel(T)
        pdp   = MDP{n};
        pdp.T = numel(T{t});
        for j = 1:pdp.T
            pdp.O(:,j) = O{n}(:,T{t}(j));
        end
        pdp   = spm_MDP_VB_XXX(pdp);

        % initial states and paths
        %------------------------------------------------------------------
        for g = 1:numel(G)
            MDP{n}.id.D{g} = 2*g - 1;
            MDP{n}.id.E{g} = 2*g - 0;

            O{n + 1}{MDP{n}.id.D{g},t} = pdp.X{g}(:,1);
            O{n + 1}{MDP{n}.id.E{g},t} = pdp.P{g}(:,end);
        end
    end
    S{n + 1}    = size(G);
    RG{n}       = G;

    % Check whether there is only one group or (generalised) object
    %----------------------------------------------------------------------
    if numel(G) == 1
        break
    end

end

return

function mdp = spm_structure_fast(O)
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

% Unique outputs
%--------------------------------------------------------------------------
o       = spm_cat(O)';
o       = fix(o*128);
[~,i,j] = unique(o,'rows','stable');

% Likelihood tensors
%--------------------------------------------------------------------------
Ng    = size(O,1);
a     = cell(Ng,1);
for g = 1:Ng
    a{g} = full(spm_cat(O(g,i)));
end

% Transition tensors
%--------------------------------------------------------------------------
Ns    = numel(i);
b     = zeros(Ns,Ns);
for t = 1:(numel(j) - 1)

    % find empty paths
    %----------------------------------------------------------------------
    if ~any(b(j(t + 1),j(t),:),'all')
        u  = find(~any(b(:,j(t),:),1),1,'first');
        if isempty(u)
            b(j(t + 1),j(t),end + 1) = 1;
        else
            b(j(t + 1),j(t),u) = 1;
        end
    end
end

% fill in empty columns of transition tensor
%--------------------------------------------------------------------------
% for u = 1:size(b,3)
%     for s = 1:Ns
%         if ~any(b(:,s,u))
%             i = find(any(squeeze(b(:,s,:)),2),1);
%             b(i,s,u) = 1;
%         end
%     end
% end

% Vectorise cell array of likelihood tensors and place in structure
%--------------------------------------------------------------------------
mdp.a    = a;
mdp.b{1} = b;

return

function g = spm_tile(Nr,Nc,m, d)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT g = spm_tile(Nr,Nc,m,[d])
%--------------------------------------------------------------------------
% Nr - Number of outcome rows (e.g., pixels)
% Nc - Number of outcome columns
% m  - Number of modalities per outcome (e.g., pixel) [default: 1]
% d  - Number of rows per tile [default: 2 or 3] 
%
% g  - Cell array of outcome indices
%
% Effectively, this leverages the conditional independencies that inherit
% from local interactions; of the kind found in metric spaces that preclude
% action at a distance.
%--------------------------------------------------------------------------

% defaults
%--------------------------------------------------------------------------
if nargin < 3, m = 1; end
if nargin < 4

    % use 3 x 3 tiles (or smaller)
    %----------------------------------------------------------------------
    if ~rem(Nr,3), dr = 3; else, dr = 2; end
    if ~rem(Nc,3), dc = 3; else, dc = 2; end

else
    dr = d;
    dc = d;
end

% deal with single row (or column) cases
%--------------------------------------------------------------------------
dr    = min(dr,Nr);
dc    = min(dc,Nc);

% Decimate rows and columns
%--------------------------------------------------------------------------
r     = 0:dr:(Nr - dr);
c     = 0:dc:(Nc - dc);
for i = 1:numel(r)
    for j = 1:numel(c)
        n = sparse((1:dr*m) + r(i)*m,1,1,Nr*m,1)*sparse((1:dc) + c(j),1,1,Nc,1)';
        g{i,j} = find(n(:));
    end
end

return


function t = spm_time(T,d)
% Grouping into a partition of non-overlapping sequences
% FORMAT t = spm_time(T,d)
% T  - total number of the timesteps
% d  - number timesteps per partition
%--------------------------------------------------------------------------
% Effectively, this enables a representation of generalised motion; that
% can also be read in terms of paths or trajectories
%--------------------------------------------------------------------------
for i = 1:floor(T/d)
    t{i} = (1:d) + (i - 1)*d;
end

return

function mdp = spm_structure(O)
% structure learning (UNUSED)
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
%
% This version tests for changes in mutual information or expected free
% energy by assigning new outcomes to the most likely previously encountered
% state and a new state and testing for increases or decreases in mutual
% information (with prior Dirichlet concentration parameters of a suitably
% small value: p)
%__________________________________________________________________________

% Unique outputs
%--------------------------------------------------------------------------
o       = spm_cat(O)';
[~,i,j] = unique(o,'rows','stable');

% Likelihood tensors
%--------------------------------------------------------------------------
Ng    = size(O,1);
a     = cell(Ng,1);
n     = sum(j == 1);
for g = 1:Ng
    a{g} = full(n*O{g,i(1)});
end

% Dirichlet concentration parameter
%--------------------------------------------------------------------------
p     = 1/32;
I     = 1;
for k = 2:numel(i)

    % find closest predecessor
    %----------------------------------------------------------------------
    No    = sum(j == k);
    Ns    = size(a{1},2);
    KL    = zeros(Ng,Ns);
    for g = 1:Ng
        KL(g,:) = O{g,i(k)}'*spm_psi(a{g});
    end
    [~,m] = max(sum(KL));

    % compressed and uncompressed likelihood mappings
    %----------------------------------------------------------------------
    EFE   = zeros(Ng,1);
    for g = 1:Ng
        a0{g} = a{g};
        a0{g}(:,end + 1) = 0;
        a0{g} = a0{g} + p;
        a1{g} = a0{g};

        % information loss
        %------------------------------------------------------------------
        da           = No*O{g,i(k)};
        a0{g}(:,m)   = a0{g}(:,m)   + da;
        a1{g}(:,end) = a1{g}(:,end) + da;

        EFE(g)       = spm_MDP_MI(a0{g}) - spm_MDP_MI(a1{g});
    end

    % test for mutual information (expected free energy)
    %----------------------------------------------------------------------
    if all(EFE > 0)
        for g = 1:Ng
            a{g} = a0{g} - p;
            I(k) = m;
        end
    else
        for g = 1:Ng
            a{g} = a1{g} - p;
            I(k) = Ns + 1;
        end
    end

end

% Transition tensors
%--------------------------------------------------------------------------
Ns    = size(a{1},2);
b     = zeros(Ns,Ns);
J     = I(j);
for t = 1:(numel(J) - 1)

    % find previous transitions
    %----------------------------------------------------------------------
    v  = find(b(J(t + 1),J(t),:),1,'first');

    % new transition
    %----------------------------------------------------------------------
    if isempty(v)

        % new path
        %------------------------------------------------------------------
        u  = find(~any(b(:,J(t),:),1),1,'first');
        if isempty(u)
            b(J(t + 1),J(t),end + 1) = 1;
        else

            % old path
            %--------------------------------------------------------------
            b(J(t + 1),J(t),u) = b(J(t + 1),J(t),u) + 1;
        end

    else

        % old transition
        %--------------------------------------------------------------
        b(J(t + 1),J(t),v) = b(J(t + 1),J(t),v) + 1;
    end
end


% Vectorise cell array of likelihood tensors and place in structure
%--------------------------------------------------------------------------
mdp.a    = a;
mdp.b{1} = b;

return

