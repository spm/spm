function [MDP,RG,S] = spm_fast_structure_learning(O,S,dx,dt)
% RG structure learning of a hierarchical POMDP
% FORMAT [MDP,RG,S] = spm_fast_structure_learning(O,S)
% O{N,T} - probabilitic exemplars of paths (cell aray)
% S      - size of pixel array (2 or 3-vector)
%
% dx     - space scaling [default: 2 or 3]
% dt     - time  scaling [default: 2]
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
OPTIONS.B = 0;
OPTIONS.N = 0;
OPTIONS.O = 0;
OPTIONS.Y = 0;

% options for model inversion (and evaluation)
%==========================================================================
S         = {S};                               % size of image
O         = {O};                               % outcomes

% scaling (RG flow)
%--------------------------------------------------------------------------
if nargin < 3,    dx = 2; end
if nargin < 4,    dt = 2; end
if numel(dx) < 2, dx = repmat(dx,1,16); end
if numel(dt) < 2, dt = repmat(dt,1,16); end

% learn the dynamics in the form of a hierarchical MDP
%==========================================================================
MDP   = {};
for n = 1:8

    % Outcomes per tile or group (e.g., TrueColor or eigenmodes)
    %----------------------------------------------------------------------
    m     = numel(O{n}(:,1))/prod(S{n});

    % Grouping into a partition of outcomes
    %----------------------------------------------------------------------
    G     = spm_tile(S{n},m,dx(n));
    T     = spm_time(size(O{n},2),dt(n));
    for g = 1:numel(G)

        % structure_learning from unique exemplars
        %------------------------------------------------------------------
        mdp  = spm_structure_fast(O{n}(G{g},:));

        % place in hierarchical structure
        %------------------------------------------------------------------
        MDP{n}.A(G{g},1) = mdp.a;
        MDP{n}.B(g,1)    = mdp.b;
        for j = 1:numel(G{g})
            MDP{n}.id.A{G{g}(j)} = uint16(g);   % state of children
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
        pdp   = spm_VB_XXX(pdp);

        % initial states and paths
        %------------------------------------------------------------------
        for g = 1:numel(G)
            MDP{n}.id.D{g} = uint16(2*g - 1);    % parents of state
            MDP{n}.id.E{g} = uint16(2*g - 0);    % parents of paths

            O{n + 1}{MDP{n}.id.D{g},t} = pdp.X{g}(:,1);
            O{n + 1}{MDP{n}.id.E{g},t} = pdp.P{g}(:,1);
        end
    end
    S{n + 1}    = size(G);
    RG{n}       = G;

    % Check whether there is only one group or (generalised) object
    %----------------------------------------------------------------------
    if (numel(G) == 1) || (numel(T) == 1)
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
o       = fix(o*4);
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
mdp.a    = spm_dir_norm(a);
mdp.b{1} = spm_dir_norm(b);

return

function g = spm_tile(N,m, d)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT g = spm_tile(N,m,[d])
%--------------------------------------------------------------------------
% N  - size of space
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
N(end + 1)       = 1;
if nargin < 2, m = 1; end
if nargin < 3

    % use 3 x 3 tiles (or smaller)
    %----------------------------------------------------------------------
    if ~rem(N(1),3), d(1) = 3; else, d(1) = 2; end
    if ~rem(N(2),3), d(2) = 3; else, d(2) = 2; end
    if ~rem(N(3),3), d(3) = 3; else, d(3) = 2; end

elseif numel(d) == 1
    d  = [d,d,d];
end

% deal with single row (or column) cases
%--------------------------------------------------------------------------
r     = cell(1,3);
s     = ones(1,3);
L     = ones(1,3);
for i = 1:3
    d(i)    = min(d(i),N(i));                 % block size
    r{i}    = 0:d(i):(N(i) - 1);              % block start
    s(i)    = numel(r{i});                    % length

end
for i = 1:3
    L(i)    = r{i}(end) + d(i);               % dim
end

% Decimate rows and columns
%--------------------------------------------------------------------------
g     = cell(s);
for i = 1:s(1)
    for j = 1:s(2)
        for k = 1:s(3)
            n{1}     = sparse((1:d(1)*m) + r{1}(i)*m,1,1,L(1)*m,1);
            n{2}     = sparse((1:d(2)  ) + r{2}(j)  ,1,1,L(2)  ,1);
            n{3}     = sparse((1:d(3)  ) + r{3}(k)  ,1,1,L(3)  ,1);

            v        = spm_cross(n{1},n{2},n{3});
            v        = v(1:N(1)*m,1:N(2),1:N(3));
            g{i,j,k} = find(v(:));
        end
    end
end

return


function t = spm_time(T,d)
% Grouping into a partition of non-overlapping sequences
% FORMAT t = spm_time(T,d)
% T  - total number of the timesteps
% d  - number timesteps per partition [default: 2]
%--------------------------------------------------------------------------
% Effectively, this enables a representation of generalised motion; that
% can also be read in terms of paths or trajectories
%--------------------------------------------------------------------------
if nargin < 2, d = 2; end
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

