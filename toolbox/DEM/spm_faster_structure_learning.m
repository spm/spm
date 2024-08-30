function [MDP] = spm_faster_structure_learning(O,S,dx,dt)
% RG structure learning of a hierarchical POMDP
% FORMAT [MDP] = spm_faster_structure_learning(O,S)
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
    m        = size(O{n},1)/prod(S{n});

    % Grouping into a partition of outcomes
    %----------------------------------------------------------------------
    G        = spm_group(S{n},m,dx(n));
    S{n + 1} = size(G);
    MDP{n}.T = dt(n);
    MDP{n}.G = G;


    % likelihood and transitions for each group
    %======================================================================
    Ng    = numel(G);
    Nt    = size(O{n},2);
    X     = cell(Ng,Nt - 1);
    P     = cell(Ng,Nt - 1);  
    for g = 1:Ng

        % structure_learning from unique exemplars
        %------------------------------------------------------------------
        mdp  = spm_structure_fast(O{n}(G{g},:));

        % likelihoods and priors for this group
        %------------------------------------------------------------------
        MDP{n}.A(G{g},1) = mdp.a;
        MDP{n}.B(g,1)    = mdp.b;

        % intial states and paths
        %------------------------------------------------------------------
        X(g,:)  = mdp.X;
        P(g,:)  = mdp.P;

        % parents of outcomes
        %------------------------------------------------------------------
        for j = 1:numel(G{g})
            MDP{n}.id.A{G{g}(j)} = uint16(g);
        end
        
    end

    % outcomes at next time scale
    %----------------------------------------------------------------------
    t        = 1:dt(n):(Nt - 1);
    O{n + 1} = [X(:,t); P(:,t)];

    % parents of states and paths
    %----------------------------------------------------------------------
    MDP{n}.id.D = num2cell(     (1:Ng));       % parents of state
    MDP{n}.id.E = num2cell(Ng + (1:Ng));       % parents of paths

    % Check whether there is only one group or (generalised) object
    %----------------------------------------------------------------------
    if (Ng == 1) || (Nt == 1), break, end

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
o       = fix(o*2);
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
Nt    = numel(j) - 1;
b     = false(Ns,Ns);
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

% Vectorise cell array of likelihood tensors and place in structure
%--------------------------------------------------------------------------
mdp.a    = a;
mdp.b{1} = b;

% add probabilities over inital states and paths
%==========================================================================
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
for t = 1:Nt
    if Nu > 1
        mdp.P{t} = squeeze(b(j(t + 1),j(t),:));
    else
        mdp.P{t} = true;
    end
end
mdp.X = mdp.X(1:Nt);

return

function g = spm_group(N,m, d)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT g = spm_group(N,m,[d])
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


