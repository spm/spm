function [MDP] = spm_faster_structure_learning(O,S,dx,dt)
% RG structure learning of a hierarchical POMDP
% FORMAT [MDP] = spm_faster_structure_learning(O,S,dx,dt)
% O{N,T} - probabilitic exemplars of paths (cell aray)
% S      - array size of streams (number of streams x 4)
%    S(:,1) - size of group (x)
%    S(:,2) - size of group (y)
%    S(:,3) - size of group (z)
%    S(:,4) - number of outcomes per group
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


% check for stream specification
%--------------------------------------------------------------------------
if size(S,2) < 4
    S(:,end + 1:4) = 1;
end

% Outcomes per group (e.g., TrueColor, eigenmodes or generalised states)
%--------------------------------------------------------------------------
if size(S,1) == 1
    S(4)  = size(O,1)/prod(S(1:3));
end

% options for model inversion (and evaluation)
%==========================================================================

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

    % for each sector or segregated stream
    %----------------------------------------------------------------------
    sg    = {};                          % squences
    N     = {};                          % outcomes for next level
    D     = [];                          % parents of intial states
    E     = [];                          % parents of intial paths
    Ng    = [];                          % number of groups
    Ns    = 0;                           % number of groups (cummulative)
    No    = 0;                           % number of outcomes (cummulative)
    for s = 1:size(S,1)

        % Grouping into a partition of outcomes
        %------------------------------------------------------------------
        G           = spm_group(S(s,:),dx(n));
        MDP{n}.T    = dt(n);
        MDP{n}.G{s} = G;

        % likelihood and transitions for each group
        %==================================================================
        Ng(s) = numel(G);
        Nt    = size(O,2);
        X     = cell(Ng(s),Nt - 1);
        P     = cell(Ng(s),Nt - 1);
        for g = 1:Ng(s)

            % structure_learning from unique exemplars
            %--------------------------------------------------------------
            gg         = No + G{g};
            fg         = Ns + g;
            [mdp,j]    = spm_structure_fast(O(gg,:));
            sg{s}(g,:) = j;

            % likelihoods and priors for this group
            %--------------------------------------------------------------
            MDP{n}.A(gg,1) = mdp.a;
            MDP{n}.B(fg,1) = mdp.b;

            MDP{n}.sA(gg)  = s;
            MDP{n}.sB(fg)  = s;

            % intial states and paths
            %--------------------------------------------------------------
            X(g,:) = mdp.X;
            P(g,:) = mdp.P;

            % parents of outcomes
            %--------------------------------------------------------------
            for j = 1:numel(gg)
                MDP{n}.id.A{gg(j)} = uint16(fg);
            end

        end

        % outcomes at next time scale
        %------------------------------------------------------------------
        t  = 1:dt(n):(Nt - 1);
        N  = [N; [X(:,t); P(:,t)]];

        % parents of states and paths
        %------------------------------------------------------------------
        iD = 2*Ns + (1:Ng(s));
        iE = 2*Ns + (1:Ng(s)) + Ng(s);
        D  = [D, iD];
        E  = [E, iE];

        % streams at next level
        %------------------------------------------------------------------
        S(s,:) = [size(G,[1,2,3]),2];
        Ns     = Ns + Ng(s);
        No     = No + G{end}(end);

    end

    % parents
    %----------------------------------------------------------------------
    MDP{n}.id.D = num2cell(D);                      % parents of state
    MDP{n}.id.E = num2cell(E);                      % parents of paths

    % Combine streams if all comprise one group
    %----------------------------------------------------------------------
    if n > 1
        si = 1;                                     % source stream
        for sj = 2:numel(Ng)                        % target stream

            % size of likelihood mappings
            %--------------------------------------------------------------
            fi = find(ismember(MDP{n}.sB,si));      % source factors
            fj = find(ismember(MDP{n - 1}.sB,sj));  % target factors
            t  = numel(sg{1}(1,:));                 % number of outcomes

            for i = 1:numel(fi)
                for j = 1:numel(fj)

                    % sizes
                    %------------------------------------------------------
                    Ni = size(MDP{n}.B{fi(i)},    1);
                    Nj = size(MDP{n - 1}.B{fj(j)},2);
                    Nu = size(MDP{n - 1}.B{fj(j)},3);

                    % prediction of inital states (D)
                    %------------------------------------------------------
                    A  = zeros(Nj,Ni);
                    for f = 1:t
                        g = MDP{n - 1}.id.D{fj(j)}(1);
                        a = MDP{n}.A{g}(:,sg{sj}(j,f));
                        A(:,sg{si}(i,f)) = A(:,sg{si}(i,f)) + a;
                    end

                    % append mappings and iD
                    %------------------------------------------------------
                    MDP{n}.A{end + 1} = spm_dir_norm(A);
                    MDP{n}.id.A{end + 1} = fi(i);
                    MDP{n - 1}.id.D{fj(j)}(end + 1) = numel(MDP{n}.A);

                    % prediction of initial paths (E)
                    %------------------------------------------------------
                    A  = zeros(Nu,Ni);
                    for f = 1:t
                        g = MDP{n - 1}.id.E{fj(j)}(1);
                        a = MDP{n}.A{g}(:,sg{sj}(j,f));
                        A(:,sg{si}(i,f)) = A(:,sg{si}(i,f)) + a;
                    end

                    % append mappings and iD
                    %------------------------------------------------------
                    MDP{n}.A{end + 1} = spm_dir_norm(A);
                    MDP{n}.id.A{end + 1} = fi(i);
                    MDP{n - 1}.id.E{fj(j)}(end + 1) = numel(MDP{n}.A);

                end
            end
        end
    end

    % if all streams comprise one group
    %======================================================================
    if true

        % remove trailing streams
        %------------------------------------------------------------------
        if max(Ng) == 1

            % remove trailing factors
            %--------------------------------------------------------------
            MDP{n}.B = MDP{n}.B(1);
            Na    = numel(MDP{n}.A);
            d     = false(1,Na);
            for g = 1:Na
                d(g) = ~any(MDP{n}.id.A{g} > 1);
            end
            i = find(d);

            % remove their children
            %--------------------------------------------------------------
            MDP{n}.A    = MDP{n}.A(d);
            MDP{n}.id.A = MDP{n}.id.A(d);

            % and update parents of subordinate factors
            %--------------------------------------------------------------
            for j = 1:numel(MDP{n - 1}.id.D)
                MDP{n - 1}.id.D{j} = find(ismember(i,MDP{n - 1}.id.D{j}));
            end
            for j = 1:numel(MDP{n - 1}.id.E)
                MDP{n - 1}.id.E{j} = find(ismember(i,MDP{n - 1}.id.E{j}));
            end

            break
        end

    else

        % Combine streams in a final level
        %------------------------------------------------------------------
        if max(Ng) == 1
            S = [1 1 1 2*sum(Ng)];

            % Check whether there is only one group
            %--------------------------------------------------------------
            if (numel(Ng) == 1) || (Nt == 1)
                break
            end
        end
    end

    % outcomes at next level
    %----------------------------------------------------------------------
    O  = N;

end

return


function [mdp,j] = spm_structure_fast(O)
% A fast form of structure learning
% FORMAT [mdp,j] = spm_structure_fast(O)
% O   - Cell array of (cells of) a sequence of probabilistic outcomes
%
% mdp - likelihood (a) and transition (b) tensors for this sequence
% j   - sequence of unique identifiers
%
% This routine emulates structure learning from deterministic sequences of
% observations or outcomes. Effectively, it identifies all unique
% combinations of outcomes and transitions among those combinations. The
% likelihood matrix then maps from all unique combinations (i.e., latent
% states) to outcomes, and the transition tensor encodes observed
% transitions among latent states.
%__________________________________________________________________________


% distance matrix (i.e., normalised vectors on a hypersphere)
%--------------------------------------------------------------------------
D       = spm_information_distance(O);

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[~,i,j] = unique(D < 1,'rows','stable');


% Likelihood tensors
%--------------------------------------------------------------------------
Ns    = numel(i);                           % number of latent causes
Ng    = size(O,1);                          % number in group
a     = cell(Ng,1);
for g = 1:Ng

    % use first mapping
    %--------------------------------------------------------------------------
    a{g} = full(spm_cat(O(g,i)));

%     % or average
%     %--------------------------------------------------------------------------
%     for s = 1:Ns
%         a{g}(:,s) = full(mean(spm_cat(O(g,ismember(j,s))),2));
%     end

end

% Transition tensors
%--------------------------------------------------------------------------
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

function g = spm_group(N,d)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT g = spm_group(N,d)
%--------------------------------------------------------------------------
% N  - size of space (Nx,Ny,Nz,m)
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
N((end + 1):4) = 1;
if nargin < 2

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
            n{1}     = sparse((1:d(1)*N(4)) + r{1}(i)*N(4),1,1,L(1)*N(4),1);
            n{2}     = sparse((1:d(2)  )    + r{2}(j)     ,1,1,L(2)     ,1);
            n{3}     = sparse((1:d(3)  )    + r{3}(k)     ,1,1,L(3)     ,1);

            v        = spm_cross(n{1},n{2},n{3});
            v        = v(1:N(1)*N(4),1:N(2),1:N(3));
            g{i,j,k} = find(v(:));
        end
    end
end

return


