function [MDP] = spm_faster_structure_learning(O,S,dx,dt)
% RG structure learning of a hierarchical POMDP
% FORMAT [MDP] = spm_faster_structure_learning(O,S,dx,dt)
% O{N,T} - probabilitic exemplars of paths (cell aray)
% S      - array size of streams (number of streams x 4)
%    S(:,1) - size (x) of group in each stream  
%    S(:,2) - size (y) of group in each stream 
%    S(:,3) - size (z) of group in each stream 
%    S(:,4) - number of outcomes per group in each stream 
% dx     - space scaling [default: 2]
% dt     - time  scaling [default: 2]
%
% This routine returns a hierarchy of generative models (MDP{n}), given a
% sequence of outcomes using RG (renormalization group) operators. Output
% modalities are segregated into distinct streams (as specified by the size
% of groups:  S). The likelihood mappings within each stream maximise
% mutual information between successive levels by retaining unique
% combinations of outcomes over groups within each stream. These unique
% combinations are identified using the information distance between
% outcomes at each level. Between stream likelihoods map from the first or
% principal stream to subsequent or trailing streams, which will usually be
% the consequences of generalised states generated in the principal stream.
%
% See also: spm_MDP_structure_learning.m
%           spm_MDP_structure_teaching.m
%           spm_merge_structure_learning
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


SPINBLOCK = false;

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
if nargin < 3, dx = 2; end
if nargin < 4, dt = 2; end
dx = [dx repmat(dx(end),1,16)];
dt = [dt repmat(dt(end),1,16)];

% learn the dynamics in the form of a hierarchical MDP
%==========================================================================
MDP   = {};
for n = 1:8

    % for each sector or segregated stream
    %----------------------------------------------------------------------
    sg    = {};                           % squences
    N     = {};                           % outcomes for next level
    Ng    = [];                           % number of groups

    MDP{n}.a    = {};                     % likelihoods
    MDP{n}.b    = {};                     % transition priors
    MDP{n}.id.A = {};                     % parents of outcomes
    MDP{n}.id.D = {};                     % parents of intial states
    MDP{n}.id.E = {};                     % parents of intial paths

    % coupling among streams
    %----------------------------------------------------------------------
    MDP{n}.ss.D  = cell(size(S,1),size(S,1));
    MDP{n}.ss.E  = cell(size(S,1),size(S,1));
    MDP{n}.ss.ID = cell(size(S,1),size(S,1));
    MDP{n}.ss.IE = cell(size(S,1),size(S,1));

    
    for s = 1:size(S,1)

        % sizes
        %------------------------------------------------------------------
        No    = numel(MDP{n}.a);          % number of groups  (cummulative)
        Ns    = numel(MDP{n}.b);          % number of factors (cummulative)
        Nt    = size(O,2);                % number of time points
        t     = 1:dt(n):(Nt - 1);         % decimated time points

        % Grouping into a partition of outcomes
        %------------------------------------------------------------------
        if SPINBLOCK

            % blocking based on location
            %--------------------------------------------------------------
            G = spm_group(S(s,:),dx(n));
            
        else

            % partition based on mutual information
            %--------------------------------------------------------------
            o = [0; prod(S,2)];
            o = o(s) + (1:prod(S(s,:)));
            d = dx(n);
            m = S(s,4);
            G = spm_rgm_group(O(o,:),d,m);

        end

        MDP{n}.T    = dt(n);
        MDP{n}.G{s} = spm_unvec(spm_vec(G) + No,G);

        % likelihood and transitions for each group
        %==================================================================
        Ng(s) = numel(G);                             % groups in stream
        for g = 1:Ng(s)

            % structure_learning from unique exemplars
            %--------------------------------------------------------------
            fg         = Ns + g;
            gg         = MDP{n}.G{s}{g};
            [mdp,j]    = spm_structure_fast(O(gg,:));
            sg{s}(g,:) = j;

            % likelihoods and priors for this group
            %--------------------------------------------------------------
            MDP{n}.a(gg,1)  = mdp.a;
            MDP{n}.b(fg,1)  = mdp.b;

            % parents
            %--------------------------------------------------------------
            iD    = 2*Ns + 2*g - 1;
            iE    = 2*Ns + 2*g - 0;
            MDP{n}.id.A(gg) = {fg};                  % parents of outcomes                  
            MDP{n}.id.D(fg) = {iD};                  % parents of state
            MDP{n}.id.E(fg) = {iE};                  % parents of paths  

            MDP{n}.sA(gg)   = s;                     % parent stream
            MDP{n}.sB(fg)   = s;                     % parent stream
            MDP{n}.sC(gg)   = s;                     % child  stream

            % initial states and paths : outcomes at next time scale
            %--------------------------------------------------------------
            N(iD,:) = mdp.X(:,t);
            N(iE,:) = mdp.P(:,t);

        end
    end

    % Link streams
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

                    % prediction of inital states (D)
                    %------------------------------------------------------
                    gi = numel(MDP{n}.a) + 1;
                    gj = MDP{n - 1}.id.D{fj}(1);
                    A  = spm_dir_norm(MDP{n}.a{gj});
                    a  = zeros(Nj,Ni);
                    for f = 1:t
                        ii = sg{si}(i,f);
                        ij = sg{sj}(j,f);
                        a(:,ii) = a(:,ii) + A(:,ij);
                    end

                    % append mappings and id
                    %------------------------------------------------------
                    MDP{n}.a{gi}       = a;
                    MDP{n}.id.A{gi}    = fi;
                    MDP{n}.sA(end + 1) = si;
                    MDP{n}.sC(end + 1) = sj;
                    MDP{n - 1}.id.D{fj}(end + 1)  = gi;
                    MDP{n - 1}.ss.D{sj,sj}(fj,fj)  = gj;
                    MDP{n - 1}.ss.D{si,sj}(fi,fj)  = gi;
                    MDP{n - 1}.ss.ID{si,sj}(fi,fj) = spm_dir_MI(a);

                    % prediction of initial paths (E)
                    %------------------------------------------------------
                    gi = numel(MDP{n}.a) + 1;
                    gj = MDP{n - 1}.id.E{fj}(1);
                    A  = spm_dir_norm(MDP{n}.a{gj});
                    a  = zeros(Nu,Ni);
                    for f = 1:t
                        ii = sg{si}(i,f);
                        ij = sg{sj}(j,f);
                        a(:,ii) = a(:,ii) + A(:,ij);
                    end

                    % append mappings and iD
                    %------------------------------------------------------
                    MDP{n}.a{gi}       = a;
                    MDP{n}.id.A{gi}    = fi;
                    MDP{n}.sA(end + 1) = si;
                    MDP{n}.sC(end + 1) = sj;
                    MDP{n - 1}.id.E{fj}(end + 1)   = gi;
                    MDP{n - 1}.ss.E{sj,sj}(fj,fj)  = gj;
                    MDP{n - 1}.ss.E{si,sj}(fi,fj)  = gi;
                    MDP{n - 1}.ss.IE{si,sj}(fi,fj) = spm_dir_MI(a);


                end
            end
        end
    end

    % Termination criteria
    %======================================================================
    if true

        % if all streams comprise one group
        %------------------------------------------------------------------
        if max(Ng) < 2 && n > 1
            break
        end

    else

        % Combine streams in a final level
        %------------------------------------------------------------------
        if max(Ng) == 1
            S = [1 1 1 2*sum(Ng)];

            % Check there is only one group
            %--------------------------------------------------------------
            if (numel(Ng) == 1) || (Nt == 1)
                break
            end
        end
    end

    % remove unitary transitions
    %======================================================================
    if ~SPINBLOCK

        % find unitary likelihood mappings at next level
        %------------------------------------------------------------------
        Ns    = numel(MDP{n}.b);          % number of factors
        d     = true(1,Ns);

        % upper bound the number of generalised states
        %------------------------------------------------------------------
        for f = 1:Ns

            % for factors in leading stream
            %--------------------------------------------------------------
            if ismember(MDP{n}.sB(f),1)
                b    = sum(MDP{n}.b{f},3);

                % remove parents of non-recurring states
                %----------------------------------------------------------
                if max(sum(b,2)) < 1 && n < 2
                    d(f) = false;
                end

                % remove parents of itinerant states
                %----------------------------------------------------------
                if min(sum(b > 0,1)/size(b,1)) > 1/2 && numel(b) > 4
                    d(f) = false;
                end

                % remove parents of singletons
                %----------------------------------------------------------
                if isscalar(MDP{n}.b{f})
                    d(f) = false;
                end

            end
        end

        % compress tensors with no parents
        %------------------------------------------------------------------
        for f = find(ismember(d,false))

            % compress prior transitions
            %--------------------------------------------------------------
            MDP{n}.b{f} = sum(MDP{n}.b{f},'all');

            % and children
            %--------------------------------------------------------------
            for g = find(ismember([MDP{n}.id.A{:}],f))
                MDP{n}.a{g} = sum(MDP{n}.a{g},2);
            end
        end

        % streams at next level
        %------------------------------------------------------------------
        sB    = MDP{n}.sB(d);
        for s = 1:size(S,1)
            S(s,:) = [sum(ismember(sB,s)),1,1,2];
        end

        % and update parents of subordinate factors
        %------------------------------------------------------------------
        d = kron(d,true(1,2));
        i = find(d);
        for j = 1:numel(MDP{n}.id.D)
            MDP{n}.id.D{j} = find(ismember(i,MDP{n}.id.D{j}));
        end
        for j = 1:numel(MDP{n}.id.E)
            MDP{n}.id.E{j} = find(ismember(i,MDP{n}.id.E{j}));
        end

        % outcomes at next level
        %------------------------------------------------------------------
        O  = N(i,:);

    else

        % streams at next level
        %------------------------------------------------------------------
        for s = 1:size(S,1)
            S(s,:) = [size(MDP{n}.G{s},[1,2,3]),2];
        end

        % outcomes at next level
        %------------------------------------------------------------------
        O  = N;
    end
    
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

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[i,j] = spm_unique(O);

% reduction matrix
%--------------------------------------------------------------------------
R     = sparse(1:numel(j),j,1,numel(j),numel(i));

% Likelihood tensors
%--------------------------------------------------------------------------
Ng    = size(O,1);                          % number in group
a     = cell(Ng,1);
for g = 1:Ng
    a{g} = spm_cat(O(g,:))*R;
end

% Transition tensors
%--------------------------------------------------------------------------
Ns    = numel(i);                           % number of latent causes
Nt    = numel(j) - 1;
b     = zeros(Ns,Ns);
for t = 1:Nt

    % Is this an existing transition?
    %----------------------------------------------------------------------
    u  = find(b(j(t + 1),j(t),:),1,'first');
    if numel(u)

        % accumulate Dirichlet counts for this transition
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
        mdp.P{t} = logical(squeeze(b(j(t + 1),j(t),:)));
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


