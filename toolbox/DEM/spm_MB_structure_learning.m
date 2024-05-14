function [MDP,RG,LG] = spm_MB_structure_learning(O,L,dt)
% RG structure learning of a hierarchical POMDP
% FORMAT [MDP,LG,RG] = spm_MB_structure_learning(O,L,dt)
% O{N,T} - probabilitic exemplars of paths (cell aray)
% L      - array of pixel locations
% dt     - time dilation [1 for images]
%
% MDP{n} - cell aray of MDPs
% RG{n}  - cell aray of group indices
% LG{n}  - cell aray of locations
%
% This routine returns a hierarchy of generative models (MDP{n}), given a
% sequence of outcomes using RG (renormalization group) operators.
%
% see: spm_MDP_structure_learning.m
%      spm_MDP_structure_teaching.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________

% options for model inversion (and evaluation)
%==========================================================================
if nargin < 3, dt = 2; end              % time scaling
O     = {O};                            % outcomes

% learn the dynamics in the form of a hierarchical MDP
%==========================================================================
MDP   = {};
RG    = {};
LG    = {};
for n = 1:8

    % locations
    %----------------------------------------------------------------------
    if size(L,2) > 1
        hold on, plot(L(:,2),L(:,1),'.r','MarkerSize',2*n + 2), axis ij image
        hold on, plot(L(:,2),L(:,1),'ow','MarkerSize',2*n + 2), axis ij image
        title(sprintf('Locations at scale %i',n)), drawnow
    end

    % Grouping into a partition of outcomes
    %----------------------------------------------------------------------
    G     = spm_space(L);
    T     = spm_time(size(O{n},2),dt);
    for g = 1:numel(G)

        % structure_learning from unique exemplars
        %------------------------------------------------------------------
        mdp  = spm_structure_fast(O{n}(G{g},:),dt);

        % place in hierarchical structure
        %------------------------------------------------------------------
        MDP{n}.A(G{g},1) = mdp.a;
        MDP{n}.B(g,1)    = mdp.b;
        for j = 1:numel(G{g})
            MDP{n}.id.A{G{g}(j)} = uint16(g);
        end
    end

    % place indices and locations in MDP
    %----------------------------------------------------------------------
    if nargout > 1
        RG{n} = G;
        LG{n} = L;
    end
    MDP{n}.RG = G;
    MDP{n}.LG = L;

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
        ig    = 1;
        L     = zeros(1,size(L,2));
        for g = 1:numel(G)

            % states, paths and average location for this goup
            %--------------------------------------------------------------
            qs = pdp.X{g}(:,1);
            qu = pdp.P{g}(:,end);
            ml = mean(MDP{n}.LG(G{g},:));

            % states (odd)
            %--------------------------------------------------------------
            MDP{n}.id.D{g} = ig;
            O{n + 1}{ig,t} = qs;
            L(ig,:)        = ml;
            ig = ig + 1;

            % paths (even)
            %--------------------------------------------------------------
            MDP{n}.id.E{g} = ig;
            O{n + 1}{ig,t} = qu;
            L(ig,:)        = ml;
            ig = ig + 1;

        end
    end

    % Check whether there is only one group or (generalised) object
    %----------------------------------------------------------------------
    if numel(G) == 1
        break
    end

end

return

function mdp = spm_structure_fast(O,dt)
% A fast form of structure learning
% FORMAT mdp = spm_structure_fast(O,dt)
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

% defaults
%--------------------------------------------------------------------------
if nargin < 2, dt = 2; end

% Unique outputs
%--------------------------------------------------------------------------
j     = spm_unique(O);

% Likelihood tensors
%--------------------------------------------------------------------------
Ng    = size(O,1);                          % number in group
Ns    = numel(unique(j));                   % number of latent causes
a     = cell(Ng,1);
for g = 1:Ng
    for s = 1:Ns
        a{g}(:,s) = full(mean(spm_cat(O(g,ismember(j,s))),2));
    end
end

% return if no dynamics
%--------------------------------------------------------------------------
if dt < 2
    mdp.a    = a;
    mdp.b{1} = eye(Ns,Ns);
    return
end

% Transition tensors
%--------------------------------------------------------------------------
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
for u = 1:size(b,3)
    for s = 1:Ns
        if ~any(b(:,s,u))
            i = find(any(squeeze(b(:,s,:)),2),1);
            b(i,s,u) = 1;
        end
    end
end

% Vectorise cell array of likelihood tensors and place in structure
%--------------------------------------------------------------------------
mdp.a    = a;
mdp.b{1} = b;

return

function G = spm_space(L)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT G = spm_space(L)
%--------------------------------------------------------------------------
% L  - Location array
%
% g  - Cell array of outcome indices
%
% Effectively, this leverages the conditional independencies that inherit
% from local interactions; of the kind found in metric spaces that preclude
% action at a distance.
%
% The implicit grouping identifies a reduced number of group centroids and
% assigns lower scale locations to the nearest centroid. The group averages
% then constitute the locations for the next scale. In this example, the
% number of locations is reduced by a factor of two in both dimensions.
%--------------------------------------------------------------------------
if isvector(L)

    % locations
    %----------------------------------------------------------------------
    L     = L(:);
    Nl    = size(L,1);
    Ng    = size(unique(L,'rows'),1);
    Ng    = ceil(Ng/2);
    if Ng == 1
        G = {1:Nl};
        return
    end

    % decimate locations
    %----------------------------------------------------------------------
    x     = linspace(min(L(:,1)),max(L(:,1)),Ng);

    % nearest (reduced) location
    %----------------------------------------------------------------------
    for i = 1:Nl
        [~,j] = min(minus(x,L(i)).^2);
        g(i)  = j;
    end

    % grouping partition
    %----------------------------------------------------------------------
    u     = unique(g,'stable');
    for i = 1:numel(u)
        G{i}  = find(ismember(g,u(i)));
    end

else

    % locations
    %----------------------------------------------------------------------
    Nl    = size(L,1);
    Ng    = size(unique(L,'rows'),1);
    Ng    = ceil(sqrt(Ng/4));
    if Ng == 1
        G = {1:Nl};
        return
    end

    % decimate locations
    %----------------------------------------------------------------------
    x     = linspace(min(L(:,1)),max(L(:,1)),Ng);
    y     = linspace(min(L(:,2)),max(L(:,2)),Ng);
    for n = 1:prod([Ng,Ng])
        i = spm_index([Ng,Ng],n);
        R(n,:) = [x(i(1)),y(i(2))];
    end

    % nearest (reduced) location
    %----------------------------------------------------------------------
    for i = 1:Nl
        [~,j] = min(sum(minus(R,L(i,:)).^2,2));
        g(i)  = j;
    end

    % grouping partition
    %----------------------------------------------------------------------
    u     = unique(g,'stable');
    for i = 1:numel(u)
        G{i} = find(ismember(g,u(i)));
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


function [j] = spm_unique(a)
% information geometry of a likelihood mapping
% FORMAT [j] = spm_unique(a)
% a{g}    - Dirichlet tensor for modality g
%
% This routine uses the information geometry inherent in a likelihood
% mapping. It computes a similarity matrix based upon an approximation to
% the information length between columns of a likelihood tensor (as
% approximated by the KL divergence of the implicit categorical
% distributions).
%__________________________________________________________________________
% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging

% Fast approximation by simply identifying unique locations in a
% multinomial statistical manifold, after discretising to probabilities of
% zero, half and one (using Matlab’s unique and fix operators).
%--------------------------------------------------------------------------
o       = spm_cat(a)';
[~,~,j] = unique(fix(2*o),'rows','stable');

return

% information geometry – divergence : likelihood mapping
%--------------------------------------------------------------------------
C     = 0;
for g = 1:size(a,1)
    q  = spm_cat(a(g,:));
    p  = spm_log(q);
    KL = q'*p;
    KL = minus(KL,diag(KL));
    C  = C + (KL + KL').^2;
end

% similarity matrix (assuming normalised vectors on a hypersphere)
%--------------------------------------------------------------------------
C     = real(C);                                % distance metric
c     = C/max(C,[],'all');                      % normalise distance
c     = 1 - 2*c;                                % correlation matrix
[u,s] = spm_svd(c,1/256);                       % eigenvectors
o     = u*s;

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[~,~,j] = unique(fix(o),'rows','stable');


% cut and paste for graphics
%==========================================================================

% Display [eigen] space
%==========================================================================
spm_figure('GetWin','Information geometry'); clf;

s     = diag(real(s));
[s,i] = sort(s,'descend');
u     = real(u(:,i));

% eigenspace
%--------------------------------------------------------------------------
subplot(2,2,2), bar(s(1:16)),
xlabel('eigenvectors'), ylabel('eigenvalues'), title('Eigenvalues')
axis square

subplot(2,2,1), imagesc(c);
xlabel('latent states'), ylabel('latent states')
title('Correlation matrix'), axis image

subplot(2,1,2)
plot3(u(:,1),u(:,2),u(:,3),'.'), title('Embedding sapce')
axis image, grid on

return