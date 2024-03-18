function [MDP,RG,LG] = spm_MB_structure_learning(O,L)
% RG structure learning of a hierarchical POMDP
% FORMAT [MDP,LG,RG] = spm_MB_structure_learning(O,L)
% O{N,T} - probabilitic exemplars of paths (cell aray)
% L      - array of pixel locations
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
dt    = 2;                              % time scaling
O     = {O};                            % outcomes

% learn the dynamics in the form of a hierarchical MDP
%==========================================================================
MDP   = {};
RG    = {};
LG    = {};
for n = 1:8

    % locations
    %----------------------------------------------------------------------
    if true
        hold on, plot(L(:,2),L(:,1),'.r'), axis ij image
        hold on, plot(L(:,2),L(:,1),'ow'), axis ij image
        title(sprintf('Locations at scale %i',n)), drawnow
    end

    % Grouping into a partition of outcomes
    %----------------------------------------------------------------------
    G     = spm_tile(L);
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
        for g = 1:numel(G)

            % check for singletons
            %------------------------------------------------------------------
            if numel(MDP{n}.B{g}) > 1
                MDP{n}.id.D{g} = 2*g - 1;
                MDP{n}.id.E{g} = 2*g - 0;
                O{n + 1}{MDP{n}.id.D{g},t} = pdp.X{g}(:,1);
                O{n + 1}{MDP{n}.id.E{g},t} = pdp.P{g}(:,end);
            end
        end
    end

    % average location of next level outcomes
    %----------------------------------------------------------------------
    L     = zeros(size(O{n + 1},1),size(L,2));
    for g = 1:numel(G)
        L(2*g - 1,:) = mean(MDP{n}.LG(G{g},:));
        L(2*g - 0,:) = L(2*g - 1,:);
    end

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

function G = spm_tile(L)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT G = spm_tile(L)
%--------------------------------------------------------------------------
% L  - Location array
%
% g  - Cell array of outcome indices
%
% Effectively, this leverages the conditional independencies that inherit
% from local interactions; of the kind found in metric spaces that preclude
% action at a distance.
%--------------------------------------------------------------------------

% locations
%--------------------------------------------------------------------------
Nl    = size(L,1);
Ng    = size(unique(L,'rows'),1);
Ng    = ceil(sqrt(Ng/4));
if Ng == 1
    G = {1:Nl};
    return
end

% decimate locations
%--------------------------------------------------------------------------
x     = linspace(min(L(:,1)),max(L(:,1)),Ng);
y     = linspace(min(L(:,2)),max(L(:,2)),Ng);
for n = 1:prod([Ng,Ng])
    i = spm_index([Ng,Ng],n);
    R(n,:) = [x(i(1)),y(i(2))];
end

% nearest (reduced) location
%--------------------------------------------------------------------------
for i = 1:Nl
    [~,j] = min(sum(minus(R,L(i,:)).^2,2));
    g(i)  = j;
end

% grouping partition
%--------------------------------------------------------------------------
u     = unique(g,'stable');
for i = 1:numel(u)
    G{i} = find(ismember(g,u(i)));
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

