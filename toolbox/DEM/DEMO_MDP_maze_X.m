function MDP = DEMO_MDP_maze_X
% Demo of mixed continuous and discrete state space modelling
%__________________________________________________________________________
%
% This demonstration of active inference focuses on navigation and
% planning. The idea is to demonstrate how epistemic foraging and goal
% (target) directed behaviour are integrated in the minimisation of
% expected free energy. In this illustration, and 8 x 8 maze is learned
% through novelty driven evidence accumulation - to learn the likelihood
% mapping between hidden states (locations in the maze) and outcomes
% (whether the current location is aversive or not). This accumulated
% experience is then used to plan a path from a start to an end (target
% location) under a task set specified by prior preferences over locations.
% These priors are based upon the distance between the current location and
% a target location.
%
% This version uses a belief propagation scheme (with deep policy searches)
% to illustrate prospective behaviour; namely, selecting policies or
% trajectories that transiently increased Bayesian risk. The code below can
% be edited to demonstrate different kinds of behaviour, under different
% preferences, policy depth and precisions.
%
% see also: DEM_MDP_maze.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEMO_MDP_maze_X.m 7766 2020-01-05 21:37:39Z karl $

% set up and preliminaries: first level
%--------------------------------------------------------------------------
rng('default')

% generative model at the sensory level (DEM): continuous states
%==========================================================================
% The generative model has two outcome modalities; namely, what (safe
% versus aversive) and where (the current location in a maze). These
% outcomes are generated from a single hidden factor (location), where the
% structure of the maze is encoded in the likelihood of observation mapping
% (that can be learned through experience). Allowable actions include four
% moves (up, down, left, right) and staying at the current location. These
% induce five transition matrices that play the role of empirical priors.
% Finally, prior preferences to avoid aversive locations while approaching
% a target location from an initial location (specified by START and END,
% respectively)
%--------------------------------------------------------------------------
label.factor     = {'where'};
label.modality   = {'what','where'};
label.outcome{1} = {'safe','aversive'};
label.action{1}  = {'up','down','left','right','stay'};

MAZE  = [...
    1 1 1 1 1 1 1 1;
    1 0 0 0 0 0 0 1;
    1 1 1 0 1 1 0 1;
    1 1 0 0 0 1 0 1;
    1 1 0 1 0 0 0 1;
    1 1 0 1 1 1 0 1;
    1 0 0 0 0 0 0 1;
    1 0 1 1 1 1 1 1];
END   = sub2ind(size(MAZE),5,5);                  % goal or target location
START = sub2ind(size(MAZE),8,2);                  % first or start location

% prior beliefs about initial states: D 
%--------------------------------------------------------------------------
D{1}  = zeros(numel(MAZE),1);
Ns    = numel(D{1});

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
A{1}  = [1 - MAZE(:), MAZE(:)]';                  % what
A{2}  = eye(Ns,Ns);                               % where
Ng    = numel(A);
for g = 1:Ng
    No(g)  = size(A{g},1);
end

% controlled transitions: B (up, down, left, right, stay)
%--------------------------------------------------------------------------
u    = [1 0; -1 0; 0 1; 0 -1; 0 0];               % allowable actions
nu   = size(u,1);                                 % number of actions
B{1} = zeros(Ns,Ns,nu);
[n,m] = size(MAZE);
for i = 1:n
    for j = 1:m
        
        % allowable transitions from state s to state ss
        %------------------------------------------------------------------
        s     = sub2ind([n,m],i,j);
        for k = 1:nu
            try
                ss = sub2ind([n,m],i + u(k,1),j + u(k,2));
                B{1}(ss,s,k) = 1;
            catch
                B{1}(s, s,k) = 1;
            end
        end
    end
end

% allowable actions: U
%--------------------------------------------------------------------------
U     = (1:nu)';

% priors: (negative cost) C: does not like what shocks and wants to be near
% the target location
%--------------------------------------------------------------------------
C{1}  = [2,-2];
[X,Y] = ind2sub(size(MAZE),END);
for i = 1:No(2)
    [x,y]   = ind2sub(size(MAZE),i);
    C{2}(i) = -sqrt((x - X)^2 + (y -Y)^2);
end

% basic MDP structure
%--------------------------------------------------------------------------
mdp.T = 8 + 1;                  % time horizon
mdp.N = 4;                      % policy depth
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model or likelihood
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states

mdp.label = label;
mdp       = spm_MDP_check(mdp);


% illustrate shortest path to target with suitable policy depth
%==========================================================================
% illustrate shortest path to target with suitable policy depth (e.g., N
% = 4), when the safe locations are known
%--------------------------------------------------------------------------
mdp.s = START;
for i = 1:4
    MDP   = mdp;
    MDP.N = i - 1;
    MDP   = spm_MDP_VB_XX(MDP);
    
    % cumulative reward
    %----------------------------------------------------------------------
    for j = 1:numel(C)
        R(i) = sum(MDP.C{j}(MDP.o(j,:),1));
    end
    
    % show results - behavioural
    %--------------------------------------------------------------------------
    str = sprintf('Figure %i',i);
    spm_figure('GetWin',str); clf
    spm_maze_plot(MDP,END)
    subplot(2,2,1), title(sprintf('Path: N = %i',i))
    
end

% show results - electrophysiological
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
spm_MDP_VB_LFP(MDP); subplot(3,2,5); delete(gca)

% evidence accumulation under task set: 
% removing knowledge about safe locations
%==========================================================================
clear MDP
mdp.a{1}   = ones(size(mdp.A{1}))/8;
mdp.a{2}   = mdp.A{2}*128;
[MDP(1:5)] = deal(mdp);
MDP = spm_MDP_VB_XX(MDP);

spm_figure('GetWin','Figure 6'); clf
spm_maze_plot(MDP,END)


% pure exploration
% removing preferences about proximity to target location
%==========================================================================
clear MDP
mdp.a{1}  = ones(size(mdp.A{1}))/64;
mdp.a{2}  = mdp.A{2}*128;
mdp.C{2}  = spm_zeros(C{2});

mdp.T = 2;
mdp.N = 2;
for i = 1:64
    
    % proceed with subsequent trial
    %----------------------------------------------------------------------
    MDP(i)   = spm_MDP_VB_XX(mdp);
    mdp      = MDP(i);
    mdp.s    = mdp.s(:,end);
    mdp.D{1} = MDP(i).X{1}(:,end);
    mdp.o    = [];
    mdp.u    = [];
    
end

% show results - behavioural
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7'); clf
spm_maze_plot(MDP)

return


function spm_maze_plot(MDP,END)
% illustrate  search graphically
%--------------------------------------------------------------------------
A  = spm_vec(MDP(1).A{1}(1,:));
ns = numel(A);
ni = sqrt(ns);
A  = reshape(A,ni,ni);
subplot(2,2,1), imagesc(A), axis image
title('Scanpath','fontsize',16);

% Cycle of the trials
%--------------------------------------------------------------------------
h     = [];
MS    = {};
MC    = {};
for p = 1:numel(MDP)
    
    %  current beliefs and preferences: A likelihood
    %----------------------------------------------------------------------
    if isfield(MDP,'a')
        Q = MDP(p).a{1};
    else
        Q = MDP(p).A{1};
    end
    Q     = Q/diag(sum(Q));
    Q     = Q(1,:);
    a     = reshape(Q(:),ni,ni);
    subplot(2,2,2), imagesc(a), axis image
    title('Likelihood','fontsize',16);
    
    %  current beliefs and preferences: B transitions
    %----------------------------------------------------------------------
    try
        b = diag(Q)*any(MDP(p).B{1},3);
    catch
        b = diag(Q)*any(MDP(p).B{1},3);
    end
    subplot(2,2,4), imagesc(-b), axis image
    title('Allowable transitions','fontsize',16);
    
    %  current beliefs and preferences: C preferences
    %----------------------------------------------------------------------
    C     = MDP(p).C{2}(:,1);
    C     = spm_softmax(C);
    C     = reshape(C,ni,ni);
    subplot(2,2,3), imagesc(C), axis image
    title('Preferences','fontsize',16);
    try
        [i,j] = ind2sub([ni,ni],MDP(p).s(1)); hold on
        plot(j,i,'.','MarkerSize',32,'Color','g');
        [i,j] = ind2sub([ni,ni],END);
        plot(j,i,'.','MarkerSize',32,'Color','r'); hold off
    end
    
    % cycle over  short-term searches
    %----------------------------------------------------------------------
    subplot(2,2,1),hold on
    s     = MDP(p).s;
    for t = 1:numel(s)
        
        % location
        %------------------------------------------------------------------
        [i,j] = ind2sub([ni,ni],s(t));
        h(end + 1) = plot(j,i,'.','MarkerSize',32,'Color','r');
        try
            set(h(end - 1),'Color','m','MarkerSize',16);
            j = [get(h(end - 1),'Xdata'), get(h(end),'Xdata')];
            i = [get(h(end - 1),'Ydata'), get(h(end),'Ydata')];
            plot(j,i,':r');
        end
        
        % save
        %------------------------------------------------------------------
        if numel(MS)
            MS(end + 1) = getframe(gca);
        else
            MS = getframe(gca);
        end
        
    end
    
    % save
    %----------------------------------------------------------------------
    subplot(2,2,3)
    if numel(MC)
        MC(end + 1) = getframe(gca);
    else
        MC = getframe(gca);
    end
    
end

% save movie
%--------------------------------------------------------------------------
subplot(2,2,1)
xlabel('click axis for movie')
set(gca,'Userdata',{MS,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

subplot(2,2,3)
xlabel('click axis for movie')
set(gca,'Userdata',{MC,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
