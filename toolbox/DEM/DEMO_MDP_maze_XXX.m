function MDP = DEMO_MDP_maze_XXX
% Demo of inductive inference based on a learned likelihood mapping
%__________________________________________________________________________
%
% This demonstration of active inference focuses on navigation and
% planning. The idea is to demonstrate how epistemic foraging and goal
% (target) directed behaviour are integrated in the minimisation of
% expected free energy. In this illustration, and 10 x 10 maze is learned
% through novelty-driven evidence accumulation - to learn the likelihood
% mapping between hidden states (locations in the maze) and outcomes
% (whether the current location is aversive or not). This accumulated
% experience is then used to plan a path from a start to an end (target
% location) under a task set specified by prior preferences over locations.
%
% This version uses a belief propagation scheme (with deep policy searches)
% to illustrate intentional behaviour. This search is constrained using
% inductive inference; namely, placing priors over the subsequent states
% that preclude regimes of state space that cannot access the final state.
% These priors rest upon a time reversed protocol under allowable (learned)
% state transitions. In effect, priors over the final state contextualise
% priors over outcomes by precluding states that lie ahead, on a path of
% least action that minimises expected free energy.The code below can be
% edited to demonstrate different kinds of behaviour, under different
% preferences, policy depth and precisions.
%
% see also: spm_MPD_VB_XXX.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up and preliminaries: first level
%--------------------------------------------------------------------------
rng('default')

% generative model at the sensory level (DEM): continuous states
%==========================================================================
% The generative model has two outcome modalities; namely, what (safe
% versus aversive) and where (the current location in a maze). These
% outcomes are generated from a single hidden factor (location), where the
% structure of the maze is encoded in the likelihood mapping (that can be
% learned through experience). Allowable actions include four moves (up,
% down, left, right) and staying at the current location. These induce five
% transition matrices that play the role of empirical priors. Finally,
% prior constraints specify aversive locations while seeking a target
% location from an initial location (i.e., START and END, respectively)
%--------------------------------------------------------------------------
label.factor     = {'where'};
label.modality   = {'what','where'};
label.outcome{1} = {'safe','aversive'};
label.action{1}  = {'up','down','left','right','stay'};

MAZE  = [...
    1 1 1 1 1 1 1 1 1 1;
    1 0 0 0 1 0 0 0 0 1;
    1 1 1 0 1 1 0 1 0 1;
    1 0 0 0 0 1 0 1 0 1;
    1 1 1 0 1 1 1 1 0 1;
    1 1 0 0 0 1 0 1 0 1;
    1 1 0 1 0 0 0 1 0 1;
    1 1 0 1 1 1 0 0 0 1;
    1 0 0 0 0 0 0 1 0 1;
    1 0 1 1 1 1 1 1 1 1];
END   = sub2ind(size(MAZE),4,7);                  % goal or target location
START = sub2ind(size(MAZE),10,2);                 % first or start location

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

% Constraints: does not like shocks (C) and wants to be at END (H)
%--------------------------------------------------------------------------
C{1}  = spm_softmax([1;0]);
C{2}  = spm_softmax(ones(No(2),1));   % path   dependent preferences
H{1}  = ones(No(2),1);                % path independent preferences

H{1}(END) = 32;                       % Target or end state

% basic MDP structure
%--------------------------------------------------------------------------
mdp.N = 0;                            % policy depth
mdp.U = U;                            % allowable policies
mdp.A = A;                            % observation model or likelihood
mdp.B = B;                            % transition probabilities
mdp.C = C;                            % preferred outcomes
mdp.D = D;                            % prior over initial states

mdp.label = label;
mdp.p     = 1/32;

% Likelihood learning: removing preferences about target location
%==========================================================================
MDP       = mdp;
MDP.a{1}  = spm_zeros(mdp.A{1}) + mdp.p;
MDP.a{2}  = mdp.A{2}*512;
MDP.s     = START;
MDP.T     = 256;

% Active learning
%----------------------------------------------------------------------
PDP       = spm_MDP_VB_XXX(MDP);
MDP       = spm_MDP_VB_update(MDP,PDP);

% show results - behavioural
%--------------------------------------------------------------------------
spm_figure('GetWin','Learning'); clf
spm_maze_plot(PDP)


% Paths of least action
%==========================================================================
% illustrate shortest path to target with suitable policy depth, Under
% greater and lesser precision over (path independent) constraints C. With
% imprecise constraints, the agent still likes to explore a bit
%--------------------------------------------------------------------------
MDP.H = H;
MDP.s = START;
MDP.T = 32;
for i = [1,8]

    % increase precision of constraints
    %----------------------------------------------------------------------
    MDP.C{1} = spm_softmax([1;0]*i);
    PDP      = spm_MDP_VB_XXX(MDP);
    
    % show results - behavioural
    %----------------------------------------------------------------------
    str = sprintf('Figure %i',i);
    spm_figure('GetWin',str); clf
    spm_maze_plot(PDP,END)
    subplot(2,2,1), title(sprintf('Path: C = %i',i))
    
end


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
for p = 1:numel(MDP)
    
    %  current beliefs and preferences: A likelihood
    %----------------------------------------------------------------------
    if isfield(MDP,'a')
        Q = MDP(p).a{1};
    else
        Q = MDP(p).A{1};
    end
    Q     = spm_dir_norm(Q);
    Q     = Q(1,:);
    a     = reshape(Q(:),ni,ni);
    subplot(2,2,2), imagesc(a), axis image
    title('Likelihood','fontsize',16);
    
    %  current beliefs and preferences: B transitions
    %----------------------------------------------------------------------
    try
        b = sum(MDP(p).b{1},3);
    catch
        b = sum(MDP(p).B{1},3);
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
    
    % cycle over short-term searches
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
    end
end

return
