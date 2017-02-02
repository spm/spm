function MDP = DEMO_MDP_maze
% Demo of mixed continuous and discrete state space modelling
%__________________________________________________________________________
%
% this demonstration of active inference focuses on navigation and planning
% in a fairly complicated maze. The idea is to demonstrate how epistemic
% foraging and goal (target) directed behaviour are integrated in the
% minimisation of expected free energy. In this illustration, and 8 x 8
% maze is learned through novelty driven evidence accumulation – to learn
% the likelihood mapping between hidden states (locations in the maze) and
% outcomes (whether the current location is open or closed). This
% accumulated experience is then used to plan a path from a start to an end
% (target location) under a task set specified by prior preferences over
% locations. These priors are based upon a simple diffusion (CF backwards
% induction) heuristic that specifies subgoals. The subgoals (i.e.,
% locations) contain the most paths from the target within the horizon of
% the current policy.
%
% we will first illustrate the novelty driven epistemic foraging that
% efficiently scans the maze to learn its structure. We then simulate
% planning of(shortest path) trajectory to the targetunder the assumption
% the maze has been previously learned. Finally, we consider exploration
% under prior preferences to simulate bbehaviour when both epistemic and
% goal directed imperatives are in play. The focus on this demo is on
% behavioural responses and electrophysiological responses over moves.
%
% a key aspect of this formulation is the  hierarchical decomposition of
% goal directed behaviour into subgoals that are within the horizon of a
% limited policy – here, to moves that correspond to a trial. The prior
% preferences then contextualise each policy or trial to ensure that the
% ultimate goal is achieved.
%
% see also: spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEMO_MDP_maze.m 7003 2017-02-02 18:22:56Z karl $

% set up and preliminaries: first level
%--------------------------------------------------------------------------
rng('default')

% generative model at the sensory level (DEM): continuous states
%==========================================================================
% the generative model has two outcome modalities; namely, what (open
% versus closed) and where (the current location in a maze). These outcomes
% are generated from a single hidden factor (location), where the structure
% of the maze is encoded in the likelihood of observation mapping (that can
% be learned through experience). Allowable actions  include for moves (up,
% down, left, right) and staying at the current location. These induce five
% transition matrices that play the role of empirical priors. Finally,
% prior preferences are based upon allowable transitions (that are function
% of learned accumulated likelihood), which are used to define attractive
% locations within the horizon of two-move policies. These priors implement
% a task set and are returned by a subfunction below: spm_maze_cost
%--------------------------------------------------------------------------
label.factor     = {'where'};
label.modality   = {'what','where'};
label.outcome{1} = {'open','closed'};

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

% allowable policies (2 moves): V
%--------------------------------------------------------------------------
V     = [];
for i = 1:nu
    for j = 1:nu
        V(:,end + 1) = [i;j];
    end
end

% priors: (negative cost) C:
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g} = zeros(No(g),1);
end

% basic MDP structure
%--------------------------------------------------------------------------
mdp.V = V;                      % allowable policies
mdp.A = A;                      % observation model or likelihood
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states

mdp.label = label;
mdp.alpha = 64;
mdp.tau   = 8;
mdp       = spm_MDP_check(mdp);


% exploratory (32 trial) sequence (no experience or task set)
%==========================================================================
% These simulations use a subroutine (spm_maze_search) to perform recursive
% variational inversions of the active inference scheme under different
% levels of experience and prior preferences
%--------------------------------------------------------------------------
SDP = spm_maze_search(mdp,32,START,END,0,0);

% show results - behavioural
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_maze_plot(SDP,END)

% show results - electrophysiological
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(SDP);

% exploratory (8 trial) sequence (with experience and task set)
%==========================================================================
MDP = spm_maze_search(mdp,8,START,END,128,1);

% show results in terms of path
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_maze_plot(MDP,END)

% and implicit subgoals
%--------------------------------------------------------------------------
for i = 1:8
    subplot(8,2,(i - 1)*2 + 2);
    c = softmax(MDP(i).C{2}(:,1));
    imagesc(spm_unvec(c,MAZE))
    axis image off, title(sprintf('%s %i','Trial',i))
end

% evidence accumulation under task set
%==========================================================================
MDP   = mdp;
for i = 1:4
    
    % return to start
    %----------------------------------------------------------------------
    MDP = spm_maze_search(MDP(end),8,START,END,0,1);
    
    % show results in terms of path
    %----------------------------------------------------------------------
    spm_figure('GetWin',sprintf('%s %i','Search',i)); clf
    spm_maze_plot(MDP,END)
    
    % and electrophysiological responses
    %----------------------------------------------------------------------
    spm_figure('GetWin',sprintf('%s %i','Responses',i)); clf
    spm_MDP_VB_game(MDP(1:8));
    subplot(6,1,3), set(gca,'YLim',[-1 1]/3);
    subplot(6,1,4), set(gca,'YLim',[0 1]);
    
end

% in silico psychophysical experiment
%==========================================================================
% Here, we return to the exploratory simulation above and probe each level
% of experience by asking the subject to execute a path to target. The
% behaviour is then assessed in terms of the latency with which the target
% will go is required  – and the number of mistakes or exploratory
% excursions into  closed locations.
%--------------------------------------------------------------------------
N     = [];
M     = [];
T     = [];
for i = 1:2:numel(SDP)
    
    % execute path to target with increasing experience
    %----------------------------------------------------------------------
    MDP = spm_maze_search(SDP(i),8,START,END,0,1);
    
    % record behaviour
    %----------------------------------------------------------------------
    s   = spm_cat({MDP.s}); s(3:3:end) = [];
    s   = [s,END];
    N   = [N,find(s == END,1)];                      % latency
    M   = [M,sum(MAZE(s))];                          % mistakes
    T   = [T,i/2];                                   % exposure to maze
    
end

% show results
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
subplot(2,1,1), bar([M;N]'), xlabel('Exposure (seconds)'), axis square 
title('Performance','FontSize',16),legend({'Mistakes','Latency'})

return



function MDP = spm_maze_search(mdp,N,START,END,alpha,beta)
% FORMAT MDP = spm_maze_search(mdp,N,START,END,alpha,beta)
% mdp   - MDP structure
% N     - number of trials (i.e., policies: default 8)
% START - index of intial state (default 1)
% END   - index of target state (default 1)
% alpha - prior concentration parameter for likelihood (default 128)
% beta  - precision of prior preference (default 0)
%
% MDP   - MDP structure array

% preliminaries
%--------------------------------------------------------------------------
try, N;     catch, N     = 8;   end
try, START; catch, START = 1;   end
try, END;   catch, END   = 1;   end
try, alpha; catch, alpha = 128; end
try, beta;  catch, beta  = 0;   end

% initialise concentration parameters: a (if unspecified)
%--------------------------------------------------------------------------
if ~isfield(mdp,'a')
    mdp.a{1} = ones(size(mdp.A{1}))/8 + mdp.A{1}*alpha;
    mdp.a{2} = mdp.A{2}*128;
end
if ~isfield(mdp,'o')
    mdp.o = [];
end
if ~isfield(mdp,'u')
    mdp.u = [];
end
mdp.s = START;

% Evaluate a sequence of eye movements
%==========================================================================

for i = 1:N
    
    % Evaluate preferred states (subgoals) on the basis of current beliefs
    %----------------------------------------------------------------------
    mdp.C{2} = spm_maze_cost(mdp,END)*beta;
    
    % proceed with subsequent trial
    %----------------------------------------------------------------------
    MDP(i)   = spm_MDP_VB_X(mdp);
    mdp      = MDP(i);
    mdp.s    = mdp.s(:,end);
    mdp.D{1} = MDP(i).X{1}(:,end);
    mdp.o    = [];
    mdp.u    = [];
    
end

return


function C = spm_maze_cost(MDP,END)
% Evaluate subgoals using
%==========================================================================
START = MDP.s(1);
if isfield(MDP,'a')
    Q = MDP.a{1};
else
    Q = MDP.A{1};
end
Q   = Q/diag(sum(Q));
Q   = Q(1,:);                                % open states
P   = diag(Q)*any(MDP.B{1},3);               % allowable transitions
ns  = length(Q);                             % number of states
X   = zeros(ns,1);X(START) = 1;              % initial state
Y   = zeros(ns,1);Y(END)   = 1;              % target state


% Preclude transitions to closed states and evaluate graph Laplacian
%--------------------------------------------------------------------------
P   = P - diag(diag(P));
P   = P - diag(sum(P));
P   = expm(P);

% evaluate (negative) cost as a path integral conjunctions
%--------------------------------------------------------------------------
for t = 1:size(MDP.V,1)
    X = P*X;
end
X     = X > exp(-3);
C     = log(X.*(P*Y) + exp(-32));

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
