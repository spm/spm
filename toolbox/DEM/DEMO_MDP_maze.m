function MDP = DEMO_MDP_maze
% Demo of mixed continuous and discrete state space modelling
%__________________________________________________________________________
%
% This routine  illustrates the combination of discrete and continuous
% state space models for active inference. In this example, the lowest
% level of a hierarchical Markov decision process (used to illustrate
% evidence accumulation during reading in related simulations) is equipped
% with a continuous time and state space dynamical model at the lowest
% level. This allows one to model both the categorical belief updates
% using belief propagation and the continuous belief updates using
% Bayesian filtering within the same model and associated inversion
% scheme.
%
% More details about each level of the model are provided in line as
% annotated descriptions.
%
% see also: spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEMO_MDP_maze.m 6977 2016-12-24 17:48:44Z karl $

% set up and preliminaries: first level
%--------------------------------------------------------------------------
rng('default')

% generative model at the sensory level (DEM): continuous states
%==========================================================================
% This level of model specification concerns the sampling of continuous
% data; here, visual stimuli encoded in a global structure (STIM). This is
% a generic specification that allows one to place various stimuli in the
% visual field – with user-specified parameters for sampling. The hidden
% causes of this generative model correspond to one attracting location
% and the content or stimulus that will be sampled at that location. The
% dynamics or hidden states in this level of the model are simple: the
% attracting location simply attracts the point of foveal fixation.
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
END   = sub2ind(size(MAZE),5,7);
START = sub2ind(size(MAZE),8,2);

% prior beliefs about initial states (in terms of counts): D and d
%--------------------------------------------------------------------------
D{1}  = zeros(numel(MAZE),1);
D{1}(START) = 1;
Ns    = numel(D{1});


% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
A{1}  = [1 - MAZE(:), MAZE(:)]';                  % what
A{2}  = eye(Ns,Ns);                               % where
Ng    = numel(A);
for g = 1:Ng
    No(g)  = size(A{g},1);
end

% controlled transitions:r
%--------------------------------------------------------------------------
B{1} = zeros(Ns,Ns,5);
[n,m] = size(MAZE);
for i = 1:n
    for j = 1:m
        
        % state - where
        %------------------------------------------------------------------
        s = (j - 1)*n + i;
        
        % up
        %------------------------------------------------------------------
        ii = i - 1; jj = j; ss = (jj - 1)*n + ii;
        try, MAZE(ii,jj); B{1}(ss,s,1) = 1; catch, B{1}(s,s,1) = 1; end
        % down
        %------------------------------------------------------------------
        ii = i + 1; jj = j; ss = (jj - 1)*n + ii;
        try, MAZE(ii,jj); B{1}(ss,s,2) = 1; catch, B{1}(s,s,2) = 1; end
        % left
        %------------------------------------------------------------------
        ii = i; jj = j - 1; ss = (jj - 1)*n + ii;
        try, MAZE(ii,jj); B{1}(ss,s,3) = 1; catch, B{1}(s,s,3) = 1; end
        % right
        %------------------------------------------------------------------
        ii = i; jj = j + 1; ss = (jj - 1)*n + ii;
        try, MAZE(ii,jj); B{1}(ss,s,4) = 1; catch, B{1}(s,s,4) = 1; end
        
        
    end
end

% stay
%------------------------------------------------------------------
B{1}(:,:,5)  = eye(Ns,Ns);


% allowable policies (2 moves) V
%--------------------------------------------------------------------------
V     = [];
for i = 1:5
    for j = 1:5
        V(:,end + 1) = [i;j];
    end
end

% priors: (utility) C: the agent does not like 'closed' outcomes
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g} = zeros(No(g),1);
end

% MDP structure for this level (and subordinate MDP level)
% including links from outcomes at the current level to states of level
% below. This  complete the specification of the mixed hierarchical model
%--------------------------------------------------------------------------
mdp.V = V;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.s = START;                  % initial state

mdp.label = label;
mdp.alpha = 64;
mdp       = spm_MDP_check(mdp);


% invert this model
%==========================================================================

% exploratory sequence
%--------------------------------------------------------------------------
alpha = 2;
for i = 1:6
    
    % Evaluate preferred states (subgoals) on the basis of current beliefs
    %----------------------------------------------------------------------
    mdp.C{2} = spm_maze_cost(mdp,END)*alpha;
    
    % proceed with subsequent trial
    %----------------------------------------------------------------------
    MDP(i)   = spm_MDP_VB_X(mdp);
    mdp      = MDP(i);
    mdp.s    = mdp.s(:,end);
    mdp.D{1} = MDP(i).X{1}(:,end);
    mdp.o    = [];
    mdp.u    = [];
    
end


% show belief updates (and behaviour) over trials
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1));

% illustrate phase-precession and electrophysiological responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP(1),[],1);

% illustrate sequence
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_game(MDP);

spm_figure('GetWin','Figure 4'); clf
spm_maze_plot(MDP,END)

return


function C = spm_maze_cost(MDP,END)
% Evaluate subgoals using 
%==========================================================================
Q   = MDP.A{1}(1,:);                         % open states
P   = diag(Q)*any(MDP.B{1},3);               % allowable transitions
ns  = length(Q);                             % number of states
X   = zeros(ns,1);X(MDP.s(1)) = 1;           % initial state
Y   = zeros(ns,1);Y(END)      = 1;           % targets state

% Preclude transitions to closed states and evaluate graph Laplacian
%--------------------------------------------------------------------------
P   = P - diag(diag(P));
P   = P - diag(sum(P));
P   = expm(P);

% evaluate (negative) cost as a path integral conjunctions
%--------------------------------------------------------------------------
C   = 0;
for t = 1:16
    X = P*X;
    Y = P*Y;
    C = C + log(X + exp(-16)) + log(Y + exp(-16));
end

C(END) = C(END) + 4;

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
    try
        a = spm_vec(MDP(p).a{1}(1,:));
    catch
        a = spm_vec(MDP(p).A{1}(1,:));
    end
    a     = reshape(a,ni,ni);
    subplot(2,2,2), imagesc(a), axis image
    title('Likelihood','fontsize',16);
    
    %  current beliefs and preferences: B transitions
    %----------------------------------------------------------------------
    try
        b = diag(MDP(p).a{1}(1,:))*any(MDP(p).B{1},3);
    catch
        b = diag(MDP(p).A{1}(1,:))*any(MDP(p).B{1},3);
    end
    subplot(2,2,4), imagesc(-b), axis image
    title('Empirical priors','fontsize',16);
    
    %  current beliefs and preferences: C preferences
    %----------------------------------------------------------------------
    C     = MDP(p).C{2}(:,1);
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
        end
        
        % save
        %----------------------------------------------------------------------
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
set(gca,'Userdata',{MS,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

subplot(2,2,3)
set(gca,'Userdata',{MC,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
