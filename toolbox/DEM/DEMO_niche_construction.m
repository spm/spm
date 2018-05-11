function DEMO_niche_construction
% Demo of niche construction using active inference
%__________________________________________________________________________
%
% The free-energy principle is an attempt to explain the structure of the
% agent and its brain, starting from the fact that an agent exists (Friston
% and Stephan, 2007; Friston et al., 2010). More specifically, it can be
% regarded as a systematic attempt to understand the ‘fit’ between an
% embodied agent and its niche, where the quantity of free-energy is a
% measure for the ‘misfit’ or disattunement (Bruineberg and Rietveld, 2014)
% between agent and environment. This paper offers a proof-of-principle
% simulation of niche construction under the free-energy principle. The key
% point of this paper is that the minimum of free-energy is not at a point
% in which the agent is maximally adapted to the statistics of a static
% environment, but can better be conceptualized an attracting manifold
% within the joint  agent-environment state-space as a whole, which the
% system tends toward through mutual interaction. We will provide a general
% introduction to active inference and the free-energy principle. Using
% Markov Decision Processes (MDPs), we then describe a canonical generative
% model and the ensuing update equations that minimize free-energy. We then
% apply these equations to simulations of foraging in an environment; in
% which an agent learns the most efficient path to a pre-specified
% location. In some of those simulations, unbeknownst to the agent, the
% environment changes as a function of the activity of the agent (i.e.
% unintentional niche construction occurs). We will show how, depending on
% the relative inertia of the environment and agent, the joint
% agent-environment system moves to different attracting sets of jointly
% minimized free-energy.
%
% see also: spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEMO_niche_construction.m 7310 2018-05-11 19:24:09Z karl $


%% Set up and preliminaries
%--------------------------------------------------------------------------
rng('default')

label.factor     = {'where'};
label.modality   = {'what','where'};
label.outcome{1} = {'open','closed'};
MAZE =  [...
    0 0 0 0 0 0 0 0;
    1 1 0 1 1 1 1 0;
    1 1 0 1 1 1 1 0;
    1 1 0 1 1 1 1 0;
    1 1 0 1 1 1 1 0;
    0 0 0 1 1 1 1 0;
    0 1 1 1 1 1 1 0;
    0 0 0 0 0 0 0 0];

END   = sub2ind(size(MAZE),1,1);                  % goal or target location
START = sub2ind(size(MAZE),8,1);

% prior beliefs about initial states: D
%--------------------------------------------------------------------------
D{1}  = zeros(numel(MAZE),1);
Ns    = numel(D{1});

% Capture generative process in form of A: Aaux
%--------------------------------------------------------------------------
Aaux{1} = [1 - MAZE(:), MAZE(:)]';                % what
Aaux{2} = eye(Ns,Ns);                             % where
Ng    = numel(Aaux);
for g = 1:Ng
    No(g)  = size(Aaux{g},1);
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

% priors: (negative cost) C:
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g} = zeros(No(g),1);
end

% allowable policies (1 move): V
%--------------------------------------------------------------------------
V = 1:nu;

% basic MDP structure
%--------------------------------------------------------------------------
mdp.V = V;                      % allowable policies
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states

mdp.label = label;
mdp.alpha = 8;
mdp.tau   = 8;


%% Simulation 1: Learning over trials
%==========================================================================
% Define concentration parameters for environment:
%--------------------------------------------------------------------------
Acp.open   = [4;0];
Acp.closed = [0;4];

% Define concentration parameters for agent:
%--------------------------------------------------------------------------
acp     = [1/8;1/8];

% Generative process: Generate probabilistic mapping from hidden states to
% outcomes: A and matrix of concentration parameters
%--------------------------------------------------------------------------
[A, As] = spm_gen_process(Aaux, Acp);

% Generative model: Generate probabilistic mapping from hidden states to
% outcomes: a matrix
a{1} = repmat(acp, 1,size(A{1},2));
a{2} = A{2}*128;

% MDP structure
%--------------------------------------------------------------------------
mdp.A  = A;
mdp.As = As;
mdp.a  = a;
mdp    = spm_MDP_check(mdp);
ntime  = 16;    % number of moves per path
ntrial = 4;     % number of paths
mdpt   = mdp;

for i = 1:ntrial
    SDP  = spm_maze_active(mdpt,ntime,START,END);
    mdpt = SDP(end);
    trial(i).MDP = SDP;
end

% show results - behavioural
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_print_agentstrials(trial, START, END, a)


%% Simulation 2: sensitivity to concentration parameters
%==========================================================================

% Show typical trial for
%--------------------------------------------------------------------------
acp    = [1/8,1/2,2;1/8,1/2,2];   % Initial concentration parameters of agent (cp1 = cp2)
ntime  = 16;                      % number of timesteps per trial
na     = length(acp);
ntrial = 4;
trial  = [];

for i = 1:na
    mdp.a{1} = repmat(acp(:,i), 1,size(A{1},2));
    mdpt = mdp;
    for j = 1:ntrial
        SDP = spm_maze_active(mdpt,ntime,START,END);
        mdpt = SDP(end);
    end
    trial(i).MDP = SDP;
end
spm_figure('GetWin','Figure 2'); clf
spm_plot_3a(trial, START, END, acp);

%% Simulation 3: sensitivity of concentration parameters agent and environment
%======================================================================
acp   = [1/8,4; 1/8,4];      % Initial concentration parameters of agent (cp1 = cp2)

Acp(1).open   = [1;0];
Acp(1).closed = [0;1];
Acp(2).open   = [16;0];
Acp(2).closed = [0;16];

ntime  = 16;                 % number of moves per path
ntrial = 32;                 % number of paths

% Inital concentration parameters
%--------------------------------------------------------------------------
na    = size(acp,2);
nenv  = size(Acp,2);
agent = [];

for i = 1:na
    mdp.a{1} = repmat(acp(:,i), 1,size(A{1},2));
    for j = 1:nenv
        [A,As] = spm_gen_process(Aaux, Acp(j));
        mdp.A  = A;
        mdp.As = As;
        mdpt   = mdp;
        for k = 1:ntrial
            SDP  = spm_maze_active(mdpt,ntime,START,END);
            mdpt = SDP(end);
            agent(i).env(j).trial(k).MDP = SDP;
        end
    end
end
spm_figure('GetWin','Figure 3'); clf
spm_plot_3b(agent, START, END);

%% Simulation 4: Singular value decomposition
%==========================================================================
spm_figure('GetWin','Figure 4'); clf
nenv  = 2;
ntt   = ntime*ntrial;
c     = linspace(10,ntt/4,ntt);
for i = 1:na
    for j = 1:nenv
        Pag  = [];
        Penv = [];
        for k = 1:ntrial
            SDP = agent(i).env(j).trial(k).MDP;
            for l = 1:ntime
                
                % get expectations
                %----------------------------------------------------------
                Pa   = SDP(l).a{1};
                Pa   = Pa/diag(sum(Pa));
                Pag  = [Pag, Pa(1,:)'];
                Pa   = SDP(l).As;
                Pa   = Pa/diag(sum(Pa));
                Penv = [Penv, Pa(1,:)'];
            end
        end
        Pag  = Pag  -  Pag(:,end)*ones(1,ntt);
        Penv = Penv - Penv(:,end)*ones(1,ntt);
        X    = [Pag'; Penv'];
        
        % eigenvariates
        %------------------------------------------------------------------
        subplot(2,2,2*(i-1) + j)
        [U,S] = svd(X);
        U   = U*S;
        V1a = U((1:ntt),1);
        V2a = U((1:ntt),2);
        V1e = U((1:ntt) + ntt,1);
        V2e = U((1:ntt) + ntt,2);
        
        scatter(V1a, V2a, [], c, 'filled'), hold on
        scatter(V1e, V2e, [], c, 'ro');
        axis([-1 1 -1 1]*2.5), axis square, hold off
        xlabel('First principle component')
        ylabel('Second principle component')
        title('phenotypic trajectories','FontSize',16)
    end
end

% get and show various free energies
%--------------------------------------------------------------------------
for i = 1:na
    for j = 1:nenv
        for k = 1:ntrial
            SDP = agent(i).env(j).trial(k).MDP;
            for l = 1:ntime;
                Gu     = log(spm_softmax(SDP(l).G));
                r(l,k) = mean(SDP(l).rt);
                f(l,k) = mean(SDP(l).H);
                g(l,k) = mean(SDP(l).P'*Gu);
            end
        end
        R(i,j,:) = sum(r);
        F(i,j,:) = sum(f);
        G(i,j,:) = sum(g);
    end
end

% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
col = {'k','r';':k',':r'};

for i = 1:na
    for j = 1:nenv
        subplot(2,2,1)
        plot(spm_vec(F(i,j,:)),col{i,j}), hold on
        xlabel('exposure'), ylabel('seconds')
        title('(negative) free energy','FontSize',16)
        axis square, spm_axis tight, %set(gca,'YLim',[-3 1])
        
        subplot(2,2,2)
        plot(spm_vec(G(i,j,:)),col{i,j}), hold on
        xlabel('exposure'), ylabel('seconds')
        title('expected free energy','FontSize',16)
        axis square, spm_axis tight, %set(gca,'YLim',[-16 2])
        
        subplot(2,1,2)
        plot(spm_vec(R(i,j,:)),col{i,j}), hold on
        xlabel('exposure'), ylabel('seconds')
        title('reaction time','FontSize',16)
        axis square, spm_axis tight
  
    end
end

legend({'A:low-E:low','A:low-E:high','A:high-E:low','A:high-E:high'})


%%
return



function MDP = spm_maze_active(mdp,N,START,END)
% FORMAT MDP = spm_maze_search(mdp,N,START,END,alpha,beta,z)
% mdp   - MDP structure
% N     - number of moves (i.e., default 8)
% START - index of intial state (default 1)
% END   - index of target state (default 1)
% MDP   - MDP structure array

% preliminaries
%--------------------------------------------------------------------------
try, N;     catch, N     = 8;   end
try, START; catch, START = 1;   end
try, END;   catch, END   = 1;   end


% initialise concentration parameters: a (if unspecified)
%--------------------------------------------------------------------------
if ~isfield(mdp,'o')
    mdp.o = [];
end
if ~isfield(mdp,'u')
    mdp.u = [];
end
mdp.s    = START;

% Evaluate a sequence moves
%==========================================================================
for i = 1:N
    
    % Evaluate preferred states (subgoals) on the basis of current beliefs
    %----------------------------------------------------------------------
    mdp.C{2} = spm_maze_cost(mdp,END);
    
    % proceed with subsequent move
    %----------------------------------------------------------------------
    MDP(i)   = spm_MDP_VB_X(mdp);
    mdp      = MDP(i);
    mdp.s    = mdp.s(:,end);
    mdp.D{1} = MDP(i).X{1}(:,end);
    mdp.o    = [];
    mdp.u    = [];
    mdp.a{1} = MDP(i).a{1};
    
    % UPDATE: Change A matrix dependent on current position of agent
    %----------------------------------------------------------------------
    pos      = MDP(i).s(:,end);
    
    % Update environment
    %----------------------------------------------------------------------
    MDP(i).As(1,pos) = MDP(i).As(1,pos) + 1;
    
    % place in MDP structure
    %----------------------------------------------------------------------
    mdp.As   = MDP(i).As;
    mdp.A{1} = MDP(i).As/diag(sum(MDP(i).As));
    
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
X     = X > exp(-3);                    % Default cut-off is -3
C     = log(X.*(P*Y) + exp(-32));



return


function spm_plot_3a (trial, START, GOAL, cp_ag)
% unpack properties
ntrial = size(trial,2);
ntime  = size(trial(1).MDP,2);
A  = spm_vec(trial(1).MDP(1).A{1}(1,:));
ns = numel(A);
ni = sqrt(ns);

for j = 1:ntrial
    
    A  = spm_vec(trial(j).MDP(end).A{1}(1,:));
    A  = reshape(A,ni,ni);
    subplot(2,ntrial,j), imagesc(A, [0,1]), axis image
    title(['CP = ' num2str(cp_ag(1,j))], 'fontsize',14)
    hold on, colormap gray
    
    % Plot trajectory
    %----------------------------------------------------------------------
    for k = 1:ntime
        step    = trial(j).MDP(k).s;
        [l1,m1] = ind2sub([ni,ni],step(1));
        [l2,m2] = ind2sub([ni,ni],step(2));
        plot([m1,m2], [l1,l2],':m');
    end
    
    
    [p,q] = ind2sub([ni,ni],START);
    plot(q,p,'.','MarkerSize',32,'Color','g');
    [p,q] = ind2sub([ni,ni],GOAL);
    plot(q,p,'.','MarkerSize',32,'Color','r');
    [p,q] = ind2sub([ni,ni],trial(j).MDP(end).s(2));
    plot(q,p,'.','MarkerSize',32,'Color','b');
    
    hold off
    
    % Plot Likelihood matrix
    %----------------------------------------------------------------------
    Q   = trial(j).MDP(end).a{1};
    Q   = Q/diag(sum(Q));
    Q   = Q(1,:);
    a   = reshape(Q(:),ni,ni);
    subplot(2,ntrial,j+3), imagesc(a, [0,1]), axis image
    title('Likelihood','fontsize',14);
    
end

return

function spm_plot_3b (agent, START, GOAL)

% unpack properties
%--------------------------------------------------------------------------
nag   = size(agent,2);
nenv  = size(agent(1).env,2);
ntime = size(agent(1).env(1).trial(1).MDP,2);
A     = spm_vec(agent(1).env(1).trial(1).MDP(1).A{1}(1,:));
ns    = numel(A);
ni    = sqrt(ns);

for i = 1:nenv
    for j = 1:nag
        A  = spm_vec(agent(j).env(i).trial(end).MDP(end).A{1}(1,:));
        A  = reshape(A,ni,ni);
        subplot(nenv,nag*2,(2*(j-1))+1+(2*nenv*(i-1)))
        imagesc(A, [0,1]), axis image
        
        hold on
        colormap gray
        % Plot trajectory
        %------------------------------------------------------------------
        for k = 1:ntime
            step    = agent(j).env(i).trial(end).MDP(k).s;
            [l1,m1] = ind2sub([ni,ni],step(1));
            [l2,m2] = ind2sub([ni,ni],step(2));
            plot([m1,m2], [l1,l2],':r');
        end
        
        [p,q] = ind2sub([ni,ni],START);
        plot(q,p,'.','MarkerSize',32,'Color','g');
        [p,q] = ind2sub([ni,ni],GOAL);
        plot(q,p,'.','MarkerSize',32,'Color','r');
        [p,q] = ind2sub([ni,ni],agent(j).env(i).trial(end).MDP(end).s(2));
        plot(q,p,'.','MarkerSize',32,'Color','b');
        
        
        % Plot Likelihood matrix
        %------------------------------------------------------------------
        Q   = agent(j).env(i).trial(end).MDP(end).a{1};
        Q   = Q/diag(sum(Q));
        Q   = Q(1,:);
        a   = reshape(Q(:),ni,ni);
        subplot(nenv,nag*2,(2*(j-1))+2+(2*nenv*(i-1))), imagesc(a, [0,1]), axis image
        
        hold off
        
    end
end
return

function [A, As] = spm_gen_process(Aaux, Acp)
% This function generates the observation matrix A of the generative
% process based on the layout of the maze Aaux and the environmental
% concentration parameters Acp
%--------------------------------------------------------------------------
initwh = Acp.open;
initbl = Acp.closed;
As     = initwh*Aaux{1}(1,:)+initbl*Aaux{1}(2,:);
Ap     = repmat(sum(As,1),2,1);

A{1} = As./Ap;    % Expected value of dirichlet is alpha1/(alpha1+alpha2);
A{2} = Aaux{2};

return

function spm_print_agentstrials(trial, START, GOAL,ainit)
%--------------------------------------------------------------------------
ntrial = size(trial,2);
ntime  = size(trial(1).MDP,2);
A  = spm_vec(trial(1).MDP(1).A{1}(1,:));
a  = ainit;
ns = numel(A);
ni = sqrt(ns);

% plot initial maze
%--------------------------------------------------------------------------
A  = spm_vec(trial(1).MDP(1).A{1}(1,:));
A  = reshape(A,ni,ni);
subplot(ntrial+1,2,1), imagesc(A, [0,1]), axis image
title('Initial maze', 'fontsize',14)
hold on
colormap gray
try
    [p,q] = ind2sub([ni,ni],START);
    plot(q,p,'.','MarkerSize',32,'Color','g');
    [p,q] = ind2sub([ni,ni],GOAL);
    plot(q,p,'.','MarkerSize',32,'Color','r');
end
hold off

% plot initial likelihood
%--------------------------------------------------------------------------
Q = a{1};
Q = Q/diag(sum(Q));
Q = Q(1,:);
a = reshape(Q(:),ni,ni);
subplot(ntrial+1,2,2), imagesc(a, [0,1]), axis image
title('Initial Likelihood','fontsize',14);

for j = 1:ntrial
    A = spm_vec(trial(j).MDP(end).A{1}(1,:));
    A = reshape(A,ni,ni);
    subplot(ntrial+1,2,1+(2*j)), imagesc(A, [0,1]), axis image
    title(['Trial ' num2str(j)], 'fontsize',14)
    hold on
    % Plot trajectory
    %----------------------------------------------------------------------
    for k = 1:ntime
        step = trial(j).MDP(k).s;
        [l1,m1]= ind2sub([ni,ni],step(1));
        [l2,m2]= ind2sub([ni,ni],step(2));
        p = plot([m1,m2], [l1,l2],':m');
    end
    
    % try
    %----------------------------------------------------------------------
    [p,q] = ind2sub([ni,ni],START);
    plot(q,p,'.','MarkerSize',32,'Color','g');
    [p,q] = ind2sub([ni,ni],GOAL);
    plot(q,p,'.','MarkerSize',32,'Color','r');
    [p,q] = ind2sub([ni,ni],trial(j).MDP(end).s(2));
    plot(q,p,'.','MarkerSize',32,'Color','b');
    
    
    hold off
    
    % Plot Likelihood matrix
    %----------------------------------------------------------------------
    Q     = trial(j).MDP(end).a{1};
    Q     = Q/diag(sum(Q));
    Q     = Q(1,:);
    a     = reshape(Q(:),ni,ni);
    subplot(ntrial+1,2,2+(2*j)), imagesc(a, [0,1]), axis image
    title(['Likelihood after trial ' num2str(j)],'fontsize',14);
end

return
