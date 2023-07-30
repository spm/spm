function MDP = DEM_Tower_of_Hanoi
% Demo of planning as sophisticated inference (Tower of Hanoi)
%__________________________________________________________________________
%
% This demonstration routine illustrates planning as sophisticated
% inference in the context of the Tower of Hanoi (or block-world) paradigm.
% This is effectively a game in which coloured blocks have to be moved
% between (usually three) pillars or towers to achieve a target
% arrangement. The generative model in this instance is based upon the
% notion of latent states as arrangements of multiple objects; where the
% physics — underwriting the generative model — prescribes allowable
% transitions among arrangements. In this example, a generative process is
% set up by associating every possible arrangement of blocks in three
% towers with latent states. These latent states generate outcomes for each
% location (i.e., a coloured block or nothing for each visual modality).
% Transitions among arrangements are modelled as seven controllable paths,
% corresponding to choosing a particular pillar and moving the top block to
% another pillar (or doing nothing).
%
% The first demonstration simply illustrates the planning as inference. We
% then demonstrate how this generative model is learned by presenting
% exemplar arrangements (and allowable transitions) to a naïve agent – that
% can then plan, given some prior preferences about a target
% distribution — and how performance depends upon the depth of planning.
%_________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%==========================================================================
rng(1)

% size of problem
%--------------------------------------------------------------------------
Nb = 3;                                          % number of blocks
Nc = 3;                                          % number of columns

%% find all arrangements of blocks
%--------------------------------------------------------------------------
% Effectively, place an ordered sequence of blocks over all combinations of
% pillars – and do this for all sequences (i.e., permutations) of coloured
% blocks. The resulting arrangement is encoded in row vectors (A), whose
% elements correspond to the colour of blocks (with zero being an empty
% location). The location space is just the number of pillars times the
% number of blocks (i.e., the maximum height of a tower).
%--------------------------------------------------------------------------
Pb    = perms(1:Nb);                             % permutations of blocks
Cc    = spm_combinations(ones(Nb,1)*Nc);         % combinations of columns
A     = [];
for b = 1:size(Pb,1)
    for c = 1:size(Cc,1)
        a     = zeros(Nb,Nc);
        for i = 1:Nb
            j = Cc(c,i);
            a(find(~a(:,j),1),j) = Pb(b,i);
        end
        A(end + 1,:)   = a(:)';
    end
end

% find unique arrangements
%--------------------------------------------------------------------------
A   = unique(A,'rows');
Ns  = size(A,1);

%% allowable transitions
%--------------------------------------------------------------------------
% Effectively, choose a pillar and move the top block to each pillar in
% turn; starting from a particular arrangement. This movement basically
% changes the top block of the chosen pillar to 0 (an empty location) and
% populates the target pillar with the moved block. The ensuing arrangement
% is then associated with the corresponding state to build a state
% transition matrix, with allowable transitions over the nine possible
% moves (including two of the three redundant transitions; namely, no
% move). This means that there are seven potential actions to be taken from
% every arrangement.
%--------------------------------------------------------------------------
B     = zeros(Ns,Ns,0);                       % transition matrices
for c = 1:Nc                                  % source column
    for t = 1:Nc                              % target column
        b     = zeros(Ns,Ns);                 % transitions for this move
        for s = 1:Ns                          % for each arrangement

            % arrangement
            %--------------------------------------------------------------
            a = reshape(A(s,:)',Nb,Nc);
            h = find(a(:,c),1,'last');        % height of source column
            g = find(a(:,t),1,'last');        % height of source column

            if isempty(h), h = 0; end
            if isempty(g), g = 0; end
            if h && c ~= t

                % move block
                %----------------------------------------------------------
                d          = a(h,c);
                a(g + 1,t) = d;
                a(h,c)     = 0;

                % new arrangement
                %----------------------------------------------------------
                i      = ismember(A,a(:)','rows');
                b(i,s) = 1;

            else
                b(s,s) = 1;
            end

        end

        % place in transition tensor
        %------------------------------------------------------------------
        B(:,:,end + 1) = b;
        
    end
end


% create likelihood tensors
%--------------------------------------------------------------------------
G     = {};
for s = 1:Ns
    a     = reshape(A(s,:)',Nb,Nc);
    for c = 1:Nc
        for h = 1:Nb
            G{h,c}(:,s)   = full(sparse(a(h,c) + 1,1,1,Nb + 1,1));
        end
    end
end
A     = G(:);
B     = {B};

% Enumerate the states and paths of the ensuing generative model
%--------------------------------------------------------------------------
Nf    = numel(B);                    % number of hidden factors
Ng    = numel(A);                    % number of outcome modalities
for f = 1:Nf
    Ns(f) = size(B{f},1);
    Nu(f) = size(B{f},3);
end
for g = 1:Ng
    No(g) = size(A{g},1);
end


%% priors: (cost) C
%--------------------------------------------------------------------------
% Finally, specify the prior preferences in terms of log probabilities over
% outcomes:

% Preferences for target
%--------------------------------------------------------------------------
% These are simply the outcomes for a particular target arrangement. The
% target itself specified in terms of the starting arrangement and applying
% n (nontrivial) transitions to create easy or difficult problems.
%--------------------------------------------------------------------------
s     = 1;
n     = 8;
r     = spm_target(s,n,B{1});         % Nb = 3, r = 5; for a difficult problem
C     = spm_cost(r,A);

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify (uninformative) prior
% beliefs about initial states (D) and paths through those states (E)
%--------------------------------------------------------------------------
for f = 1:numel(B)
    D{f} = ones(Ns(f),1);
    E{f} = ones(Nu(f),1);
    H{f} = ones(Ns(f),1);
end
H{f}(r)  = 32;                         % 32 for inductive inference 


% specify controllable factors; here, the first and only factor
%--------------------------------------------------------------------------
U     = 1;                            % controllable factors

% MDP Structure, specifying 64 epochs (i.e., 16 seconds of active vision)
%==========================================================================
MDP.T = 8;                            % numer of moves
MDP.U = U;                            % controllable factors
MDP.A = A;                            % likelihood probabilities
MDP.B = B;                            % transition probabilities
MDP.C = C;                            % prior preferences
MDP.D = D;                            % prior over initial states
MDP.H = H;                            % prior over final states
MDP.E = E;                            % prior over initial paths
MDP.N = 3;                            % planning depth (-1)

% Solve an example with known (veridical) structure
%==========================================================================

% specify a trial, with initial conditions (s)
%--------------------------------------------------------------------------
MDP.s = s;                            % initial state
PDP   = spm_MDP_VB_XXX(MDP);

% illustrate behavioural responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(PDP);

% illustrate action and inference about states and paths
%--------------------------------------------------------------------------
spm_figure('GetWin','optimal behaviour'); clf
spm_report(PDP,Nb,Nc,r)


% Illustrate structure learning via assimilation of epochs of observations
%==========================================================================
OPTIONS.N = 0;       % suppress neuronal responses
OPTIONS.P = 0;       % suppress plotting
OPTIONS.B = 1;       % replay
OPTIONS.G = 1;       % suppress graphics

OPTIONS.NF = 1;      % maxmium number of factors
OPTIONS.NS = Ns;     % maxmium number of states
OPTIONS.NU = Nu;     % maxmium number of paths

% Generate training sequence of paths(length T)
%--------------------------------------------------------------------------
O     = spm_MDP_structure_teaching(MDP);

% Structure learning
%--------------------------------------------------------------------------
mdp.p = 1/512;      % adjust to make the agent more or less impressionable
mdp   = spm_MDP_structure_learning(mdp,O,OPTIONS);

% Illustrate allowable and discovered transitions (i.e., the physics)
%--------------------------------------------------------------------------
spm_figure('GetWin','generative process'); clf
spm_MDP_params(MDP)
spm_figure('GetWin','before motor babbling'); clf

spm_MDP_params(mdp)

% active learning (c.f., motor babbling)
%==========================================================================
% So far, structure learning is incomplete because there are certain
% transitions that have not been experienced. We can now use the expected
% information gain about model parameters to drive motor explanation; under
% uninformative prior preferences.
%
% Note that the specification of the generative process and the structure
% learning are not necessarily isomorphic. This means explicit action needs
% to be modelled; namely, the action associated with the generative process
% is chosen to minimise the free energy of the outcome. This is just the
% predictive accuracy, based inference under the generative model.
%--------------------------------------------------------------------------
OPTIONS.A = 1;                           % enable explicit action
OPTIONS.B = 0;                           % suppress replay
OPTIONS.N = 0;                           % suppress neuronal responses
OPTIONS.G = 0;                           % suppress graphics

mdp.A = MDP.A;                           % Generative process (likelihood)
mdp.B = MDP.B;                           % Generative process (transition)
mdp.U = 1;
mdp.N = 1;

% engage active learning (i.e., accept parameter updates if G increases)
%--------------------------------------------------------------------------
mdp.beta = 512;

% Illustrate exploration of states and paths
%--------------------------------------------------------------------------
spm_figure('GetWin','motor learning'); clf
mdp.T = 64;
pdp   = spm_MDP_VB_XXX(mdp,OPTIONS);             % invert
spm_MDP_VB_trial(pdp);                           % report

% Illustrate the accompanying increase in G or mutual information
%--------------------------------------------------------------------------
spm_figure('GetWin','motor babbling'); clf
NT    = 256;
mdp.T = 8;
MI    = [];
for t = 1:NT

    % active inference and learning
    %----------------------------------------------------------------------
    pdp  = spm_MDP_VB_XXX(mdp,OPTIONS);          % invert
    mdp  = spm_MDP_VB_update(mdp,pdp,OPTIONS);   % update parameters


    % report graphically
    %======================================================================
    if t == 1
        spm_report(pdp,Nb,Nc);
    end

    % expected free energy (mutual information)
    %----------------------------------------------------------------------
    MI(t,1) = spm_MDP_MI(mdp.a);
    MI(t,2) = spm_MDP_MI(mdp.b);

    subplot(2,2,3)
    plot(1:t,MI(:,1)), set(gca,'XLim',[1,NT])
    title({'Mutual information (a) ','(likelihood)'})
    xlabel('time'), ylabel('nats'), box off, axis square

    subplot(2,2,4)
    plot(1:t,MI(:,2)), set(gca,'XLim',[1,NT])
    title({'Mutual information (b) ','(priors)'})
    xlabel('time'), ylabel('nats'), box off, axis square
    drawnow

end

% consolidate learning, using Bayesian model reduction    
%--------------------------------------------------------------------------
for g = 1:Ng
     mdp.a{g} = spm_MDP_VB_prune(mdp.a{g},mdp.p,0,0,0,'SIMPLE')*512;
end
for f = 1:Nf
     mdp.b{f} = spm_MDP_VB_prune(mdp.b{f},mdp.p,0,0,0,'SIMPLE')*512;
end

% Illustrate learned transition priors
%--------------------------------------------------------------------------
spm_figure('GetWin','after motor babbling'); clf

spm_MDP_params(pdp)

subplot(3,2,4)
imagesc(~~sum(MDP.B{1} > 1/2,3)),axis square
title('Allowable transitions','FontSize',14), axis square

subplot(3,2,6)
imagesc(~~sum(pdp.b{1} > 1/2,3)),axis square
title('Discovered transitions','FontSize',14), axis square
xlabel('latent states'), drawnow


% Illustrate solutions to increasingly difficult problems, with increasing
% depth of planning
%==========================================================================
mdp.T = 8;                                % maximum number of moves
NT    = 100;                              % number of trials
N     = [0 1 2 3 4];                      % planning depths
C     = zeros(NT,numel(N));               % number of inferred moves
for t = 1:NT

    % target
    %----------------------------------------------------------------------
    s     = randperm(Ns,1);               % random source configuration
    n     = randperm(5,1);                % random number of moves
    r     = spm_target(s,n,B{1});         % target configuration
    mdp.C = spm_cost(r,A);
    mdp.s = s;

    for n = 1:numel(N)

        % depth of planning
        %------------------------------------------------------------------
        mdp.N = N(n);

        % active inference and learning
        %------------------------------------------------------------------
        pdp  = spm_MDP_VB_XXX(mdp,OPTIONS);

        % spm_report(pdp,Nb,Nc); record inferred moves
        %------------------------------------------------------------------
        try
            C(t,n) = find(pdp.s == r,1);
        end
        if C(t,n)
            C(t,n:end) = C(t,n);     
            break
        end
    end
end

% report graphically
%==========================================================================
spm_figure('GetWin','performance'); clf

% Maximum number of moves
%--------------------------------------------------------------------------
c     = C';
c(~c) = 32;
NM    = min(c);

for n = 1:numel(N)
    c = C(:,n);
    for m = 1:max(C,[],'all')
        i      = find(NM == m);
        Nm(m)  = numel(i);
        P(n,m) = 100*sum(c(i) > 1)/Nm(m);
        str{m} = sprintf('%i moves',m - 1);
    end
end

% plot Performance and Number of trials
%--------------------------------------------------------------------------
subplot(2,2,1)
bar(P(:,2:end))
title('Performance','FontSize',14)
xlabel('depth of planning'), ylabel('success (%)'), box off, axis square
legend(str(2:end))

subplot(2,4,3)
bar(Nm(2:end))
title(sprintf('Number of trials (%i/%i)',sum(Nm),NT),'FontSize',14)
xlabel('number of moves'), ylabel('incidence'), box off, axis square


return


% subroutines
%==========================================================================

function spm_report(MDP,Nb,Nc,r)
% Plots the inferred sequence of moves
%--------------------------------------------------------------------------
% If there are more than eight moves this subroutine will plot a movie

%% preliminaries
%==========================================================================
col = {'w','r','g','b','y','m','c','k'};
n   = 4;

% Plot target
%--------------------------------------------------------------------------
subplot(n,3,n*3 - 1), hold off
for c = 1:Nc,line([c,c],[0,(Nb + 1)],'LineWidth',4),end
axis([0 (Nc + 1) 0 (Nb + 1)]), hold on

% balls
%--------------------------------------------------------------------------
if nargin < 4
    C = MDP.C;
else
    C = spm_cost(r,MDP.A);
end
C     = reshape(C,Nb,Nc);
for h = 1:Nb
    for c = 1:Nc
        [d,o] = max(C{h,c});
        if o > 1
            plot(c,h,'.','MarkerSize',64,'Color',col{o}), hold on
        end
    end
end
title('Target','FontSize',12)
axis square, drawnow


%% show sequence of moves
%==========================================================================
for t = 1:MDP.T

    % movie
    %----------------------------------------------------------------------
    if MDP.T > 8
        subplot(4,3,2), cla
    else
        subplot(n,n,t), cla
    end

    % plot
    %----------------------------------------------------------------------
    for c = 1:Nc,line([c,c],[0,(Nb + 1)],'LineWidth',4),end
    axis([0 (Nc + 1) 0 (Nb + 1)]), hold on
    o     = reshape(MDP.o(:,t),Nb,Nc);
    for h = 1:Nb
        for c = 1:Nc
            if o(h,c) > 1
                plot(c,h,'.','MarkerSize',64,'Color',col{o(h,c)}), hold on
            end
        end
    end
    title(sprintf('trial %i',t),'FontSize',12)
    axis square, hold off, drawnow

    % movie
    %----------------------------------------------------------------------
    if MDP.T > 8
        I(t) = getframe(gca); pause(1/4)
    end
end

% Place movie in graphic subject
%--------------------------------------------------------------------------
if MDP.T > 8
    subplot(4,3,2)
    set(gca,'Userdata',{I,32})
    set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
end




return

function t = spm_target(s,n,B)

% Returns a target state for any given starting state
%--------------------------------------------------------------------------
% s - initial state
% n - number of moves
% B - transition tensor
%
% Note that there may be shorter paths to the target than computed with
% this subroutine.
%__________________________________________________________________________

% Compose a sequence of non-trivial transitions
%==========================================================================
t     = s;
for m = 1:n
    [i,j] = find(squeeze(B(:,t,:)));  % find next state
    i     = i(~ismember(i,s));        % prohibit return to a previous state
    if numel(i)
        t = i(randperm(numel(i),1));  % choose a transition at random
        s = [s,t];
    else
        break
    end
end

return

function C = spm_cost(r,A)

% Returns a cost for any given target state
% r - target state
% A - likelihood tensor
% B - transition tensor
%--------------------------------------------------------------------------
% The cost in outcome space is just the likelihood mapping from the
% (target) latent state
%__________________________________________________________________________

% Cost
%==========================================================================
for g = 1:numel(A)
    C{g} = A{g}(:,r);
end

return

%% routines that call spm_MDP_VB_XXX
%--------------------------------------------------------------------------
% Routine:                   Demonstrating
%--------------------------------------------------------------------------
DEM_demo_MDP_XXX            % context learning
DEMO_MDP_maze_X             % no inductive inference
DEMO_MDP_maze_XXX           % Expected information gain (parameters)
DEM_surveillance            % Factorial problem
DEM_sharingX                % Active selection with hyperlinks
DEM_dSprites                % Structure learning with dynamics
DEM_Tower_of_Hanoi          % Active inference and structure learning
DEM_MNIST                   % Active selection without dynamics

%%% DEM_syntax




