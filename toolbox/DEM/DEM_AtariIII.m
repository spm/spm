function MDP = DEM_AtariIII
% de novo Structure learning — of pullback attractors —  afrom pixels 
%__________________________________________________________________________
%
% This routine illustrates the de novo learning of a state action policy —
% and goal-directed behaviour using inductive inference under the ensuing
% model. The approach is based upon the fast structure learning of dynamics
% (i.e., transitions among generalised or deep states) that possess a
% nonequilibrium steady-state distribution (NESS). This is guaranteed when
% all states are transient states and there are no absorbing states. In
% this example, random play is generated for 10,000 frames of a simplified
% Atari game (Pong). The first thousand frames are then used to learn the
% structure of a renormalising generative model; namely, the recursive
% composition of states and paths specified by edges or dependencies (i.e.,
% indices of parents). The resulting structure is then used to accumulate
% sequences of play that, by chance, include rewarded outcomes (precluding
% costly or constrained outcomes). Following this active selection of
% training data, the training data are used to accumulate transitions to,
% and only to, the basins of attraction of goal states. When all goal
% states lie in a basin of attraction — of other goal states — the
% accumulation terminates. At this point, the dynamics possess a
% nonequilibrium steady-state, which is used to reduce the model by
% successively removing states with a low NESS density, until the emergence
% of an absorbing state.
%
% The ensuing model is then demonstrated in generative mode by generating
% (perfect) fantasy play. The model is then used for real gameplay by
% installing the generative process in the model and using explicit action
% to learn from any mistakes: inductive inference is used to find the path
% of least action from one goal state to the next and new paths are
% retained, provided they do not incur any costs. Finally, the model is
% demonstrated after a final reduction or compression by retaining deep
% (generalised) states that predict current and future actions.
%
% This form of de novo structure learning uses a generative model that
% generates predictions in multiple streams. The principal (leading) stream
% generates pixel data, while subsequent streams generate rewards, costs
% and the consequences of action (e.g., proprioception or telemetry). By
% having separate streams, one can find an efficient or reduced model that
% preserves the mutual information between generalised states at the
% highest level and future action. These deep states are effectively
% multimodal and serve to couple pixel level states to action; thereby
% affording an efficient state action policy — that is realised through
% inductive inference.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

% set up and preliminaries
%--------------------------------------------------------------------------
rng(2)

% Get game: i.e., generative process (as a partially observed MDP)
%==========================================================================
Nr = 12;                                     % number of rows
Nc = 9;                                      % number of columns
Sc = 3;                                      % scaling
Nd = 4;                                      % random initial conditions
C  = 32;                                     % log cost

% get game in MDP form
%--------------------------------------------------------------------------
[GDP,~,~,~,RGB] = spm_MDP_pong(Nr,Nc,Nd,true,0);

% size of streams
%--------------------------------------------------------------------------
S      = ones(4,3);
S(1,:) = [Nr,Nc,1];                          % sensory stream
S(2,:) = [1 1 1];                            % reward  stream                 
S(3,:) = [1 1 1];                            % cost    stream
S(4,:) = [1 1 1];                            % policy  stream

% reward and cost functions [of outcomes]
%--------------------------------------------------------------------------
spm_get_hits = @(o,id) find(o(id.reward,:)    > 1);
spm_get_miss = @(o,id) find(o(id.contraint,:) > 1);

% Generate (probabilistic) outcomes under random actions
%==========================================================================
spm_figure('GetWin','Gameplay'); clf

GDP.tau = 2;                                 % smoothness of random paths
GDP.T   = 10000;                             % training length
PDP     = spm_MDP_generate(GDP);             % generate play

% illustrate sequence of random play
%---------------------------------------------------------------------------
con   = PDP.id.control;
for t = 1:128
    subplot(2,1,1)
    imshow(spm_O2rgb(PDP.O(:,t),RGB))
    subplot(4,3,8)
    imshow(PDP.O{con,t}')
    drawnow
end

% initial structure learning: grouping operators (iA,iB,iC,...)
%==========================================================================
MDP = spm_faster_structure_learning(PDP.O(:,1:1000),S,Sc);

% forget parameters to select rewarded episodes
%--------------------------------------------------------------------------
Nm    = numel(MDP);
Ne    = max(2^(Nm - 1),1);
for n = 1:Nm
    for g = 1:numel(MDP{n}.a)
        MDP{n}.a{g} = [];
    end
    for f = 1:numel(MDP{n}.b)
        MDP{n}.b{f} = [];
    end
end

% find rewarded and costly events
%--------------------------------------------------------------------------
r     = spm_get_hits(PDP.o,GDP.id);
c     = spm_get_miss(PDP.o,GDP.id);
for i = 1:numel(r)

    % for each sequence ending with an intended outcome
    %----------------------------------------------------------------------
    s  = c(find(c < r(i),1,'last'));
    t  = (s + Ne):(r(i) + Ne);
    if numel(t)
       
        % assimilate this sequence
        %------------------------------------------------------------------
        for s = 1:Ne
            MDP = spm_merge_structure_learning(PDP.O(:,t + s),MDP);
        end
    end
end

% Step through training data to enlarge basins of attraction to goal states
%==========================================================================
spm_figure('GetWin','Attractors'); clf

NT = 100;                                    % number of outcomes
NS = [];                                     % number of states  
NU = [];                                     % number of paths 
NA = [];                                     % number of absorbing states
for i = 1:256

    % Accumulate these states under random play
    %----------------------------------------------------------------------
    q     = rem(i,100 - 1);
    t     = (0:(NT + Ne)) + q*NT;
    for s = 1:Ne
        MDP = spm_merge_structure_learning(PDP.O(:,t + s),MDP);
    end

    % test for NESS
    %----------------------------------------------------------------------
    [MDP,d] = spm_RDP_basin(MDP,[2,3],[C,-C]);

    NS(end + 1) = size(MDP{Nm}.b{1},2);    
    NU(end + 1) = size(MDP{Nm}.b{1},3);
    NA(end + 1) = sum(~d);

    subplot(2,3,1), plot(NS), title('Deep states'),      axis square
    subplot(2,3,2), plot(NU), title('Deep paths'),       axis square
    subplot(2,3,3), plot(NA), title('Absorbing states'), axis square
    drawnow

    % break if all (deep) states are transient (i.e., no absorbing states)
    %----------------------------------------------------------------------
    if all(d), break, end

end

% Retain (and sort) states with a high NESS probability
%--------------------------------------------------------------------------
for q = 1:4
    MDP = spm_RDP_sort(MDP);
end

% Illustrate transitions in deep (generalised) state space 
%-=========================================================================
MDP   = spm_set_goals(MDP,[2,3],[C,-C]);
hid   = MDP{Nm}.id.hid;

subplot(2,2,3)
spm_dir_orbits(MDP{Nm}.b{1},hid,128);

% paths to hits
%--------------------------------------------------------------------------
subplot(2,2,4)
B     = sum(MDP{Nm}.b{1},3) > 0;
Ns    = size(B,1);
Nt    = 32;
h     = sparse(1,hid,1,1,Ns);
P     = zeros(Nt,Ns);
for t = 1:Nt
    P(t,:) = h;
    h      = (h + h*B) > 0;
end
imagesc(P), hold on 
plot(hid,zeros(size(hid)) + 1/2,'or','MarkerSize',8), hold off
title('Paths to hits','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square

% Generate play from recursive generative model
%==========================================================================

% assemble RGM
%--------------------------------------------------------------------------
RDP   = spm_set_goals(MDP,[2,3],[C,-C]);   % set intended states (hid/cid)
RDP   = spm_set_costs(RDP,[2,3],[C,-C]);   % set contraints (C)
RDP   = spm_mdp2rdp(RDP);                  % get nested model
RDP.T = 64;
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Generative AI'); clf
spm_show_RGB(PDP,RGB);

% and hits
%--------------------------------------------------------------------------
h     = spm_get_hits(PDP.Q.o{1},GDP.id);
subplot(Nm + 3,2,2*(Nm + 1))
plot(h,zeros(size(h)),'.r','MarkerSize',16)
drawnow


%% Active inference
%==========================================================================
% In what follows, we engage inductive planning as inference, with explicit
% action; namely, the generative process is used to engage actions that
% reproduce predicted outcomes (here, bat position) at the lowest level.
%--------------------------------------------------------------------------

% install generative process in MDP structure
%--------------------------------------------------------------------------
MDP{1}.GA  = GDP.A;
MDP{1}.GB  = GDP.B;
MDP{1}.GU  = GDP.U;
MDP{1}.GD  = GDP.D;
MDP{1}.ID  = GDP.id;
MDP{1}.chi = 512;                               % sticky action/shaky hand

% enable learning of transition parameters to induce exploration
%--------------------------------------------------------------------------
NT    = 256;                                    % length of each game
NR    = 32;                                     % number of replications
NS    = 256;                                    % concentration parameter
F     = NaN(6,NR);
for i = 1:NR

    % set intended states and constraints
    %----------------------------------------------------------------------
    RDP   = spm_set_goals(MDP,[2,3],[C,-C]);
    RDP   = spm_set_costs(RDP,[2,3],[C,-C]);

    % assemble RGM and reduce prior precision
    %----------------------------------------------------------------------
    RDP   = spm_mdp2rdp(RDP,0,1/NS);

    % play
    %----------------------------------------------------------------------
    RDP.T = fix(NT/Ne);
    PDP   = spm_MDP_VB_XXX(RDP);
    h     = spm_get_hits(PDP.Q.o{1},GDP.id);

    if true

        % Illustrate intentional play
        %------------------------------------------------------------------
        spm_figure('GetWin','Active inference'); clf
        spm_show_RGB(PDP,RGB,4,0);

        % and hits
        %------------------------------------------------------------------
        subplot(Nm + 3,2,2*(Nm + 1))
        plot(h,zeros(size(h)) - 2,'.r','MarkerSize',16)
        drawnow

    end

    % record structure learning and behaviour
    %---------------------------------------------------------------------
    F(1,i) = size(PDP.B{1},2);                  % number of (deep) states
    F(2,i) = size(PDP.B{1},3);                  % number of (deep) paths
    F(3,i) = PDP.Q.F + sum(PDP.F);              % ELBO
    F(4,i) = numel(h);                          % number of hits
    F(5,i) = size(MDP{end}.b{1},2);             % number of (deep) states
    F(6,i) = size(MDP{end}.b{1},3);             % number of (deep) paths

    % learn from mistakes
    %======================================================================

    % accumulate self-generated sequence (i.e., Dirichlet counts)
    %----------------------------------------------------------------------
    O     = PDP.Q.O{1};
    t     = 0:(NT - Ne);
    for s = 1:Ne
        MDP = spm_merge_structure_learning(O(:,t + s),MDP);
    end
    MDP   = spm_RDP_basin(MDP,[2,3],[C,-C]);    % basins of attraction

    % plot results
    %----------------------------------------------------------------------
    spm_figure('GetWin','Structure learning');
    str = {...
        'Latent states',...
        'Latent paths',...
        'ELBO',...
        'Reward count'};
    for f = 1:numel(str)
        subplot(3,2,f)
        plot(F(f,:)), axis square
        title(str{f},'FontSize',14), xlabel('games')
    end
    t  = (1:i)*NT;
    plot(t,F(f,1:i)), axis square
    title(str{f},'FontSize',14), xlabel('frames')
    drawnow

end

% Illustrate performance : before compression
%-=========================================================================

% assemble RGM
%--------------------------------------------------------------------------
RDP   = spm_RDP_sort(MDP);
RDP   = spm_set_goals(RDP,[2,3],[C,-C]);
RDP   = spm_set_costs(RDP,[2,3],[C,-C]);
RDP   = spm_mdp2rdp(RDP,0,1/NS);

% play
%--------------------------------------------------------------------------
RDP.T = 128;
PDP   = spm_MDP_VB_XXX(RDP);

spm_figure('GetWin','Active inference (before compression)'); clf
spm_show_RGB(PDP,RGB,8,false);

% and hits
%--------------------------------------------------------------------------
subplot(Nm + 3,2,2*(Nm + 1))
h   = spm_get_hits(PDP.Q.o{1},GDP.id);
plot(h,zeros(size(h)) - 2,'.r','MarkerSize',16)
drawnow

% Illustrate attractor 
%--------------------------------------------------------------------------
spm_figure('GetWin','Orbits'); clf
subplot(2,2,1)
HID   = PDP.id.hid;
spm_dir_orbits(PDP.B{1},HID,64);

% and paths to hits
%--------------------------------------------------------------------------
subplot(2,2,3)
B     = sum(PDP.B{1},3) > 1/32;
Ns    = size(B,1);
h     = sparse(1,HID,1,1,Ns);
I     = [];
for t = 1:32
    I(t,:) = h;
    h      = (h + h*B) > 0;
end
imagesc(I), hold on 
plot(HID,zeros(size(HID)) + 1/2,'or','MarkerSize',8), hold off
title('Paths to hits (before)','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square


% Illustrate performance : after compression
%-=========================================================================

% assemble RGM
%--------------------------------------------------------------------------
RDP   = spm_RDP_sort(MDP);                 % eigenreduction
RDP   = spm_RDP_MI(RDP);                   % merge conserving MI 
RDP   = spm_set_goals(RDP,[2,3],[C,-C]);   % set intended states (hid/cid)
RDP   = spm_set_costs(RDP,[2,3],[C,-C]);   % set contraints (C)
RDP   = spm_mdp2rdp(RDP,0,1/NS);           % get nested model

% active inference and learning
%--------------------------------------------------------------------------
RDP.T = 128;
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference (with compression)'); clf
spm_show_RGB(PDP,RGB,4,0);

% and hits
%--------------------------------------------------------------------------
h = spm_get_hits(PDP.Q.o{1},GDP.id);
subplot(Nm + 3,2,2*(Nm + 1))
plot(h,zeros(size(h)) - 2,'.r','MarkerSize',16)
drawnow

% Illustrate attractor 
%--------------------------------------------------------------------------
spm_figure('GetWin','Orbits');
subplot(2,2,2)
HID   = PDP.id.hid;
spm_dir_orbits(PDP.B{1},HID,64);

% paths to hits
%--------------------------------------------------------------------------
subplot(2,2,4)
B     = sum(PDP.B{1},3) > 1/32;
Ns    = size(B,1);
h     = sparse(1,HID,1,1,Ns);
I     = [];
for t = 1:32
    I(t,:) = h;
    h      = (h + h*B) > 0;
end
imagesc(I), hold on 
plot(HID,zeros(size(HID)) + 1/2,'or','MarkerSize',8), hold off
title('Paths to hits (after)','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square

% number of parameters
%--------------------------------------------------------------------------
PDP   = spm_RDP_MI(MDP);
np    = 0;
for n = 1:Nm
    for g = 1:numel(PDP{n}.a)
        np = np + nnz(PDP{n}.a{g});
    end
    for f = 1:numel(PDP{n}.b)
        np = np + nnz(PDP{n}.b{f});
    end
end
disp('number of parameters'), disp(np)

return


% NOTES
%==========================================================================
% NT    = 100;
% NS    = [];
% NU    = [];
% for q = 0:NT:(GDP.T - NT)
% 
%     % get attracting paths
%     %--------------------------------------------------------------------
%     S   = spm_get_sequences(MDP);
% 
%     % augment MDP
%     %--------------------------------------------------------------------
%     t   = (1:NT) + q;
%     MDP = spm_daisy_chain(PDP.O(:,t),S,MDP,GDP);
% 
%     disp(100*q/GDP.T)
%     NS(end + 1) = size(MDP{Nm}.b{1},2);
%     NU(end + 1) = size(MDP{Nm}.b{1},3);
% 
% end