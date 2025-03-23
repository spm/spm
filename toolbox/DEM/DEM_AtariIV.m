function MDP = DEM_AtariIV
% de novo Structure learning — of pullback attractors —  afrom pixels 
%__________________________________________________________________________
%
% This routine illustrates the de novo learning of a state action policy —
% and goal-directed behaviour using inductive inference under the ensuing
% model. The approach is based upon the fast structure learning of dynamics
% (i.e., transitions among generalised or deep states) that possess a
% nonequilibrium steady-state distribution (NESS). This is guaranteed when
% all states are transient states on a (pullback) attractor. In this
% example, random play is generated for 10,000 frames of a simplified
% arcade game. These are used to learn the structure of a renormalising
% generative model; namely, the recursive composition of states and paths
% specified by edges or dependencies (i.e., indices of parents). The
% resulting structure is then used to accumulate sequences of random play
% that, by chance, include rewarded outcomes. Following this active
% selection the training data are used to accumulate
% transitions to, and only to, the basins of attraction of goal states.
% When all goal states lie in a basin of attraction — of other goal states
% — the accumulation terminates. At this point, the dynamics possess a
% nonequilibrium steady-state, which is used to reduce the model by
% pruning away redundant states.
%
% This demonstration features an efficient self-[Eigen]reduction, using
% the nonequilibrium steady-state distribution to evaluate the expected
% reward. Generalised states are removed if the (analytic) gradient of
% expected cost — with respect to the probability of transitions to the
% state in question — is positive.
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
rng(1)

EIGEN     = true;                            % eigenreduction
INDUCTION = true;                            % inductive inference

% Get game: i.e., generative process (as a partially observed MDP)
%==========================================================================
Nr = 12;                                     % number of rows
Nc = 9;                                      % number of columns
Sc = 32;                                     % spatial scaling
Nd = 1;                                      % random initial conditions
C  = 32;                                     % log cost
T  = 2;                                      % time scaling
L  = [16,4];                                 % path lengths

% get game in MDP form
%--------------------------------------------------------------------------
[GDP,~,~,~,RGB] = spm_MDP_arcade(Nr,Nc,Nd,true);

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

GDP.tau = 1;                                  % smoothness of random paths
GDP.T   = 10000;                              % training length
PDP     = spm_MDP_generate(GDP);              % generate play

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

% structure learning: grouping operators (iA,iB,iC,...)
%==========================================================================
MDP = spm_faster_structure_learning(PDP.O,S,Sc,T);

% forget parameters
%--------------------------------------------------------------------------
Nm  = numel(MDP);
for n = 1:Nm
    for g = 1:numel(MDP{n}.a)
        MDP{n}.a{g} = [];
    end
    for f = 1:numel(MDP{n}.b)
        MDP{n}.b{f} = [];
    end
end

% Parameter learning: basins of attraction of goal states
%==========================================================================
spm_figure('GetWin','Discovery'); clf

NT    = 100;                         % number of frames per game 
NS    = [];                          % number of states  
NU    = [];                          % number of paths 
NA    = [];                          % number of childless states
NO    = [];                          % number of orphan states
NH    = [];                          % number of goal states
NG    = [];                          % number of goal states with children
NF    = [];                          % number of goal states in orbit
GDP.T = NT;
for i = 1:200

    % Accumulate states under random play
    %======================================================================
    if i < 1

        % reuse continuous training data
        %------------------------------------------------------------------
        t = (0:(NT + T)) + rem(i,100 - 1)*NT;
        
    else

        % regenerate training data from initial conditions
        %------------------------------------------------------------------
        t = (0:(NT - T));
        for j = 1:128
            PDP = spm_MDP_generate(GDP);

            % break, if there is a least one hit
            %--------------------------------------------------------------
            if any(spm_get_hits(PDP.o,GDP.id))
                break
            end
        end
    end


    % Accumulate these states
    %----------------------------------------------------------------------
    for s = 1:T
        MDP = spm_merge_structure_learning(PDP.O(:,t + s),MDP);
    end

    % Retain states in basin of attraction to goal states
    %----------------------------------------------------------------------
    [MDP,d,o,h,c] = spm_RDP_basin(MDP,[2,3],[C,-C],L);
    [r,R,g,f]     = spm_RDP_radius(MDP,h,c);

    % report implcit physics discovery
    %======================================================================
    NS(end + 1) = size(MDP{Nm}.b{1},2);    
    NU(end + 1) = size(MDP{Nm}.b{1},3);
    NA(end + 1) = sum(~d);
    NO(end + 1) = sum(~o);
    NH(end + 1) = numel(h);
    NF(end + 1) = numel(f);
    NG(end + 1) = sum(g);

    subplot(4,2,1), plot(NS), title('Generalised states'), axis square
    subplot(4,2,2), plot(NU), title('Generalised paths'),  axis square
    subplot(4,2,3), plot(NA), title('Childless states'),   axis square, hold on
    subplot(4,2,3), plot(NO), title('Childless states'),   axis square, hold off
 
    legend({'childless','orphans'},'Location','northwest','Box','off')

    subplot(4,2,4), plot(NH), title('Goal states'), axis square, hold on
    subplot(4,2,4), plot(NG), title('Goal states'), axis square, hold on
    subplot(4,2,4), plot(NF), title('Goal states'), axis square, hold off

    legend({'goal states','with children','on attractor'},'Location','northwest','Box','off')

    [~,j] = sort(sum(R,1),'descend');
    subplot(4,2,5); spy(sum(MDP{Nm}.b{1},3),'.k')
    title('Transition discovery'), xlabel('states'), axis square
    subplot(4,2,6); imagesc(1 - R(j,j))
    title('Goal adjacency'),       xlabel('games'),  axis square
    subplot(4,1,4); plot(sort(r,'descend'))
    title('Radius of goals'),      xlabel('goals'),  axis square
    drawnow
    
end

% model reduction
%==========================================================================

% prune childless states
%--------------------------------------------------------------------------
MDP   = spm_RDP_prune(MDP);
MDP   = spm_set_goals(MDP,[2,3],[C,-C]);
hid   = MDP{end}.id.hid;
cid   = MDP{end}.id.cid;

% Eigenreduction to maximize expected reward under NESS density
%--------------------------------------------------------------------------
if EIGEN
    for i = 1:2
        MDP = spm_RDP_ness(MDP,hid,cid);
        MDP = spm_set_goals(MDP,[2,3],[C,-C]);
        hid = MDP{end}.id.hid;
        cid = MDP{end}.id.cid;
    end
end

% illustrate pruning
%--------------------------------------------------------------------------
spm_figure('GetWin','Attractors'); clf
subplot(2,2,1); imagesc(1 - R(j,j)), title('Goal adjacencies'), axis square

[~,R] = spm_RDP_radius(MDP,hid,cid);
[~,j] = sort(sum(R,1),'descend');

subplot(2,2,2); imagesc(1 - R(j,j)), title('After pruning'),    axis square


% Illustrate transititions and paths to goals
%-=========================================================================
subplot(2,2,3)
spm_dir_orbits(MDP{Nm}.b{1},hid,64);

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
plot(hid,zeros(size(hid)) + 1/2,'.g','MarkerSize',16), hold off
title('Paths to goals','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square
drawnow

% Illustrate generative (i.e., fictive) play
%==========================================================================

% assemble RGM
%--------------------------------------------------------------------------
RDP   = spm_set_goals(MDP,[2,3],[C,-C]);   % set intended states (hid/cid)
RDP   = spm_set_costs(RDP,[2,3],[C,-C]);   % set contraints (C)
RDP   = spm_mdp2rdp(RDP);                  % get nested model
RDP.T = fix(512/T);
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Generative AI'); clf
spm_show_RGB(PDP,RGB);

% and hits
%--------------------------------------------------------------------------
subplot(Nm + 3,2,2*(Nm + 1))
h = spm_get_hits(PDP.Q.o{1},GDP.id);
c = spm_get_miss(PDP.Q.o{1},GDP.id);
plot(h,zeros(size(h)) - 2,'.g','MarkerSize',16)
plot(c,zeros(size(c)) - 2,'.r','MarkerSize',16)
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
MDP{1}.GE  = GDP.E;

MDP{1}.ID  = GDP.id;
MDP{1}.chi = 512;                           % sticky action/shaky hand
NS         = 512;                           % concentration parameter

% Continual learning
%--------------------------------------------------------------------------
NT    = 200;                                 % length of each game
NR    = 100;                                 % number of replications
F     = NaN(5,NR);
t     = 0:(NT - T);
for i = 1:NR

    % set intended states and constraints
    %----------------------------------------------------------------------
    RDP   = spm_set_goals(MDP,[2,3],[C,-C]);
    RDP   = spm_set_costs(RDP,[2,3],[C,-C]);

    % assemble RGM and reduce prior precision
    %----------------------------------------------------------------------
    RDP   = spm_mdp2rdp(RDP,0,1/NS);

    % preclude inductive inference for illustrative purposes
    %----------------------------------------------------------------------
    if ~INDUCTION, RDP.U = spm_zeros(RDP.U); end

    % active inference
    %----------------------------------------------------------------------
    RDP.T = fix(NT/T);
    PDP   = spm_MDP_VB_XXX(RDP);
    
    % learn from mistakes
    %======================================================================

    % accumulate self-generated sequence (i.e., Dirichlet counts)
    %----------------------------------------------------------------------
    O     = PDP.Q.O{1};
    for s = 1:T
        MDP = spm_merge_structure_learning(O(:,t + s),MDP);
    end

    % retain states in basins of attraction
    %----------------------------------------------------------------------
    MDP   = spm_RDP_basin(MDP,[2,3],[C,-C],L);


    % Illustrate intentional play
    %======================================================================
    spm_figure('GetWin','Active inference'); clf
    spm_show_RGB(PDP,RGB,4,0);

    % and hits
    %----------------------------------------------------------------------
    subplot(Nm + 3,2,2*(Nm + 1))
    h = spm_get_hits(PDP.Q.o{1},GDP.id);
    c = spm_get_miss(PDP.Q.o{1},GDP.id);
    plot(h,zeros(size(h)) - 2,'.g','MarkerSize',16)
    plot(c,zeros(size(c)) - 2,'.r','MarkerSize',16)
    drawnow

    % record structure learning and behaviour
    %----------------------------------------------------------------------
    F(1,i) = size(PDP.B{1},2);                  % number of (deep) states
    F(2,i) = size(PDP.B{1},3);                  % number of (deep) paths
    F(3,i) = PDP.Q.F + sum(PDP.F);              % ELBO
    F(4,i) = numel(h);                          % number of hits
    F(5,i) = numel(c);                          % number of misses

    % plot results
    %----------------------------------------------------------------------
    spm_figure('GetWin','Structure learning');
    str = {...
        'Latent states',...
        'Latent paths',...
        'ELBO',...
        'Rewards and losses'};
    for f = 1:numel(str)
        subplot(3,2,f)
        plot(F(f,:)), axis square
        title(str{f},'FontSize',14), xlabel('games')
    end
    plot((1:i)*NT,F(4,1:i)), hold on
    plot((1:i)*NT,F(5,1:i)), hold off
    legend({'rewards','losses'},'Location','northwest','Box','off')
    title(str{4},'FontSize',14), xlabel('frames'), axis square
    drawnow

end

% Illustrate performance : before compression
%-=========================================================================

% assemble RGM
%--------------------------------------------------------------------------
RDP   = spm_set_goals(MDP,[2,3],[C,-C]);
RDP   = spm_set_costs(RDP,[2,3],[C,-C]);
RDP   = spm_mdp2rdp(RDP,0,1/NS);

% play
%--------------------------------------------------------------------------
RDP.T = fix(512/T);
PDP   = spm_MDP_VB_XXX(RDP);

spm_figure('GetWin','Active inference (before compression)'); clf
spm_show_RGB(PDP,RGB,8,false);

% and hits
%--------------------------------------------------------------------------
subplot(Nm + 3,2,2*(Nm + 1))
h = spm_get_hits(PDP.Q.o{1},GDP.id);
c = spm_get_miss(PDP.Q.o{1},GDP.id);
plot(h,zeros(size(h)) - 2,'.g','MarkerSize',16)
plot(c,zeros(size(c)) - 2,'.r','MarkerSize',16)
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
plot(HID,zeros(size(HID)) + 1/2,'.g','MarkerSize',16), hold off
title('Paths to goals (before)','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square


% Illustrate performance : after compression
%-=========================================================================

% assemble RGM
%--------------------------------------------------------------------------
RDP   = spm_RDP_MI(MDP);                   % lossless merge 
RDP   = spm_set_goals(RDP,[2,3],[C,-C]);   % set intended states (hid/cid)
RDP   = spm_set_costs(RDP,[2,3],[C,-C]);   % set contraints (C)
RDP   = spm_mdp2rdp(RDP,0,1/NS);           % get nested model

% active inference and learning
%--------------------------------------------------------------------------
RDP.T = fix(512/T);
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference (with compression)'); clf
spm_show_RGB(PDP,RGB,4,0);

% and hits
%--------------------------------------------------------------------------
subplot(Nm + 3,2,2*(Nm + 1))
h = spm_get_hits(PDP.Q.o{1},GDP.id);
c = spm_get_miss(PDP.Q.o{1},GDP.id);
plot(h,zeros(size(h)) - 2,'.g','MarkerSize',16)
plot(c,zeros(size(c)) - 2,'.r','MarkerSize',16)
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
plot(HID,zeros(size(HID)) + 1/2,'.g','MarkerSize',16), hold off
title('Paths to goals (after)','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square

% number of parameters
%--------------------------------------------------------------------------
np    = 0;
for n = 1:Nm
    for g = 1:numel(MDP{n}.a)
        np = np + nnz(MDP{n}.a{g});
    end
    for f = 1:numel(MDP{n}.b)
        np = np + nnz(MDP{n}.b{f});
    end
end
disp('number of parameters'), disp(np)

return





% UTILITY CODE
%==========================================================================
%     % attracting paths
%     %--------------------------------------------------------------------
%     S   = spm_get_sequences(MDP);
%     for j = 1:numel(h)
%         imshow(spm_O2rgb(S{:,h(j)},RGB))
%         drawnow,pause
%     end

% structure learning from self-generated play
%==========================================================================
if SELF

    % forget parameters
    %----------------------------------------------------------------------
    for n = 1:Nm
        for g = 1:numel(MDP{n}.a)
            MDP{n}.a{g} = [];
        end
        for f = 1:numel(MDP{n}.b)
            MDP{n}.b{f} = [];
        end
    end

    % self-teaching
    %----------------------------------------------------------------------
    NT    = 200;                                    % length of each game
    NR    = 100;                                    % number of games
    F     = NaN(5,NR);
    t     = 0:(NT - T);
    RDP.T = fix(NT/T);
    for i = 1:NR

        % accumulate self-generated sequence (i.e., Dirichlet counts)
        %------------------------------------------------------------------
        PDP   = spm_MDP_VB_XXX(RDP);
        O     = PDP.Q.O{1};
        for s = 1:T
            MDP = spm_merge_structure_learning(O(:,t + s),MDP);
        end

        % retain states in basins of attraction
        %------------------------------------------------------------------
        MDP   = spm_RDP_basin(MDP,[2,3],[C,-C],L);

        % record structure learning and behaviour
        %------------------------------------------------------------------
        h      = spm_get_hits(PDP.Q.o{1},GDP.id);
        c      = spm_get_miss(PDP.Q.o{1},GDP.id);
        F(1,i) = size(MDP{end}.b{1},2);             % number of (deep) states
        F(2,i) = size(MDP{end}.b{1},3);             % number of (deep) paths
        F(3,i) = PDP.Q.F + sum(PDP.F);              % ELBO
        F(4,i) = numel(h);                          % number of hits
        F(5,i) = numel(c);                          % number of misses

        % plot results
        %------------------------------------------------------------------
        spm_figure('GetWin','Self-teaching');
        str = {...
            'Latent states',...
            'Latent paths',...
            'ELBO',...
            'Rewards and losses'};
        for f = 1:numel(str)
            subplot(3,2,f)
            plot(F(f,:)), axis square
            title(str{f},'FontSize',14), xlabel('games')
        end
        plot(F(4,:)), hold on
        plot(F(5,:)), hold off
        legend({'rewards','losses'},'Location','northwest','Box','off')
        drawnow
    end
end

