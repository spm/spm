function MDP = DEM_AtariIII
% Structure learning from pixels
%__________________________________________________________________________
%
% This version of the Atari demos illustrates the use of multimodal streams
% in a renormalising generative model. The first stream derives from pixels
% reporting the state of the game that generates consequences in training
% streams; here, proprioception (e.g., paddle location) rewards and
% punishments. Fast structure learning is based on a long run of
% pseudorandom play, picking out episodes that end in a reward, which are
% appended or merged into the renormalising model (RGM). The second stage
% of structure learning engages purposeful behaviour to learn the
% transitions among rewarded episodes. Here, purposeful behaviour is
% induced using expected cost at each level of the RGM and inductive
% inference in the space of generalised states (i.e., episodes at the
% highest level. Finally, the mapping from the leading stream (pixels) to
% training streams is used to merge states at the highest level to compress
% or reduce the model in a way that preserves the mutual information (i.e.,
% expected free energy) between cause (i.e., first stream episodes) and
% consequences (i.e., policies, rewards and cost).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

% set up and preliminaries
%--------------------------------------------------------------------------
rng(1)

% Get game: i.e., generative process (as a partially observed MDP)
%==========================================================================
Nr = 12;                                   % number of rows
Nc = 9;                                    % number of columns
Sc = 16;                                   % Spatial scaling

% specify game
%--------------------------------------------------------------------------
G  = @(Nr,Nc) spm_MDP_breakout(Nr,Nc);     % game
G  = @(Nr,Nc) spm_MDP_pong(Nr,Nc,4,true);  % game

[GDP,~,~,~,RGB] = G(Nr,Nc);                % get game in MDP form

% size of streams
%--------------------------------------------------------------------------
S      = ones(4,3);
S(1,:) = [Nr,Nc,1];                        % sensory stream
S(2,:) = [1 1 1];                          % reward  stream                 
S(3,:) = [1 1 1];                          % cost    stream
S(4,:) = [1 1 1];                          % policy  stream

% reward and cost functions [of outcomes]
%--------------------------------------------------------------------------
spm_get_hits = @(o,id) find(o(id.reward,:) > 1);
spm_get_miss = @(o,id) find(o(id.contraint,:) > 1);

% Generate (probabilistic) outcomes under random actions
%==========================================================================
spm_figure('GetWin','Gameplay'); clf

GDP.tau = 2;                               % smoothness of random paths
GDP.T   = 5000;                            % (half) training length
PDP     = spm_MDP_generate(GDP);           % generate play

% find rewarded and costly events
%--------------------------------------------------------------------------
r     = spm_get_hits(PDP.o,GDP.id);
c     = spm_get_miss(PDP.o,GDP.id);
MDP   = {};
for i = 1:numel(r)

    % for each sequence ending with a reward
    %----------------------------------------------------------------------
    Nm = numel(MDP);
    Ne = max(2^(Nm - 1),1);
    s  = c(find(c < r(i),1,'last'));
    t  = (s + Ne):(r(i) + Ne);
    if numel(t)
       
        % assimilate this sequence
        %------------------------------------------------------------------
        O = PDP.O(:,t);
        if Nm
            MDP = spm_merge_structure_learning(O,MDP);
        else
            MDP = spm_faster_structure_learning(O,S,Sc);
        end

        % add delayed sequences
        %------------------------------------------------------------------
        for t = 2:Ne
            MDP = spm_merge_structure_learning(O(:,t:end),MDP);
        end

        % illustrate sequence
        %------------------------------------------------------------------
        subplot(2,1,1)
        for t = 1:size(O,2)
            imshow(spm_O2rgb(O(:,t),RGB)), drawnow
        end
    end
end


% Generate play from recursive generative model
%==========================================================================

% rewarded events
%--------------------------------------------------------------------------
C     = 32;
MDP   = spm_set_goals(MDP,[2,3], [C,-C]);
MDP   = spm_set_costs(MDP,[2,3], [C,-C]);
RDP   = spm_mdp2rdp(MDP);
RDP.T = 32;
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Generative AI'); clf
spm_show_RGB(PDP,RGB);

% and hits
%--------------------------------------------------------------------------
h   = spm_get_hits(PDP.Q.o{1},GDP.id);
subplot(Nm + 3,2,2*(Nm + 1))
plot(h,zeros(size(h)),'.r','MarkerSize',16)
drawnow


%% Active inference
%==========================================================================
% In what follows, we engage inductive planning as inference, with explicit
% action; namely, the generative process is used to engage actions that
% reproduce predicted outcomes (here, bat position) at the lowest level.
%--------------------------------------------------------------------------

% create hierarchical model with prior concentration parameters
%--------------------------------------------------------------------------
MDP{1}.GA  = GDP.A;
MDP{1}.GB  = GDP.B;
MDP{1}.GU  = GDP.U;
MDP{1}.GD  = GDP.D;
MDP{1}.ID  = GDP.id;

MDP{1}.chi = 512;                          % sticky action/shaky hand

% enable active learning (with minimal forgetting)
%--------------------------------------------------------------------------
for m = 1:numel(MDP)
    MDP{m}.beta = 4;
    MDP{m}.eta  = 512;
end

% enable learning of transition parameters to induce exploration
%--------------------------------------------------------------------------
FIX.A = 1;                                 % fix likelihood
FIX.B = 0;                                 % learn transitions
NT    = 200;                               % length of each game
NR    = fix(GDP.T/NT);                     % number of replications
F     = NaN(6,NR);
for i = 1:NR

    % assemble RGM
    %----------------------------------------------------------------------
    MDP   = spm_set_goals(MDP,[2,3],[C,-C]);
    MDP   = spm_set_costs(MDP,[2,3],[C,-C]);
    RDP   = spm_mdp2rdp(MDP,0,1/512,2,FIX);
    RDP.T = NT/Ne;

    % play
    %----------------------------------------------------------------------
    PDP   = spm_MDP_VB_XXX(RDP);
    h     = spm_get_hits(PDP.Q.o{1},GDP.id);

    if true

        % Illustrate recursive model
        %------------------------------------------------------------------
        spm_figure('GetWin','Active inference'); clf
        spm_show_RGB(PDP,RGB,4,0);

        % and hits
        %------------------------------------------------------------------
        subplot(Nm + 3,2,2*(Nm + 1))
        plot(h,zeros(size(h)),'.r','MarkerSize',16)
        drawnow

    end

    % record behaviour
    %----------------------------------------------------------------------
    F(1,i) = mean(PDP.F);                 % ELBO (last level - states)
    F(2,i) = mean(PDP.Z);                 % ELBO (last level - paths)

    F(3,i) = PDP.Q.F;                     % ELBO (lower levels)
    F(4,i) = size(PDP.B{1},2);            % number of (last level) states
    F(5,i) = size(PDP.B{1},3);            % number of (last level) paths
    F(6,i) = numel(h);                    % number of hits
    
    % learn from mistakes
    %======================================================================
    R     = PDP.Q.O{1}(:,(Ne*Ne):end);
    for q = 1:2

        % new attracting paths
        %------------------------------------------------------------------
        S   = spm_get_sequences(MDP);

        % augment MDP
        %------------------------------------------------------------------
        MDP = spm_daisy_chain(R,S,MDP,GDP);

    end

    % results
    %----------------------------------------------------------------------
    spm_figure('GetWin','Structure learning');
    str = {'ELBO - states',...
        'ELBO - paths',...
        'ELBO - both',...
        'Latent states',...
        'Latent paths',...
        'Reward count'};
    for f = 1:numel(str)
        subplot(4,2,f)
        plot(F(f,:)), axis square
        title(str{f},'FontSize',14), xlabel('games')
    end

end

% Illustrate in latent state space 
%-=========================================================================

% assemble RGM
%--------------------------------------------------------------------------
MDP   = spm_RDP_sort(MDP);
MDP   = spm_set_goals(MDP,[2,3],[C,-C]);
MDP   = spm_set_costs(MDP,[2,3],[C,-C]);
RDP   = spm_mdp2rdp(MDP,0,1/512);
RDP.T = 128;

% play
%--------------------------------------------------------------------------
PDP   = spm_MDP_VB_XXX(RDP);

spm_figure('GetWin','Active inference (before compression)'); clf
spm_show_RGB(PDP,RGB,8,false);

% and hits
%--------------------------------------------------------------------------
subplot(Nm + 3,2,2*(Nm + 1))
h   = spm_get_hits(PDP.Q.o{1},GDP.id);
plot(h,zeros(size(h)),'.r','MarkerSize',16)
drawnow

% Illustrate in latent state space 
%-=========================================================================
spm_figure('GetWin','Orbits'); clf

subplot(2,2,1)
HID   = PDP.id.hid;
spm_dir_orbits3(PDP.B{1},HID);

% paths to hits
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
plot(HID,zeros(size(HID)) + 1/2,'.r','MarkerSize',16), hold off
title('Paths to hits (before)','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square


% Bayesian model reduction: illustrate play under a reduced model
%==========================================================================

% assemble RGM
%--------------------------------------------------------------------------
RDP   = spm_RDP_MI(MDP);                   % merge conserving MI 
RDP   = spm_RDP_sort(RDP);                 % sort states by NESS density
RDP   = spm_set_goals(RDP,[2,3],[C,-C]);   % set intended states (h/cid)
RDP   = spm_set_costs(RDP,[2,3],[C,-C]);   % set contraints (C)
RDP   = spm_mdp2rdp(RDP);                  % get nested model

B        = RDP.B{1}*64 + 1;                % reduce prior precision
RDP.B{1} = spm_dir_norm(B);                % and normalise

RDP.T = 128;                               % length of each game                          
NR    = 8;                                 % number of replications
F     = NaN(4,NR);
for i = 1:NR

    % active inference and learning
    %----------------------------------------------------------------------
    PDP  = spm_MDP_VB_XXX(RDP);

    % Illustrate recursive model
    %----------------------------------------------------------------------
    spm_figure('GetWin','Active inference (with compression)'); clf
    spm_show_RGB(PDP,RGB,4,0);

    % and hits
    %----------------------------------------------------------------------
    h = spm_get_hits(PDP.Q.o{1},GDP.id);
    subplot(Nm + 3,2,2*(Nm + 1))
    plot(h,zeros(size(h)) - 2,'.r','MarkerSize',16)
    drawnow

    % record behaviour
    %----------------------------------------------------------------------
    F(1,i) = PDP.Q.F;                      % ELBO (lower levels)
    F(2,i) = numel(h);                     % number of hits
 
end

% results
%--------------------------------------------------------------------------
spm_figure('GetWin','Performance');
str = {'ELBO','Number of rewards'};
for f = 1:numel(str)
    subplot(4,2,f)
    plot(F(f,:)), axis square
    title(str{f},'FontSize',14), xlabel('games')
end

% Illustrate in latent state space 
%-=========================================================================
spm_figure('GetWin','Orbits');

subplot(2,2,2)
HID   = RDP.id.hid;
spm_dir_orbits3(RDP.B{1},HID);

% paths to hits
%--------------------------------------------------------------------------
subplot(2,2,4)
B     = sum(RDP.B{1},3) > 1/32;
Ns    = size(B,1);
h     = sparse(1,HID,1,1,Ns);
I     = [];
for t = 1:32
    I(t,:) = h;
    h      = (h + h*B) > 0;
end
imagesc(I), hold on 
plot(HID,zeros(size(HID)) + 1/2,'.r','MarkerSize',16), hold off
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

