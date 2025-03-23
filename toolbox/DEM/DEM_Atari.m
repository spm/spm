function MDP = DEM_Atari
% Structure learning from pixels
%__________________________________________________________________________
%
% This routine addresses the problem of learning a generative model from
% pixels, under the constraints supplied by sparse rewards. The problem is
% solved in a fast and frugal way using a structure learning approach based
% upon a generalised Markov decision process (that treats states and their
% paths as pairs of random variables, whose dynamics are encoded in
% transition tensors).
% 
% A key architectural aspect of these deep structures is an appeal to the
% renormalisation group; in the sense that there is a recursive application
% of grouping and reduction operators, as we move from one hierarchical
% level to the next. These operate over both space and time: in the sense
% that groups of pixels are jointly inferred over a short period of time;
% such that the state any given level of the model generates the
% instantaneous state and path (i.e., velocity) of a group of latent states
% at the lower level. In these examples, the group operator effectively
% tiles a lattice of latent states (and instantaneous paths) into little
% contiguous squares. This necessary induces a separation of temporal
% scales such as the states at the highest level see all the pixels at the
% lowest level and, effectively,, encode a trajectory over to the 2^n
% timesteps, where n is the depth of the model.
% 
% Note that in this illustration, there is no learning of model parameters.
% The structure is identified using a simple form of supervised structure
% learning. This simple form conforms to the approach in machine learning,
% in which random play is sporadically rewarded. By selecting episodes of
% play (i.e., exemplar training data) that end with a reward, one can
% emulate reinforcement learning with active model selection. This can be
% read as replacing a scheme in which 'everything is presented but rewarded
% sequences are learned' with the equivalent protocol in which 'rewarded
% sequences are presented and everything is learned'. The resulting
% structure is then sufficient to reproduce expert play; because this is
% the kind of play that the model knows.
% 
% One can finesse expert play by the use of inductive inference; provided
% initial conditions or novel trajectories to be recognised as such, and
% actively reconfigured to recognisable (expert trajectories) through
% explicit action. Explicit action refers to (merely) reflexive active
% inference; namely, selecting those actions that minimise free energy that
% produce outcomes that match predicted outcomes. Here, the outcomes are
% pixel-based outcomes.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%--------------------------------------------------------------------------
rng(1)

% Get game: i.e., generative process (as a partially observed MDP)
%==========================================================================
Nr = 12;                               % number of rows
Nc = 9;                                % number of columns
G  = @(Nr,Nc) spm_MDP_breakout(Nr,Nc); % game
NT = 2048;                             % exposures (training)
G  = @(Nr,Nc) spm_MDP_pong(Nr,Nc);     % game
NT = 256;                              % exposures (training)



%% Simulate learning and subsequent performance
%--------------------------------------------------------------------------
[GDP,hid,cid,con,RGB] = G(Nr,Nc); 

% generate (probabilistic) outcomes under random actions
%==========================================================================
spm_figure('GetWin','Gameplay'); clf

GDP.tau = 2;                          % smoothness of random paths
GDP.T   = (Nr + 2)*4;                 % uppr bound on length of paths

O     = {};
for n = 1:(NT*32)

    % initialise this batch of training exemplars
    %----------------------------------------------------------------------
    PDP = spm_MDP_generate(GDP);

    % smart data selection
    %----------------------------------------------------------------------
    t   = find(ismember(PDP.s',hid','rows'),1,'first');
    c   = find(ismember(PDP.s',cid','rows'),1,'first');
    c   = min([c; GDP.T]);

    % Maxwell's Daemon
    %----------------------------------------------------------------------
    ACCEPT = numel(t) && t < GDP.T && c > t;

    if ACCEPT

        % accumulate (rewarded) sequences
        %------------------------------------------------------------------
        if numel(O)
            O = [O PDP.O(:,1:t)];
            h(end + 1) = h(end) + t;
        else
            O = PDP.O(:,1:t);
            h = t;
        end

        % start after we left off
        %------------------------------------------------------------------
        GDP.s = PDP.s(:,t + 1);
        GDP.u = PDP.u(:,t + 1);

        % illustrate outcomes
        %------------------------------------------------------------------
        if size(O,2) < 256
            subplot(2,1,1)
            for i = 1:t
                imshow(spm_O2rgb(PDP.O(:,i),RGB)), drawnow
            end
        end

    end

    % break if a sufficient number of episodes have been accumulated
    %----------------------------------------------------------------------
    clc; fprintf('Number of samples %i (%i)\n',size(O,2),(n*GDP.T))
    if size(O,2) > NT
        O = O(:,1:NT);
        break
    end

end

% illustrate intial outcomes
%--------------------------------------------------------------------------
spm_figure('GetWin','Gameplay'); clf
for t = 1:4
    subplot(4,4,t)
    imshow(spm_O2rgb(O(:,t),RGB)), drawnow
end

% illustrate orbits
%--------------------------------------------------------------------------
subplot(3,2,3)
u     = spm_svd(spm_cat(O'));
plot(u(:,1),u(:,2),'c'), hold on
plot(u(:,1),u(:,2),'ob','MarkerSize',4)
plot(u(h,1),u(h,2),'or','MarkerSize',8), hold off
title('Orbits or paths'), xlabel('first'), ylabel('second'), axis square

subplot(3,2,4)
u     = spm_svd(spm_cat(O'));
plot(u(:,2),u(:,3),'c'), hold on
plot(u(:,2),u(:,3),'ob','MarkerSize',4)
plot(u(h,2),u(h,3),'or','MarkerSize',8), hold off
title('Orbits or paths'), xlabel('second'), ylabel('third'), axis square

subplot(3,2,5)
plot(u(:,1),u(:,2),'c'), hold on
plot(u(:,1),u(:,2),'ob','MarkerSize',4)
plot(u(h,1),u(h,2),'or','MarkerSize',8), hold off
title('Orbits or paths'), xlabel('first'), ylabel('second'), axis square
axis([-1 1 -1 1]/1e3)
drawnow


% RG structure learning
%==========================================================================
tic, MDP = spm_faster_structure_learning(O,[Nr,Nc],8); toc

% rewarded events
%--------------------------------------------------------------------------
[HID,CID,HITS] = spm_get_episodes(hid,cid,GDP,MDP);

% Illustrate orbits
%==========================================================================
spm_figure('GetWin','Orbits'); clf

subplot(2,2,1)
spm_dir_orbits(MDP{end}.b{1},HID,24);
title('Orbits and goals')

% priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Priors'); clf
spm_RDP_params(MDP)

% Illustrate episodic dynamics (at deepest level)
%==========================================================================
Nm    = numel(MDP);
spm_figure('GetWin',sprintf('Paramters: level %i',Nm)); clf
spm_MDP_params(MDP{Nm})

% paths to hits
%--------------------------------------------------------------------------
subplot(2,2,3)
B     = sum(MDP{Nm}.b{1},3) > 0;
Ns    = size(B,1);
h     = sparse(1,HID,1,1,Ns);
for t = 1:16
    I(t,:) = h;
    h      = (h + h*B) > 0;
end

imagesc(I), hold on 
plot(HID,0*HID + 1/2,'.r','MarkerSize',16), hold off
title('Paths to hits','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square


% Generate play by sampling from the resulting deep generative model
%==========================================================================
spm_figure('GetWin','Generative AI - I'); clf

% sample outcomes
%--------------------------------------------------------------------------
pdp   = MDP{end};
pdp.T = 32;
pdp.s = 1;
pdp.u = 1;
pdp   = spm_MDP_VB_XXX(pdp);

Q     = cell(Nm,1);
Q{Nm} = pdp.O;
for n = Nm:-1:2

    % set empirical priors over states and paths
    %----------------------------------------------------------------------
    for t = 1:size(Q{n},2)
        pdp   = MDP{n - 1};
        pdp.T = 2;
        for g = 1:numel(pdp.id.D)
            try
                pdp.D{g} = Q{n}{pdp.id.D{g},t};
                pdp.E{g} = Q{n}{pdp.id.E{g},t};
            catch
                pdp.D{g} = ones(size(pdp.b{g},2),1);
                pdp.E{g} = ones(size(pdp.b{g},3),1);
            end
        end
        pdp      = spm_MDP_VB_XXX(pdp);
        Q{n - 1} = [Q{n - 1} pdp.O];
    end

end

spm_show_outcomes(Q,Nr,Nc)

% Generate play from recursive generative model
%==========================================================================

% Create deep recursive model
%--------------------------------------------------------------------------
RDP       = spm_mdp2rdp(MDP);
[~,Ns,Nu] = spm_MDP_size(RDP);

RDP.T    = 8;
RDP.D{1} = sparse(1,1,1,Ns(1),1);
RDP.E{1} = sparse(1,1,1,Nu(1),1);
PDP      = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Generative AI - II'); clf
spm_show_RGB(PDP,RGB);


% Now repeat but engage action and active inference
%==========================================================================
% The above illustrations are simply generating outcomes from a
% hierarchical or recursive generative model. In what follows, we engage
% inductive planning as inference, with explicit action; namely, the
% generative process is used to engage actions that reproduce predicted
% outcomes (here, purely visual) at the lowest level.
%--------------------------------------------------------------------------

%% create hierarchical model with prior concentration parameters
%--------------------------------------------------------------------------
MDP{1}.GA  = GDP.A;
MDP{1}.GB  = GDP.B;
MDP{1}.GD  = GDP.D;
MDP{1}.GE  = GDP.E;
MDP{1}.GU  = GDP.U;

for g = 1:numel(GDP.A)
    MDP{1}.ID.A{g} = 1:(ndims(GDP.A{g}) - 1);
end

MDP{1}.ID.control = con;                    % controlled outcomes
MDP{1}.chi = 1;                             % sticky action/shaky hand


% enable active learning (with minimal forgetting)
%--------------------------------------------------------------------------
for m = 1:numel(MDP)
    MDP{m}.beta = 4;
    MDP{m}.eta  = 512;
end

% enable inductive inference at final level
%--------------------------------------------------------------------------
MDP{Nm}.U      = 1;
MDP{Nm}.id.hid = HID;
MDP{Nm}.id.cid = CID;

% train: with small concentration parameters
%--------------------------------------------------------------------------
FIX.A      = 0;                             % learn likelihood
FIX.B      = 0;                             % but not transitions
RDP        = spm_mdp2rdp(MDP,1/128,1/512,2,FIX);
RDP.T      = 128;
PDP        = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference'); clf
spm_show_RGB(PDP,RGB,4,false);

% and hits
%--------------------------------------------------------------------------
subplot(Nm + 3,2,2*(Nm + 1))
h     = find(ismember(PDP.Q.o{1}',HITS','rows'));
plot(h,zeros(size(h)),'.r','MarkerSize',16)


% Illustrate in latent state space 
%-=========================================================================
spm_figure('GetWin','Orbits');

subplot(2,2,2)
spm_dir_orbits(PDP.B{1},HID,24);

return



% subroutines
%==========================================================================

function RGB = spm_colour(O)
% subfunction: returns an RGB rendering of a multinomial distribution
%--------------------------------------------------------------------------
MAP = [1 0 1;
       0 1 0;
       0 0 1;
       1 1 0;
       1 1 1;
       1 1 0;
       0 1 1;
       1 0 1];
MAP = MAP(1:numel(O),:)';
RGB = MAP*O;

return

function spm_show_outcomes(O,Nr,Nc,Y)
% Plots the inferred sequence of moves
% FORMAT spm_show_outcomes(O,Nr,Nc,[Y])
% O{n}   - Cell array of probabilistic, hierarchal outcomes
% Nr     - Number of rows (x dimension)
% Nc     - Number of columns (x dimension)
% Y{1}   - Predicted outcomes
%
%--------------------------------------------------------------------------
% This auxiliary routine plots the hierarchal outcomes from a deep
% generative model according to the grouping of outcomes specified by the
% shape of RG (i.e., cell array containing lists of outcomes for each
% group)
%
% The separation of temporal scales is illustrated by showing the outcomes
% as a movie saved as the user data in the current figure.
%__________________________________________________________________________


% show sequence of moves
%==========================================================================
T     = size(O{1},2);
Nn    = numel(O);

% Show predictive posteriors over hierarchical outcomes
%--------------------------------------------------------------------------
for n = 1:Nn
    subplot(Nn + 2,1,2 + n)
    imagesc(1 - spm_cat(O{n})), axis xy
    title(sprintf('Predictive posteriors at level %i',n),'FontSize',12)
end

% And repeat in image format using outcome groups
%--------------------------------------------------------------------------
for t = 1:T

    % observed outcomes
    %--------------------------------------------------------------
    X     = zeros(Nr*Nc,3);
    o     = O{1}(:,t);
    for j = 1:numel(o)
        X(j,:) = spm_colour(o{j});
    end

    % predicted outcomes
    %--------------------------------------------------------------
    if nargin > 4

        P     = zeros(Nr*Nc,3);
        o     = Y{1}(:,t);
        for j = 1:numel(o)
            P(j,:) = spm_colour(o{j});
        end

        % plot first level outcomes as a coloured image
        %----------------------------------------------------------
        subplot(Nn + 2,6,min(t,6))
        imagesc(reshape(X,Nr,Nc,3)), axis image, axis xy
        title(sprintf('Outcomes: t =  %i',t),'FontSize',12)

        % plot first level predictions as a coloured image
        %----------------------------------------------------------
        subplot(Nn + 2,6,min(t,6) + 6)
        imagesc(reshape(P,Nr,Nc,3)), axis image, axis xy
        title(sprintf('Predictions: t = %i',t),'FontSize',12)
        drawnow

    else

        % plot first level outcomes as a coloured image
        %----------------------------------------------------------
        subplot(Nn + 2,6,min(t,6))
        imagesc(reshape(X,Nr,Nc,3)), axis image, axis xy
        title(sprintf('Outcomes: t = %i',t),'FontSize',12)
        drawnow

    end

end

return


%% NOTES on formulation with multiple MPDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code illustrates the composition of a meta-model of MDP structures
% or agents using the same (renormalisation group) approach above. Here,
% the conditional independence of factors pertaining to each group of
% outputs is automatically assured by assigning each group to its own MDP.
%==========================================================================
mdp.id.hid = hid;                          % Intended (latent) states
mdp.T = 128;                               % duration of game
GDP   = spm_MDP_VB_XXX(mdp);

% illustrate action and action selection
%--------------------------------------------------------------------------
spm_figure('GetWin','Inductive inference'); clf
spm_report(GDP,Nr,Nc,hid)

% Create image array
%--------------------------------------------------------------------------
S{1}  = [Nr,Nc];
I     = GDP.O;
for t = 1:size(I,2)
    for m = 1:size(I,1)
        O{1}{m,t} = I(m,t);
    end
end

% Create hierarchy of models
%--------------------------------------------------------------------------
S{1}  = [Nr,Nc];
dx    = 3; 
dt    = 2;     
for n = 1:8

    % Grouping into a partition of outcomes
    %----------------------------------------------------------------------
    g     = spm_tile(S{n}(1),S{n}(2),dx);
    t     = spm_time(size(O{n},2),dt);

    mdp   = cell(size(g));
    X     = {};
    P     = {};
    N     = {};
    for m = 1:numel(g)

        % mdp = spm_MDP_structure_learning(mdp,O)
        %------------------------------------------------------------------
        mdp{m} = spm_structure(O{n}(g{m},:));

        % Invert model to recover sequences of initial states and paths
        %------------------------------------------------------------------
        for i = 1:numel(t)
            pdp   = mdp{m};
            pdp.T = numel(t{i});
            for j = 1:pdp.T
                o = [O{n}{g{m},t{i}(j)}];
                pdp.O(:,j) = o(:);
            end
            pdp   = spm_MDP_VB_XXX(pdp);

            % initial states and paths
            %--------------------------------------------------------------
            X{m,i} = pdp.X{1}(:,1);
            P{m,i} = pdp.P{1}(:,end);
            N{m,i} = {X{m,i}; P{m,i}};

        end
    end

    MDP{n}   = mdp;
    G{n}     = g;
    if numel(g) == 1
        break,
    else
        O{n + 1} = N;
        S{n + 1} = size(g);
    end

end

% Generate expert play by sampling from deep generative Meta-model
%==========================================================================

% sample outcomes
%--------------------------------------------------------------------------
pdp     = MDP{end}{1};
pdp.T   = 32;
pdp.s   = 1;
pdp.u   = 1;
pdp     = spm_MDP_VB_XXX(pdp);
Q{n}{1} = pdp.O;

for n = Nm:-1:2

    % set empirical priors over states and paths
    %----------------------------------------------------------------------
    for m = 1:numel(Q{n})
        for g = 1:numel(G{n}{m})
            Qm    = {};
            mm    = G{n}{m}(g);
            for t = 1:size(Q{n}{m},2)
                pdp     = MDP{n - 1}{mm};
                pdp.T   = dt;
                pdp.D   = Q{n}{m}(2*g - 1,t);
                pdp.E   = Q{n}{m}(2*g - 0,t);
                pdp     = spm_MDP_VB_XXX(pdp);
                Qm      = [Qm pdp.O];
            end
            Q{n - 1}{mm} = Qm;

        end
    end
end

spm_show_outcomes(Q,G,Nr,Nc)

return

