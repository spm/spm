function MDP = DEM_Atari_learning
% learning from pixels
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
% read as replacing a scheme in which “everything is presented but rewarded
% sequences are learned” with the equivalent protocol in which “rewarded
% sequences are presented and everything is learned”. The resulting
% structure is then sufficient to reproduce expert play; because this is
% the kind of play that the model knows.
% 
% One can finesse expert play by the use of inductive inference through
% explicit action. Explicit action refers to (merely) reflexive active
% inference; namely, selecting those actions that minimise free energy that
% produce outcomes that match predicted outcomes. Here, the outcomes are
% pixel-based outcomes.
%
% This example accumulates sequences of rewarded (random) play by appending
% them to outcomes generated through gameplay using inductive inference. In
% other words, structure learning proceeds on the basis of what the agent
% has previously learned enabling the composition of long sequences of
% rewarded behaviour.

%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%--------------------------------------------------------------------------
rng(1)

% Get game: i.e., generative process (as a partially observed MDP)
%==========================================================================

% number of rows and columns
%--------------------------------------------------------------------------
Nr = 12; Nc = 9; [GDP,hid,cid] = spm_MDP_breakout(Nr,Nc);
Nr = 12; Nc = 9; [GDP,hid,cid] = spm_MDP_pong(Nr,Nc);


%% Simulate learning and subsequent performance
%--------------------------------------------------------------------------

% generate (probabilistic) outcomes under random actions
%==========================================================================
[O,o] = spm_generate_next(GDP,hid,cid);

spm_figure('GetWin','Gameplay'); clf
spm_report(o,Nr,Nc), drawnow

% Accumulate training data (O) by appending (rewarded) random play to
% self-generated outcomes
%--------------------------------------------------------------------------
for i = 1:64

    % RG structure learning
    %----------------------------------------------------------------------
    tic, [MDP,RG]  = spm_fast_structure_learning(O,[Nr,Nc]); toc

    % Get rewarding episodes
    %----------------------------------------------------------------------
    [HID,CID,HITS] = spm_get_rewards(hid,cid,GDP,MDP);

    % Equip model with generative process and initialise states
    %======================================================================
    Nm        = numel(MDP);
    MDP{1}.GA = GDP.A;
    MDP{1}.GB = GDP.B;
    MDP{1}.GD = GDP.D;
    MDP{1}.GE = GDP.E;
    MDP{1}.GU = GDP.U;

    % create hierarchical model with small prior concentration parameters
    %----------------------------------------------------------------------
    RDP    = spm_mdp2rdp(MDP);
    T      = size(O,2);
    RDP.U  = 1;
    RDP.T  = ceil(T/(2^(Nm - 1))) + 16;

    % Specify rewarded states for inductive inference
    %----------------------------------------------------------------------
    RDP.id.hid = HID;

    % invert 
    %======================================================================
    PDP    = spm_MDP_VB_XXX(RDP);

    % Mutual information (expected free energy)
    %----------------------------------------------------------------------
    MIA(i) = spm_MDP_MI(RDP.A);
    MIB(i) = spm_MDP_MI(RDP.B{1}(:,:,1));
    LEN(i) = T;

    % break if MI decreases
    %----------------------------------------------------------------------
    if MIB(i) < max(MIB) && LEN(i) == max(LEN)
        break
    end

    % or append next training sequence
    %----------------------------------------------------------------------
    t     = find(ismember(PDP.Q.o{1}(:,1:T)',HITS','rows'),1,'last');
    GDP.s = PDP.Q.s{1}(:,t + 1);
    GDP.u = PDP.Q.u{1}(:,t + 1);
    [O,o] = spm_generate_next(GDP,hid,cid);
    O     = [PDP.Q.O{1}(:,1:t) O];
    o     = [PDP.Q.o{1}(:,1:t) o];

    % Illustrate structure
    %----------------------------------------------------------------------
    spm_figure('GetWin','Paramters'); clf
    spm_MDP_params(MDP{end})

    subplot(4,2,3), plot(LEN,MIA)
    title('Mutual information (likelihood)','FontSize',14)
    xlabel('number of samples (frames)')
    ylabel('Mutual information'), drawnow
    
    subplot(4,2,4), plot(LEN,MIB)
    title('Mutual information (priors)','FontSize',14)
    xlabel('number of samples (frames)')
    ylabel('Mutual information'), drawnow

end

% Illustrate belief updating
%--------------------------------------------------------------------------
spm_figure('GetWin','Inference'); clf
spm_MDP_VB_trial(PDP);

spm_figure('GetWin','Active inference'); clf
spm_show_outcomes(spm_get_O(PDP),RG,Nr,Nc,PDP.Q(1).Y)

h = find(ismember(PDP.Q.o{1}',HITS','rows'));
subplot(2*Nm,1,Nm + 1), hold on
plot(h,ones(size(h)),'.r','MarkerSize',32), drawnow

return



% subroutines
%==========================================================================

function [O,o] = spm_generate_next(GDP,hid,cid)
% generates selected outcomes
% FORMAT [O,o] = spm_generate_next(GDP,hid,cid)
%
% This auxiliary routine generates outcomes from the random play, from an
% initial state, until a reward is encountered (indexed by hid in latent
% state space) in the absence of any constraints (indexed by cid latent
% state space).
%__________________________________________________________________________

% sample outcomes from initial state
%--------------------------------------------------------------------------
o     = [];
O     = {};
GDP.T = 32;
for i = 1:1024

    % generate outputs
    %----------------------------------------------------------------------
    PDP   = spm_MDP_generate(GDP);

    % smart data selection
    %----------------------------------------------------------------------
    t = find(ismember(PDP.s',hid','rows'),1,'first');
    c = find(ismember(PDP.s',cid','rows'),1,'first');
    c = min([c; GDP.T]);
    if numel(t) && c > t && t < GDP.T

        % return selected sequences
        %------------------------------------------------------------------
        O   = PDP.O(:,1:(t + 1));
        o   = PDP.o(:,1:(t + 1));

        break
    end
end


return

function Q = spm_get_O(RDP)
% hierarchical outcomes from a recursive MDP
% FORMAT Q = spm_get_O(RDP)
%
% Extracts hierarchal outcomes from a recursive MDP that have been
% accumulated during recursive inversion
%__________________________________________________________________________

% recover outcomes
%--------------------------------------------------------------------------
Q          = RDP.Q.O;
Q{end + 1} = RDP.O;

return


function spm_report(o,Nr,Nc)
% Plots sequence of moves in a MDP
% FORMAT spm_report(o,Nr,Nc)
% Nr - Number of rows (x dimension)
% Nc - Number of columns (x dimension)
%--------------------------------------------------------------------------
% If there are more than eight moves this subroutine will plot a movie

%% show sequence of moves
%==========================================================================
Nt    = 4;
T     = size(o,2);
for t = 1:T

    % movie
    %----------------------------------------------------------------------
    if T > Nt
        subplot(4,3,2), cla
    else
        subplot(4,4,t), cla
    end

    % plot
    %----------------------------------------------------------------------
    imagesc(reshape(o(:,t),Nr,Nc)), axis image, axis xy
    title(sprintf('Time %i',t),'FontSize',12)
    drawnow, pause(1/16)

end

return


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

function spm_show_outcomes(O,RG,Nr,Nc,Y)
% Plots the inferred sequence of moves
% FORMAT spm_show_outcomes(O,RG,Nr,Nc,[Y])
% O{n}   - Cell array of probabilistic, hierarchal outcomes
% RG{n}  - Cell array of hierarchical outcome groups
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
    subplot(2*Nn,1,n + Nn)
    imagesc(1 - spm_cat(O{n})), axis xy
    title(sprintf('Predictive posteriors at level %i',n),'FontSize',12)
end

% And repeat in image format using outcome groups
%--------------------------------------------------------------------------
for t = 1:T

    % Cycle over hierarchical levels (n)
    %----------------------------------------------------------------------
    for n = 1:Nn

        % time scaling for this level
        %------------------------------------------------------------------
        tn   = ceil(t*size(O{n},2)/T);
        if n == 1

            % observed outcomes
            %--------------------------------------------------------------
            X     = zeros(Nr*Nc,3);
            o     = O{n}(:,tn);
            for j = 1:numel(o)
                X(j,:) = spm_colour(o{j});
            end

            % predicted outcomes
            %--------------------------------------------------------------
            if nargin > 4
                P     = zeros(Nr*Nc,3);
                o     = Y{n}(:,tn);
                for j = 1:numel(o)
                    P(j,:) = spm_colour(o{j});
                end

                % plot first level outcomes as a coloured image
                %----------------------------------------------------------
                subplot(2*Nn,2,1)
                imagesc(reshape(X,Nr,Nc,3)), axis image, axis xy
                title(sprintf('Outcomes at time %i',t),'FontSize',12)

                % plot first level predictions as a coloured image
                %----------------------------------------------------------
                subplot(2*Nn,2,2)
                imagesc(reshape(P,Nr,Nc,3)), axis image, axis xy
                title(sprintf('Predictions at time %i',t),'FontSize',12)
                
            else

                % plot first level outcomes as a coloured image
                %----------------------------------------------------------
                subplot(2*Nn,1,1)
                imagesc(reshape(X,Nr,Nc,3)), axis image, axis xy
                title(sprintf('Outcomes at time %i',t),'FontSize',12)
            end

        else

            % higher level
            %--------------------------------------------------------------
            [nr,nc] = size(RG{n - 1});
            X       = cell(nr,nc);
            P       = cell(nr,nc);

            o     = O{n}(:,tn);
            m     = 1;
            for j = 1:numel(o)
                m = max(m,numel(o{j}));
            end
            m     = ceil(sqrt(m));

            for j = 1:numel(X)
                x    = zeros(m,m);
                p    = zeros(m,m);
                xo   = o{2*j - 1};
                po   = o{2*j - 0};
                x(1:numel(xo)) = xo;
                p(1:numel(po)) = po;
                X{j} = x;
                P{j} = p;
            end

            % And subsequent levels as (grouped) images of posteriors
            %--------------------------------------------------------------
            subplot(2*Nn,2,n*2 - 1)
            imagesc(1 - spm_cat(X)), axis image, axis xy
            title(sprintf('Latent states at time %i',t),'FontSize',12)

            subplot(2*Nn,2,n*2 - 0)
            imagesc(1 - spm_cat(P)), axis image, axis xy
            title(sprintf('Latent paths at time %i',t),'FontSize',12)

        end

    end

    % save movie
    %----------------------------------------------------------------------
    drawnow
    I(t) = getframe(gcf);

end

% Place movie in graphic subject
%--------------------------------------------------------------------------
set(gcf,'Userdata',[])
set(gcf,'Userdata',{I,8})
set(gcf,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

return


function [HID,CID,HITS,MISS] = spm_get_rewards(hid,cid,GDP,MDP)
% Gets rewarded (and restricted) states at the deepest level of an MDP
% FORMAT [HID,CID,HITS,MISS] = spm_get_rewards(hid,cid,GDP,MDP)
% hid  - List of intended (i.e., rewarded) goal states at the lowest level
% cid  - List of constrained (i.e., restricted) states at the lowest level
% GDP  - Generative process
% MDP  - Generative models (hierarchical)
%
% HID  - List of intended (i.e., rewarded) goal states at highest level
% CID  - List of constrained (i.e., restricted) states at highest level
% HITS - List of intended (i.e. rewarded) goal outcomes at lowest level
% MISS - List of constrained (i.e. restricted) outcomes at lowest level
%--------------------------------------------------------------------------
% This auxiliary routine identifies the episodes (i.e., paths) encoded by
% states at the deepest or highest level of a hierarchical MDP that entail
% one or more intended states at the lowest level.
%__________________________________________________________________________

% Find outcomes generated by rewarded and restricted states
%--------------------------------------------------------------------------
Ng    = numel(GDP.A);
Nh    = size(hid,2);
Nc    = size(cid,2);
HITS  = zeros(Ng,Nh);
MISS  = zeros(Ng,Nh);
for g = 1:Ng
    for h = 1:Nh
        ind       = num2cell(hid(:,h));
        HITS(g,h) = find(GDP.A{g}(:,ind{:}),1);
    end
    for h = 1:Nc
        ind       = num2cell(cid(:,h));
        MISS(g,h) = find(GDP.A{g}(:,ind{:}),1);
    end
end

% Convert hierarchical MDP to recursive RDP
%--------------------------------------------------------------------------
RDP       = spm_mdp2rdp(MDP);
[~,Ns,Nu] = spm_MDP_size(RDP);

% Solve for an episode encoded by the states of the highest level
%--------------------------------------------------------------------------
RDP.T = 1;
HID   = [];
CID   = [];
for s = 1:Ns(1)

    RDP.D{1} = sparse(s,1,1,Ns(1),1);
    RDP.E{1} = sparse(1,1,1,Nu(1),1);
    PDP      = spm_MDP_VB_XXX(RDP);

    % does this path elcit a reward?
    %----------------------------------------------------------------------
    o       = PDP.Q.o{1};
    
    % Uncomment to show paths
    %------------------------------------------------------------------
    % spm_report(o,6,9)
    if sum(ismember(o',HITS','rows'))
        HID = [HID,s];
    end
    if sum(ismember(o',MISS','rows'))
        CID = [CID,s];
    end

end

return



%% NOTES on alternative formulations
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

function O = spm_enrich(GDP,hid,cid)
% Enriches outcomes from a generative process (unused)
% FORMAT O = spm_enrich(GDP,hid,cid)
%
% This auxiliary routine generates outcomes from the random play and then
% selects sequences that intervene between repetitions of the initial
% state. A subset of sequences is then retained that has the greatest
% number of its (as indexed by hid in latent state space).
%__________________________________________________________________________

% repeated instances of first state
%--------------------------------------------------------------------------
i     = find(ismember(GDP.s',GDP.s(:,1)','rows'));

% for sequences between the most common output
%--------------------------------------------------------------------------
t     = {};
for k = 2:numel(i)

    %  indices of k-th sequence
    %----------------------------------------------------------------------
    t{k} = i(k - 1):(i(k) - 1);
    s    = GDP.s(:,t{k});

    % accept if there is a hit
    %----------------------------------------------------------------------
    h    = find(ismember(s',hid','rows'));
    H(k) = numel(h)/numel(t{k});
   
end

[i,j] = sort(H,'descend');
t     = t(j(1:round(numel(j)/4)));

% enriched outputs
%--------------------------------------------------------------------------
O = GDP.O(:,spm_cat(t));

return

% NOTES:
%==========================================================================
%     [Nf,Ns]  = spm_MDP_size(RDP);
%     H        = zeros(Ns,1);
%     H(HID)   =  4;
%     H(CID)   = -4;
%     RDP.H{1} = spm_softmax(H);
