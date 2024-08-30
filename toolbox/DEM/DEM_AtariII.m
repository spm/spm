function MDP = DEM_AtariII
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
Sc = 3;                                % Spatial scaling

G  = @(Nr,Nc) spm_MDP_breakout(Nr,Nc); % game
NT = 2048;                             % exposures (training)
G  = @(Nr,Nc) spm_MDP_pong(Nr,Nc);     % game
NT = 1024;                             % exposures (training)

NT = 256;                             % exposures (training)



%% Simulate learning and subsequent performance
%--------------------------------------------------------------------------
[GDP,hid,cid,con,RGB] = G(Nr,Nc);
GDP.tau = 2;                          % smoothness of random paths
GDP.T   = (Nr + 2)*4;                 % upper bound on length of paths


% Determine length of episodes Ne
%--------------------------------------------------------------------------
PDP   = spm_MDP_generate(GDP);
O     = PDP.O;
MDP   = spm_faster_structure_learning(O,[Nr,Nc],Sc);
Nm    = numel(MDP);
Ne    = 2^(Nm - 1);

% generate (probabilistic) outcomes under random actions
%==========================================================================
spm_figure('GetWin','Gameplay'); clf

dcat  = @(O) spm_cat([O(:,1:end - 1); O(:,2:end)]);

O     = {};
R     = {};
s     = {};
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
        else
            O = PDP.O(:,1:t);
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

    else

        % unused (random) paths
        %------------------------------------------------------------------
        R{end + 1} = PDP.O;
        s{end + 1} = PDP.s;

    end

    % report
    %----------------------------------------------------------------------
    Nt   = size(O,2);
    clc; fprintf('Number of samples %i (%i)\n',Nt,(n*GDP.T))

    % break if a sufficient number of episodes have been accumulated
    %----------------------------------------------------------------------
    if Nt > NT
        O = O(:,1:NT);
        break
    end

    % break if recurrence
    %----------------------------------------------------------------------
    if Nt > Ne
        V       = dcat(O(:,(Nt - Ne + 1):Nt));
        U       = dcat(O(:,1:(Nt - Ne)));
        e       = find(ismember(V',U','rows'));
        if numel(e)
            i     = find(ismember(U',V(:,e)','rows'),1,'first');
            inset = 1:(i - 1);
            orbit = i:(Nt - Ne + (e - 1));

            % loop learning data
            %------------------------------------------------------------------
            i     = inset;
            for e = 1:log2(Ne)
                i = [i,orbit];
            end
            O = O(:,i);
            break
        end
    end

end

% RG structure learning
%==========================================================================
MDP = spm_faster_structure_learning(O,[Nr,Nc],3);
Ns  = size(MDP{end}.B,1);

% more structure learning
%==========================================================================
for t = 2:Ne
    MDP = spm_merge_structure_learning(O(:,t:end),MDP);
end

% rewarded events
%--------------------------------------------------------------------------
[HID,CID,HITS,~,S] = spm_get_episodes(hid,cid,GDP,MDP);


% Illustrate orbits
%==========================================================================
spm_figure('GetWin','Flows'); clf

subplot(2,2,1)
spm_dir_orbits(MDP{end}.B{1},HID,[],Ns);
title('Flows: orbit')

Q     = {};
for q = 2:3

    % start of attracting paths
    %----------------------------------------------------------------------
    U     = [];
    for m = 1:numel(S)
        U(:,m) = dcat(S{m}(:,1:2));
    end

    % for each realization
    %----------------------------------------------------------------------
    for n = 1:numel(R)
        
        % for each path in attracting set
        %------------------------------------------------------------------
        V     = dcat(R{n});
        for m = find(ismember(U',V','rows')')
            t     = find(ismember(V',U(:,m)','rows'));
            t     = t(t > Ne);
            for i = 1:numel(t)

                % pre-path
                %----------------------------------------------------------
                j   = (1:Ne) - Ne + t(i) - 1 ;
                if ~any(ismember(s{n}(:,j)',cid','rows'))

                    % select new path
                    %------------------------------------------------------
                    Q{end + 1} = R{n}(:,j);

                    % and append
                    %------------------------------------------------------
                    MDP = spm_merge_structure_learning([Q{end},S{m}],MDP);

                end
            end
        end
    end

    % new attracting paths
    %----------------------------------------------------------------------
    S = Q;
    Q = {};

    % Illustrate orbits
    %----------------------------------------------------------------------
    subplot(2,2,q)
    spm_dir_orbits(MDP{end}.B{1},HID);
    title(sprintf('Flows: %i-order',q))
    drawnow

end

% Illustrate orbits
%==========================================================================
spm_figure('GetWin','Orbits'); clf

subplot(2,2,1)
spm_dir_orbits(MDP{end}.B{1},HID);
title('Orbits & goals')

% paths to hits
%--------------------------------------------------------------------------
subplot(2,1,2)
B     = sum(MDP{Nm}.B{1},3) > 0;
Ns    = size(B,1);
h     = sparse(1,HID,1,1,Ns);
I     = [];
for t = 1:64
    I(t,:) = h;
    h      = (h + h*B) > 0;
end
imagesc(I), hold on 
plot(HID,0*HID + 1/2,'or'), hold off
title('Paths to hits','FontSize',14)
xlabel('latent states'), ylabel('time steps'), axis square

% priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Priors'); clf
spm_RDP_params(MDP)


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



% Active inference
%==========================================================================
% In what follows, we engage inductive planning as inference, with explicit
% action; namely, the generative process is used to engage actions that
% reproduce predicted outcomes (here, purely visual) at the lowest level.
%--------------------------------------------------------------------------

%% create hierarchical model with prior concentration parameters
%--------------------------------------------------------------------------
MDP{1}.GA  = GDP.A;
MDP{1}.GB  = GDP.B;
MDP{1}.GD  = GDP.D;
MDP{1}.GE  = GDP.E;
MDP{1}.GU  = GDP.U;

MDP{1}.s   = [3;5]; %%%%



for g = 1:numel(GDP.A)
    MDP{1}.ID.A{g} = 1:(ndims(GDP.A{g}) - 1);
end

MDP{1}.ID.control = con;                    % controlled outcomes
MDP{1}.chi = 512;                             %%%% sticky action/shaky hand


% enable active learning (with minimal forgetting)
%--------------------------------------------------------------------------
for m = 1:numel(MDP)
    MDP{m}.beta = 4;
    MDP{m}.eta  = 32;
end

% train: with small concentration parameters
%--------------------------------------------------------------------------
FIX.A      = 1;                             % learn likelihood
FIX.B      = 1;                             % but not transitions
RDP        = spm_mdp2rdp(MDP,0,1/1024,2,FIX);    %%%%%%

RDP.U      = 1;
RDP.T      = 128;
RDP.id.hid = HID;
RDP.id.cid = CID;

PDP        = spm_MDP_VB_XXX(RDP);

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference'); clf
spm_show_RGB(PDP,RGB,8,false);

% add ELBOs
%--------------------------------------------------------------------------
subplot(Nm + 3,2,2*(Nm + 1))
T     = numel(PDP.Q.E{1});
t     = linspace(1,T,RDP.T);
plot(t,PDP.F), hold on
for n = 1:numel(PDP.Q.E)
    t = linspace(1,T,numel(PDP.Q.E{n}));
    plot(t,PDP.Q.E{n})
end

% and hits
%--------------------------------------------------------------------------
h     = find(ismember(PDP.Q.o{1}',HITS','rows'));
plot(h,ones(size(h)),'.r','MarkerSize',16)
title('ELBO'), spm_axis tight



% Illustrate in latent state space 
%-=========================================================================
spm_figure('GetWin','Orbits');

subplot(2,2,2)
spm_dir_orbits(PDP.B{1},HID,PDP.X{1});

return
