function MDP = DEM_sharing
% Demo of active (visual) scene-construction
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
% inference (with belief propagation) to model active scene construction.
% It focuses on a discrete state space representation of a dynamic scene;
% generating visual snapshots at about the same frequency of saccadic eye
% movements. The generative model starts with latent states that correspond
% to natural kinds (e.g., objects) subject to natural laws (e.g., object
% invariance, classical mechanics, occlusion, and so on). A second latent
% factor (e.g., a 'where' stream) generates the fixation points in visual
% space for saccadic eye movements. The factors corresponding to multiple
% objects are themselves a Kronecker tensor product of attributes that
% depend upon each other; for example, position, velocity, pose, and
% non-spatial attributes that depend on spatial attributes. This
% interdependence means that object-specific attributes cannot be
% factorised; hence their combination as a tensor product (e.g., a 'what'
% stream).
%
% In what follows, we build a generative model, starting from state
% transitions that entail natural laws. Position refers to radial
% coordinates in egocentric space, implying a distinction between angular
% and radial (depth) states - and likewise for motion. This allows us to
% incorporate head orientation; in the form of head movements that
% reorientate the direction of gaze - that also depends upon the deployment
% of saccades in a head-centred frame of reference. Head movements are
% implemented, in the generative model, as moving objects in the egocentric
% frame of reference. This means that head movement is implemented via
% action-dependent transitions in location, while saccades are implemented
% via transitions among the latent states representing where gaze is
% deployed (in a head-centred frame of reference).
%
% Equipped with all of these hidden states, one can then complete a
% relatively simple generative model by specifying the likelihood mapping
% from hidden states to observations. This likelihood mapping is a high
% dimensional tensor - encoding all the high order dependencies generating
% visual input for the epoch in question. High order here refers to
% dependencies such as the interaction between two objects in the same line
% of sight that depends upon their relative depth to model occlusions.
%
% These outcomes are themselves discrete and multimodal. a high acuity
% modality models the parvocellular stream, with a restricted (central)
% field of view. This is complemented by two other modalities with a more
% peripheral field of view reporting contrast and motion energy, that is
% not spatially resolved (cf, the magnocellular stream). Note that in this
% construction (designed to generate the outputs of a computer vision
% scheme) motion is converted into a categorical (present versus absent)
% variable over discrete epochs of time. Note that the kind of scene
% construction and representation is implemented in egocentric and head
% centric frames of reference throughout. There is no part of the
% generative model that requires an allocentric representation - and yet,
% the agent can skilfully navigate a relatively complicated moving
% environment. in the example here, there are two inanimate objects (that
% play the role of landmarks) and an inanimate object (namely, a person who
% occasionally engages the agent with eye contact). This setup allows the
% simulation of reciprocal gaze and a primitive form of dyadic interaction.
% In other words, the prior preferences of this agent are to position
% itself and its direction of gaze to find someone who is looking at her.
%
% The code below is briefly annotated to illustrate how to build a
% generative model and then simulate active inference under that model, to
% produce relatively realistic sampling of a visual scene; namely, active
% scene construction. This inversion uses a sophisticated active inference
% scheme based upon a recursive estimation of expected free energy. This
% finesses the numerics because it uses belief propagation into the future
% - as opposed to marginal (variational) message passing. The numerical
% complexity of these models is a nontrivial issue: this is because most of
% the heavy lifting in the generative model is in the connectivity encoding
% dependencies that corresponds to high-dimensional tensors. In these
% simulations, the connectivity tensors are represented in working memory;
% whereas, in the brain or analogue (neuromorphic) implementations they
% would be simpler to instantiate.
%_________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


rng('default')

%% set up and preliminaries
%==========================================================================

% Latent states
%--------------------------------------------------------------------------
% First we specify the hidden states for each factor (i.e., object). The
% object has a number of attributes: line of sight, depth (foreground,
% background and distant), movement (none, right, left, retreat, approach),
% and disposition (e.g., looking in my direction or not).
%--------------------------------------------------------------------------

% number of: lines of sight, depth, motion and disposition
%--------------------------------------------------------------------------
los   = 9;                          % number of lines of sight
cm    = [3 5 7];                    % central lines of sight for each agent

% associate each object (i.e., natural kind) with a latent factor
%--------------------------------------------------------------------------
% This sets up a label structure with the names and cardinality of each
% factor and the accompanying latent states. These are further equipped
% with a number of actions that require a probability transition matrix for
% each allowable action. Here, the actions correspond to turning the head
% in one direction or another - or maintaining head position (the action
% 'stay')
%--------------------------------------------------------------------------
for i = 1:los
    Attribute{1}{i} = sprintf('%i',i);
end
Attribute{2} = {'close','near','gone'};
Attribute{3} = {'slowly','right','left','away','closer'};
Attribute{4} = {'Girl','Boy'};

% number of levels for each attribute
%--------------------------------------------------------------------------
for i = 1:numel(Attribute)
    N(i) = numel(Attribute{i});
end

for j = 1:prod(N)

    % Unpack object attributes
    %----------------------------------------------------------------------
    [a1,a2,a3,a4]    = spm_ind2sub(N,j);    % attributes of object
    label.name{1}{j} = [Attribute{4}{a4} ', ' ...
        Attribute{2}{a2} ' moving ' ...
        Attribute{3}{a3} ' at ' ...
        Attribute{1}{a1}];
end

% and a controllable factor (line of sight): saccadic eye movements
% and a syntax factor for language
%--------------------------------------------------------------------------
label.factor = [{'scene'} {'gaze'}];
label.name   = [label.name   {{'right','centre','left'}} ];
label.action = [{'none'} {{'right','centre','left'}}];

% with central and peripheral output modalities
%--------------------------------------------------------------------------
label.modality = {...
    'what',...
    'contrast-left',...
    'contrast-centre',...
    'contrast-right',...
    'motion-left',...
    'motion-centre',...
    'motion-right',...
    'spoken'};

% Central ('what') outcomes generated by object attributes
%--------------------------------------------------------------------------
label.outcome{1} = {  ...
    'person-near-right-yes',...
    'person-near-right-no',...
    'person-near-front-yes',...
    'person-near-front-no',...
    'person-near-left-yes',...
    'person-near-left-no',...
    'person-near-back',...
    'person-far-right',...
    'person-far-front',...
    'person-far-left',...
    'person-far-back',...
    'background'};

% and the label contrast and motion energy in peripheral field of vision
%--------------------------------------------------------------------------
for i = 2:7
    label.outcome{i} = {'near','far','none'};
end

% spoken outcome
%--------------------------------------------------------------------------
label.outcome{8} = label.name{1};


%% Transitions: B
%==========================================================================
% Next, we specify the probabilistic transitions of hidden states for each
% factor (i.e., object).
%--------------------------------------------------------------------------
for i = 1:size(N,1)
    for j = 1:size(N,2)
        I{i,j} = speye(N(i,j),N(i,j));
    end
end

% spatial and non-spatial attributes of objects
%==========================================================================
% Transitions among these states are characterised by movement that induces
% conditional dependencies between location and motion, the location
% factorises into lines of sight and depth (i.e., egocentric polar
% coordinates)
%--------------------------------------------------------------------------
disp('specifying generative model (c.f., training)'), disp(' ')
for f = 1:size(N,1)

    % for each head motion (disabled in this demo), shift the scene
    %----------------------------------------------------------------------
    for u = 1

        % transitions along lines of sight that depend up (angular) movement
        %------------------------------------------------------------------
        b     = cell(N(f,3),N(f,3));
        c     = 0;                            % head movements
        d     = [0 -1  1  0  0];              % object movements

        % shift angle appropriately by combining head and object movements
        %------------------------------------------------------------------
        for i = 1:N(f,3)
            b{i,i} = spm_speye(N(f,1),N(f,1),d(i) + c(u),1);
        end

        % and place in a Kronecker tensor product
        %------------------------------------------------------------------
        b     = spm_kron({spm_cat(b),I{f,2},I{f,4}});
        B{1}  = spm_permute_kron(b,N(f,[1,3,2,4]),[1,3,2,4]);

        % transitions between depth that depend upon (radial) movement
        %------------------------------------------------------------------
        b     = cell(N(f,3),N(f,3));
        d     = [0  0  0 -1  1];

        % shift depth appropriately (see spm_speye.m)
        %------------------------------------------------------------------
        for i = 1:N(f,3)
            b{i,i} = spm_speye(N(f,2),N(f,2),d(i),2);
        end
        b     = spm_kron({I{f,1}, spm_cat(b),I{f,4}});
        B{2}  = spm_permute_kron(b,N(f,[1 2 3 4]),[1 2 3 4]);

        % transitions between movement that depend upon depth
        %------------------------------------------------------------------
        b      = cell(3,3);
        b{1,1} = [...     % when near...
            1 0 0 1 1;    % when near and standing tend to move or withdraw
            4 1 0 0 0;    % when near continue moving or withdraw
            4 0 1 0 0;    % when near continue moving or withdraw
            1 1 1 0 0;    % when near and withdrawing stand still
            0 0 0 0 0];   % when near and approaching stand still
        b{2,2} = [...     % when far...
            1 0 0 1 1;    % when Far and standing tend to move or approach
            4 4 0 0 0;    % ...
            4 0 4 0 0;
            1 0 0 0 0;
            4 1 1 0 0];
        b{3,3} = [...     % when distant...
            1 1 1 1 1;
            0 4 0 0 0;
            0 0 4 0 0;
            0 0 0 0 0;
            4 0 0 0 0];

        for i = 1:3
            b{i,i} = b{i,i}(1:N(f,3),1:N(f,3));
        end
        b     = b(1:N(f,2),1:N(f,2));
        b     = spm_kron({I{f,1}, spm_cat(b),I{f,4}});
        b     = bsxfun(@rdivide,b,sum(b));
        B{3}  = spm_permute_kron(b,N(f,[1 3 2 4]),[1 3 2 4]);


        % transitions between disposition
        %------------------------------------------------------------------
        B{4}  = spm_kron({I{f,1},I{f,2},I{f,3},I{f,4}});

        % compose transitions over object attributes
        %------------------------------------------------------------------
        b  = 1;
        for i = 1:numel(B)
            b = full(b*B{i});
        end

        % transitions for this object and head motion
        %------------------------------------------------------------------
        T{f}(:,:,u) = b;

    end
end

%% Specify with controllable transitions among gaze directions
%--------------------------------------------------------------------------
nx    = numel(label.name{end});
nu    = numel(label.action{end});
b     = zeros(nx,nx,nu);
for u = 1:nu
    b(u,:,u) = 1;
end

% supplement B{:} and record the number of states for each factor
%--------------------------------------------------------------------------
B     = [T b];
for f = 1:numel(B)
    Nf(f) = size(B{f},1);
end


%% outcome probabilities: A
%==========================================================================
% Next, we specify the probabilistic mappings between latent states and
% outcomes with a tensor for each outcome modality, which models the high
% order interactions among the causes of outcomes (e.g., occlusion).
%--------------------------------------------------------------------------
for i = 1:numel(label.modality)
    mA{i} = zeros([numel(label.outcome{i}),Nf]);
end

% loop over every combination of object attributes to specify outcomes
%--------------------------------------------------------------------------
% Here, we will assume that objects in the distance cannot be seen and that
% objects in the foreground occlude background objects. By construction,
% the dispositional (non-spatial) attributes of an object can only be
% observed when an object is in the foreground.
%
% Attributes:
%  line of sight: 1,2,...N(1)
%  depth:         foreground, background and distant
%  movement:      none, right, left, withdraw, approach
%  disposition:   looking or not
%--------------------------------------------------------------------------

% for each agent
%--------------------------------------------------------------------------
nm    = numel(cm);                    % number of agents
for m = 1:nm
    A      = mA;                      % reset likelihood denser
    for o1 = 1:Nf(1)
        for u = 1:Nf(2)

            % agent's line of sight
            %==============================================================
            c        = cm(m) - 2 + u;

            % object attributes
            %==============================================================

            % Unpack object attributes
            %--------------------------------------------------------------
            [a1,a2,a3,a4] = spm_ind2sub(N,o1);   % attributes of object
            a        = [a1,a2,a3,a4];            % object x attribute array


            % {'what'}: outcomes
            %--------------------------------------------------------------
            %     'person-near-right-yes',...1
            %     'person-near-right-no',... 2
            %     'person-near-front-yes',...3
            %     'person-near-front-no',... 4
            %     'person-near-left-yes',... 5
            %     'person-near-left-no',...  6
            %     'person-near-back',...     7
            %     'person-far-right',...     8
            %     'person-far-front',...     9
            %     'person-far-left',...     10
            %     'person-far-back',...     11
            %     'distance',...            12
            %--------------------------------------------------------------


            % generate outcome from i-th object
            %--------------------------------------------------------------
            o            = spm_what(a,c);
            A{1}(o,o1,u) = 1;


            % {'contrast-left'}:
            %--------------------------------------------------------------
            % near...1
            % far... 2
            % none...3
            %--------------------------------------------------------------

            % and get contrast and motion energy
            %--------------------------------------------------------------
            o     = spm_contrast_energy(a,c - 1);
            A{2}(o,o1,u) = 1;
            o     = spm_motion_energy(a,c - 1);
            A{5}(o,o1,u) = 1;

            % and get contrast and motion energy
            %--------------------------------------------------------------
            o     = spm_contrast_energy(a,c);
            A{3}(o,o1,u) = 1;
            o     = spm_motion_energy(a,c);
            A{6}(o,o1,u) = 1;

            % and get contrast energy
            %--------------------------------------------------------------
            o     = spm_contrast_energy(a,c + 1);
            A{4}(o,o1,u) = 1;
            o     = spm_motion_energy(a,c + 1);
            A{7}(o,o1,u) = 1;

            % add linguistic likelihood model
            %--------------------------------------------------------------
            A{8}(o1,o1,u) = 1;

        end

    end

    % save this agent's likelihood tensor
    %----------------------------------------------------------------------
    Am{m} = A;

end


%% priors: (utility) C
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of log
% probabilities over outcomes. Here, the agent prefers eye contact but
% finds staring at someone with averted gaze aversive.
%--------------------------------------------------------------------------
% {'what'}: outcomes
%--------------------------------------------------------------------------
%     'person-near-right-yes',...1
%     'person-near-right-no',... 2
%     'person-near-front-yes',...3
%     'person-near-front-no',... 4
%     'person-near-left-yes',... 5
%     'person-near-left-no',...  6
%     'person-near-back',...     7
%     'person-far-right',...     8
%     'person-far-front',...     9
%     'person-far-left',...     10
%     'person-far-back',...     11
%     'background',...          12
%--------------------------------------------------------------------------
C{1}  = [4 4 4 4 4 4 4 2 2 2 2 0];

% and uninformative preferences over peripheral vision
%--------------------------------------------------------------------------
for i = 2:numel(A)
    C{i}  = zeros(1,size(A{i},1));
end

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify prior beliefs about
% initial states
%--------------------------------------------------------------------------
D{1} = ones(1,prod(N));                                    % uninformative
D{2} = [0,1,0];                                            % centre gaze

% allowable actions (with an action for each controllable state)
%--------------------------------------------------------------------------
U     = [ ...
    1 1;          % saccade to right line of sight
    1 2;          % saccade to centre line of sight
    1 3];         % saccade to left line of sight

% specify which states are shared and which are agent specific
%--------------------------------------------------------------------------
T      = 10;
m      = [1,0];                     % first (scene) is factor shared
n      = zeros(numel(A),1);
n(end) = -1;


% MDP Structure, specifying 64 epochs (i.e., 16 seconds of active vision)
%==========================================================================
mdp.T = T;                        % numer of moves
mdp.U = U;                        % actions
mdp.A = A;                        % likelihood probabilities
mdp.B = B;                        % transition probabilities
mdp.C = C;                        % prior preferences
mdp.D = D;                        % prior over initial states
mdp.N = 0;                        % policy depth
mdp.m = m;                        % shared outcomes
mdp.n = n;                        % shared states

mdp.label = label;


% create a cohort of agents, each with their unique likelihood mapping
%--------------------------------------------------------------------------
for i = 1:numel(Am)
    MDP(i,1)   = mdp;
    MDP(i,1).A = Am{i};
end

% Solve - an example with multiple epochs to illustrate how the agent
% resolves uncertainty about where she is looking and comes to track
% anybody who might be looking at her
%==========================================================================
disp('inverting generative model (c.f., active inference)'), disp(' ')

MDP  = spm_MDP_VB_XX(MDP);

% illustrate scene construction and perceptual synthesis
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference'); clf
spm_surveillance_percept(MDP,N,cm)

% now compare behaviour in the absence of language
%--------------------------------------------------------------------------
for i = 1:numel(Am)
    NDP(i,1)        = mdp;
    NDP(i,1).A      = Am{i};
end
for i = 1:numel(Am)
    NDP(i,1).s      = MDP(i,1).s;                  % reproduce scene
    NDP(i,1).s(2,:) = 0;                           % but not gaze
    NDP(i,1).A{end} = ones(size(MDP(i,1).A{end})); % disable language
end

NDP  = spm_MDP_VB_XX(NDP);

%% compare MDP and NDP
%--------------------------------------------------------------------------
spm_figure('GetWin','comparison'); clf

for m = 1:numel(MDP)

    % illustrate posterior beliefs and action selection
    %----------------------------------------------------------------------
    for t = 1:T
        q = MDP(m).X{1}(:,t);
        q = spm_sum(reshape(q,N),[2,3,4]);
        Q(:,t) = q;
    end
    subplot(6,2,(m*2) - 1), image(64*(1 - Q))
    xlabel('time'), ylabel('location')
    hold on
    i   = find(MDP(m).s(2,:) ~= NDP(m).s(2,:));
    plot(1:T,MDP(m).s(2,:) + cm(m) - 2,'.c','MarkerSize',16)
    plot(i,  MDP(m).s(2,i) + cm(m) - 2,'.r','MarkerSize',16)

    % repeat for naïve agent
    %----------------------------------------------------------------------
    for t = 1:T
        q = NDP(m).X{1}(:,t);
        q = spm_sum(reshape(q,N),[2,3,4]);
        Q(:,t) = q;
    end
    subplot(6,2,(m*2) - 0), image(64*(1 - Q))
    xlabel('time'), ylabel('location')
    hold on
    plot(1:T,NDP(m).s(2,:) + cm(m) - 2,'.c','MarkerSize',16)

    subplot(4,1,3)
    i   = find(MDP(m).s(2,:) == NDP(m).s(2,:));

    plot(i,MDP(m).F(i) - NDP(m).F(i)), hold on

end

subplot(4,1,3)
plot([0,T],[0,0],'k'), hold on
set(gca,'XLim',[0 T - 1]), box off
xlabel('time'), ylabel('natural units')
title('Free energy differences','FontSize',14)

% illustrate behavioural responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1));

% illustrate physiological responses (first unit at the fourth epoch)
%--------------------------------------------------------------------------
[i,j] = max(MDP(3).X{1}(:,3));
spm_figure('GetWin','Figure 2a with language'); clf
spm_MDP_VB_LFP(MDP(3),[j;3],1);

spm_figure('GetWin','Figure 2b without language'); clf
spm_MDP_VB_LFP(NDP(3),[j;3],1);

% focus on first two seconds
%--------------------------------------------------------------------------
for i = [2,3,4]
    spm_figure('GetWin','Figure 2a with language');
    subplot(3,2,i), set(gca,'XLim',[0 2])
    a = axis;
    spm_figure('GetWin','Figure 2b without language');
    subplot(3,2,i), set(gca,'XLim',[0 2])
    axis(a);
end


%% now we turn to learning the likelihood tensor modelling language
% acquisition. We will render the middle agent naïve to language and see if
% it can learn linguistic mappings from heard modalities to concepts
%==========================================================================
rng('default')
spm_figure('GetWin','belief updating'); clf
spm_figure('GetWin','language acquisition'); clf
spm_figure('GetWin','generations'); clf

clear MDP
ndp   = mdp;
NG    = 3*2;                               % number of generations
NT    = 128;                               % number of games
ndp.T = 16;                                % number of epochs per game
ndp.a = ndp.A;                             % initialise Dirichlet likelihood
ndp   = rmfield(ndp,'A');                  % remove likelihood array

BMR.g = 8;                    % - modality [default: 1]
BMR.f = 2;                    % - hidden factors to sum over [default: 0]
BMR.T = 1;                    % - Occam's window [default: 2]

% set up parents
%--------------------------------------------------------------------------
for m = 1:3                                % for each agent

    % likelihood mappings
    %----------------------------------------------------------------------
    MDP(m,1) = ndp;                        % initialise
    for i = 1:numel(Am{m})
        MDP(m,1).a{i} = Am{m}{i}*128 + 1;
    end
end

for g = 1:0%NG                               % for each generation

    % add a child of agent 1 + rem(g,3)
    %======================================================================
    parent   = 1 + rem(g,3);
    MDP(4,1) = MDP(parent,1);

    % initialise learnable language mapping for the child
    %----------------------------------------------------------------------
    MDP(4).a{8}  = MDP(4).a{8}*0 + 1/256;
    MDP(4).a0{8} = MDP(4).a{8};

    % experience-dependent learning over NT exposures
    %======================================================================
    KL    = [];
    FQ    = [];
    FA    = [];
    for t = 1:NT                           % for each exposure

        % active inference and learning
        %------------------------------------------------------------------
        PDP  = spm_MDP_VB_XX(MDP);

        % illustrate belief updating and learning
        %------------------------------------------------------------------
        if g == 1 && t < 4
            DDP(1,t) = PDP(4);
            if t == 3
                spm_figure('GetWin','belief updating');
                spm_MDP_VB_game(DDP);
            end
        end

        spm_figure('GetWin','language acquisition');
        %------------------------------------------------------------------
        for j = 1:numel(N)
            if t ~= 7 && g == 1
                subplot(12,8,(j - 1)*8 + min(t,8))
                a = spm_contract_a(PDP(4).a{8},N,j);
                a = spm_norm(a);
                imagesc(a), axis image
            end
        end, drawnow

        PDP(4) = spm_MDP_VB_sleep(PDP(4),BMR);

        % update Dirichlet parameters
        %------------------------------------------------------------------
        MDP(4).a{8}  = PDP(4).a{8};
        MDP(4).a0{8} = PDP(4).a0{8};

        % report acquisition in terms of divergence
        %------------------------------------------------------------------
        q     = MDP(4).a{8}(:,:);
        p     = MDP(parent).a{8}(:,:);


        % report acquisition in terms of correlations and free energies
        %----------------------------------------------------------------------
        subplot(3,3,7)
        q     = MDP(4).a{8}(:,:);
        p     = MDP(parent).a{8}(:,:);
        c     = corrcoef(q(:),p(:));
        KL(t) = c(1,2);
        plot(1:t,KL), set(gca,'XLim',[1,NT]);
        title({'Correlations','(likelihood)'}),xlabel('time'),ylabel('correlation')
        box off, axis square

        FQ(:,t) = sum(spm_cat({PDP.F}'),2);
        FA(:,t) = sum(spm_cat({PDP.Fa}'),2);
        FQ      = spm_conv(FQ,0,16);
        FA      = spm_conv(FA,0,16);

        subplot(3,3,8)
        plot(min(FQ(:)) - FQ'), set(gca,'XLim',[1,NT])
        title({'Free energy','(inference)'}),xlabel('time'),ylabel('natural units')
        box off, axis square

        subplot(3,3,9)
        plot(min(FA(:)) - FA'), set(gca,'XLim',[1,NT])
        title({'Free energy','(learning)'}),xlabel('time'),ylabel('natural units')
        box off, axis square

    end

    % illustrate language conservation over parents
    %----------------------------------------------------------------------
    spm_figure('GetWin','generations');
    for i = 1:3
        for j = 1:numel(N)
            subplot(12,8,(j - 1)*8 + i + 5*(g > 1))
            a = spm_contract_a(MDP(i).a{8},N,j);
            a = spm_norm(a);
            imagesc(a), axis image
        end
    end, drawnow

    % replace parent with a child and continue the next generation
    %----------------------------------------------------------------------
    MDP(parent,1) = MDP(4,1);

end



%% finally, see if a language emerges via free energy minimisation
%==========================================================================
rng('default')
spm_figure('GetWin','language emergence'); clf

ndp    = mdp;
NT     = 1024;                               % number of games
ndp.T  = 32;                                % number of epochs per game
ndp.a  = ndp.A;                             % initialise Dirichlet likelihood
ndp.a0 = ndp.A;                             % initialise Dirichlet likelihood
ndp    = rmfield(ndp,'A');                  % remove likelihood array

BMR.g  = 8;                                 % modality [default: 1]
BMR.f  = 2;                                 % hidden factors to sum over
BMR.T  = 2;
clear OPTIONS
OPTIONS.BMR = BMR;
OPTIONS.eta = 1;

% set up naïve population
%--------------------------------------------------------------------------
clear MDP
for m = 1:3                                 % for each agent

    % likelihood mappings
    %----------------------------------------------------------------------
    MDP(m,1) = ndp;                         % initialise
    for i = 1:numel(Am{m})
        MDP(m,1).a{i} = Am{m}{i}*128;
    end

    % initialise learnable language mapping for all agents
    %----------------------------------------------------------------------
    dA            = 1/32;
    MDP(m,1).a{8} = abs(randn(size(MDP(m).a{8}))*dA/4 + dA);

end


% illustrate language conservation over parents
%----------------------------------------------------------------------
spm_figure('GetWin','emergence'); clf
for i = 1:numel(MDP)
    for j = 1:numel(N)
        subplot(12,8,(j - 1)*8 + i)
        a = spm_contract_a(MDP(i).a{8},N,j);
        a = spm_norm(a);
        imagesc(a), axis image
    end
end, drawnow


%% experience-dependent learning over NT exposures
%==========================================================================
k     = 1:64;
KL    = [];
FQ    = [];
FA    = [];
EF    = [];
for t = 1:NT                                 % for each exposure

    spm_figure('GetWin','language emergence');
    %----------------------------------------------------------------------
    [a,j,n] = spm_dir_sort(MDP(1).a{8}(:,:,1));
    for i = 1:numel(MDP)

        % plot likelihood tensor (unfolded)
        %------------------------------------------------------------------
        subplot(6,4,(i - 1)*4 + min(t,4))
        a = spm_dir_norm(MDP(i).a{8}(:,:,1));
        a = a(j,n);
        imagesc(a(k,k)), axis image

        % plot correlations among Dirichlet counts
        %------------------------------------------------------------------
        subplot(6,4,numel(MDP)*4 + min(t,4))
        for m = (i + 1):numel(MDP)
            plot(MDP(i).a{8}(:),MDP(m).a{8}(:),'.'), axis square, hold on
        end

    end, drawnow
    set(gca,'ColorOrderIndex',1);

    % active inference and learning; i.e., exposure
    %----------------------------------------------------------------------
    PDP = spm_MDP_VB_XX(MDP);

    % update Dirichlet parameters
    %----------------------------------------------------------------------
    MDP = spm_MDP_VB_update(MDP,PDP,OPTIONS);

    % report acquisition in terms of correlations and free energies
    %----------------------------------------------------------------------
    subplot(3,4,9), hold off
    for i = 1:numel(MDP)
        for j = (i + 1):numel(MDP)
            q         = MDP(i).a{8}(:,:);
            p         = MDP(j).a{8}(:,:);
            c         = corrcoef(q(:),p(:));
            KL(t,i,j) = c(1,2);
            plot(1:t,KL(:,i,j)), set(gca,'XLim',[1,NT]); hold on
        end
        EF(t,i) = spm_MDP_MI(MDP(i).a{8}(:,:));
    end
    title({'Correlations','(likelihood)'}),xlabel('time'),ylabel('correlation')
    box off, axis square

    FQ(:,t) = sum(spm_cat({PDP.F}'),2);
    FA(:,t) = sum(spm_cat({PDP.Fa}'),2);
    FQ      = spm_conv(FQ,0,16);
    FA      = spm_conv(FA,0,16);

    subplot(3,4,10)
    plot(min(FQ(:)) - FQ'), set(gca,'XLim',[1,NT])
    title({'Free energy','(inference)'}),xlabel('time'),ylabel('natural units')
    box off, axis square

    subplot(3,4,11)
    plot(min(FA(:)) - FA'), set(gca,'XLim',[1,NT])
    title({'Free energy','(learning)'}),xlabel('time'),ylabel('natural units')
    box off, axis square

    subplot(3,4,12)
    plot(EF), set(gca,'XLim',[1,NT])
    title({'expected Free energy','(learning)'}),xlabel('time'),ylabel('natural units')
    box off, axis square

end

%% illustrate language conservation over parents
%--------------------------------------------------------------------------
spm_figure('GetWin','emergence');
for i = 1:numel(MDP)
    for j = 1:numel(N)
        subplot(12,8,(j - 1)*8 + i + 5)
        a = spm_contract_a(MDP(i).a{8},N,j);
        a = spm_norm(a);
        imagesc(a), axis image
    end
end, drawnow

% repeat but align with original syntax
%--------------------------------------------------------------------------
a     = 0;
for i = 1:numel(MDP)
    a = a + MDP(i).a{8};
end
a     = spm_dir_norm(a);
K     = a(:,:)*A{8}(:,:)';
[m,k] = max(K');
for i = 1:numel(MDP)
    for j = 1:numel(N)
        subplot(12,8,(j - 1)*8 + i + 3 + 40)
        a = MDP(i).a{8};
        k = 1:size(a,1);
        a = spm_contract_a(a(k,:,:,:,:),N,j);
        a = spm_norm(a);
        imagesc(a), axis image
    end
end, drawnow

return






% subroutines
%==========================================================================

function o = spm_what(a,c)
% returns an outcome vector from a list of attributes
% FORMAT o = spm_what(a,c)
% a   - attribute
% c   - agents line of sight
%
% {'what'}: outcomes
%--------------------------------------------------------------------------
%     'person-near-right-yes',...1
%     'person-near-right-no',... 2
%     'person-near-front-yes',...3
%     'person-near-front-no',... 4
%     'person-near-left-yes',... 5
%     'person-near-left-no',...  6
%     'person-near-back',...     7
%     'person-far-right',...     8
%     'person-far-front',...     9
%     'person-far-left',...     10
%     'person-far-back',...     11
%     'background',...          12
%--------------------------------------------------------------------------

% if there is an object or natural kind
%--------------------------------------------------------------------------
if a(1) == c

    % animate object:
    %      depth        direction    disposition
    %----------------------------------------------------------------------
    if     a(2) == 1 && a(3) == 2 && a(4) == 1
        o = 1;
    elseif a(2) == 1 && a(3) == 2 && a(4) == 2
        o = 2;
    elseif a(2) == 1 && (a(3) == 1 || a(3) == 5) && a(4) == 1
        o = 3;
    elseif a(2) == 1 && (a(3) == 1 || a(3) == 5) && a(4) == 2
        o = 4;
    elseif a(2) == 1 && a(3) == 3 && a(4) == 1
        o = 5;
    elseif a(2) == 1 && a(3) == 3 && a(4) == 2
        o = 6;
    elseif a(2) == 1 && a(3) == 4
        o = 7;
    elseif a(2) == 2 && a(3) == 2
        o = 8;
    elseif a(2) == 2 && (a(3) == 1 || a(3) == 5)
        o = 9;
    elseif a(2) == 2 && a(3) == 3
        o = 10;
    elseif a(2) == 2 && a(3) == 4
        o = 11;
    else

        %  background object
        %------------------------------------------------------------------
        o = 12;
    end

else

    % nothing to see
    %----------------------------------------------------------------------
    o = 12;

end


function o = spm_contrast_energy(a,c)
% returns an outcome vector from a list of attributes
% FORMAT o = spm_contrast_energy(a,c)
% a   - attribute
% c   - agents line of sight
%
% {'where'}: outcomes
%--------------------------------------------------------------------------
% near...1
% far... 2
% none...3
%--------------------------------------------------------------------------

% if there is an object or natural kind
%--------------------------------------------------------------------------
if a(1) == c
    if  a(2) == 1
        o = 1;                % if an object is in foreground
    elseif a(2) == 2
        o = 2;                % if an object in background
    else
        o = 3;                % object is in the distance
    end
else
    o = 3;                    % or none in this line of sight
end

return


function o = spm_motion_energy(a,c)
% returns an outcome vector from a list of attributes
% FORMAT o = spm_motion_energy(a,c)
% a   -  attribute
% c   - agents line of sight
%
% {'where'}: outcomes
%--------------------------------------------------------------------------
% near...1
% far... 2
% none...3
%--------------------------------------------------------------------------

% if there is an object or natural kind
%--------------------------------------------------------------------------
if a(1) == c
    if a(3) > 1                % if there is motion
        if a(2) == 1
            o = 1;             % in foreground
        elseif a(2) == 2
            o = 2;             % or background
        else
            o = 3;             % in the distance
        end
    else
        o = 3;                 % not moving
    end
else
    o = 3;                     % no object in line of sight
end

return


function spm_surveillance_percept(MDP,N,c)
%% illustrates visual search graphically
%--------------------------------------------------------------------------
% number of: lines of sight, depth, object motion and disposition
%--------------------------------------------------------------------------
% N     = [7,3,5,2];               % an animate object   (e.g., person)


% load images
%--------------------------------------------------------------------------
load DEM_scene
outcomes{1}(:,:,[12 13]) = [];
[k,l,m] = size(outcomes{1});

% {'what'}: outcomes{1}
%--------------------------------------------------------------------------
%     'person-near-right-yes',...1
%     'person-near-right-no',... 2
%     'person-near-front-yes',...3
%     'person-near-front-no',... 4
%     'person-near-left-yes',... 5
%     'person-near-left-no',...  6
%     'person-near-back',...     7
%     'person-far-right',...     8
%     'person-far-front',...     9
%     'person-far-left',...     10
%     'person-far-back',...     11
%     'landmark-near',...       12
%     'landmark-far',...        13
%     'background',...          14
%--------------------------------------------------------------------------
spoken = MDP(1).label.outcome{end};

% loop over agents
%--------------------------------------------------------------------------
Nm    = size(MDP,1);
for m = 1:Nm

    % loop over time
    %----------------------------------------------------------------------
    for t = 1:MDP(m).T

        % what the agent actually sees
        %==================================================================
        % peripheral vision, or magnocellular (contrast energy)
        % foveal (central) vision, or parvocellular (contrast energy)
        % peripheral vision, or magnocellular (motion energy)
        %------------------------------------------------------------------
        seen  = { ...
            outcomes{2}(:,:,MDP(m).o(2,t)) outcomes{2}(:,:,MDP(m).o(3,t)) outcomes{2}(:,:,MDP(m).o(4,t));
            outcomes{1}(:,:,end)           outcomes{1}(:,:,MDP(m).o(1,t)) outcomes{1}(:,:,end);
            outcomes{2}(:,:,MDP(m).o(5,t)) outcomes{2}(:,:,MDP(m).o(6,t)) outcomes{2}(:,:,MDP(m).o(7,t))};
        seen  = spm_cat(seen);

        subplot(4,Nm,m)
        image(spm_cat(seen))
        axis image, box off
        text(32,32,spoken{MDP(m).o(8,t)},'Color','r','FontSize',10)
        vision{m}(t) = getframe(gca);

        % actual scene
        %==================================================================
        subplot(4,1,2)

        % Unpack object attributes for this combination of objects
        %------------------------------------------------------------------
        [a1,a2,a3,a4] = spm_ind2sub(N,MDP(m).s(1,t)); % attributes of ith object
        a             = [a1,a2,a3,a4];                % attribute array

        % nearest (foreground or background) object in line of sight
        %------------------------------------------------------------------
        for j = 1:N(1)
            o    = spm_what(a,j);
            scene{j} = outcomes{1}(:,:,o);
        end

        % add direction of gaze and display
        %------------------------------------------------------------------
        s    = MDP(m).s(end,t);
        s    = c(m) - 2 + s;
        image(spm_cat(scene))
        text(l*(s - 1/2),k/2,'+','FontSize',32,'Color','r')
        axis image, box off
        actual(t) = getframe(gca);


        % what the agent thinks she sees
        %==================================================================
        subplot(2*Nm,1,Nm + m)

        % find the most likely state of each object
        %------------------------------------------------------------------
        [q,s] = max(MDP(m).X{1}(:,t));


        % Unpack object attributes for this combination of objects
        %------------------------------------------------------------------
        [a1,a2,a3,a4] = spm_ind2sub(N,s);         % attributes of ith object
        a             = [a1,a2,a3,a4];            % object attribute array

        % nearest (foreground or background) object in line of sight
        %------------------------------------------------------------------
        for j = 1:N(1)
            o    = spm_what(a,j);
            scene{j} = outcomes{1}(:,:,o);
        end

        % add direction of gaze and display
        %------------------------------------------------------------------
        [q,s] = max(MDP(m).X{end}(:,t));
        s    = c(m) - 2 + s;
        image(spm_cat(scene))
        text(l*(s - 1/2),k/2,'+','FontSize',32,'Color','r')
        axis image, box off
        percept{m}(t) = getframe(gca);

    end
end

% assign movies to each graph object
%--------------------------------------------------------------------------
for m = 1:Nm
    subplot(4,Nm,m)
    set(gca,'Userdata',{vision{m},4})
    set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
    title('Visual samples','FontSize',14)
end

subplot(4,1,2)
set(gca,'Userdata',{actual,4})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Actual scene','FontSize',14)

for m = 1:Nm
    subplot(2*Nm,1,Nm + m)
    set(gca,'Userdata',{percept{m},4})
    set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
    title('Perceived scene','FontSize',14)
end

return

function a  = spm_contract_a(A,N,i)
% sums the elements of a likelihood mapping for attribute i
%--------------------------------------------------------------------------
a = zeros(N(i));                                  % contracted tensor
A = bsxfun(@rdivide,A,sum(A,1));                  % normalise Dirichlet counts
f = 2:6; f(i) = [];
A = squeeze(spm_sum(reshape(A,[prod(N) N 3]),f)); % marginalise (states)
f = 1:4; f(i) = [];
for i = 1:size(A,2)
    a(:,i) = spm_sum(reshape(A(:,i),N),f);        % marginalise (outcomes)
end

return



function A  = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A           = bsxfun(@rdivide,A,sum(A,1));
A(isnan(A)) = 1/size(A,1);
