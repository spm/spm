function MDP = DEM_sharingX
% Demo of active (visual) scene-construction
%__________________________________________________________________________
%
% This demonstration illustrates the emergence of belief sharing and
% distributed (i.e., federated) inference when posterior beliefs are
% broadcast among multiple agents. It can be read from a number of
% perspectives: from a systemic perspective, this emergence can be read as
% the minimisation of joint free energy that ensues when inference,
% learning and selection are simulated as a free energy minimising
% processes. A more anthropomorphic interpretation is in terms of language
% and communication; namely, the acquisition of language, and its
% transmission over generations, in the spirit of cultural niche
% construction.
% 
% The particular example here uses three agents that can be thought of as
% three sentinels maintaining surveillance for potential predators.
% Crucially, each of the three agents has a restricted field of view, such
% that each can only see about 1/3 of the horizon. In addition to visual
% and proprioceptive (gaze-related) observation modalities, each agent can
% also hear the others (but not itself). The agents have a simple
% (identity) likelihood mapping between posterior beliefs — about latent
% states that are shared (or conserved) over agents — and output modalities
% or communication channels.
% 
% Here, the birds report the location of a potential predator in terms of
% its radial location in an allocentric frame of reference and its
% proximity (near or far). In addition, the agents can report whether the
% predator (i.e., subject of surveillance) is friendly or not. This active
% inference is deliberately made difficult (in addition to the limited
% field of view afforded each agent): first, the agents cannot see motion.
% This means movement of the subject has to be inferred from successive
% observations of location. This inference is crucial because it
% underwrites predictions of location at the next time step, which are
% broadcast to other agents. In the absence of this prediction, there would
% be less predictive value for the other agents. Second, the disambiguation
% between friend and foe is only possible when the predator is sufficiently
% close to each bird.
% 
% This scheme uses a generalised MDP, in which the paths are assumed to be
% unchanging during each epoch (unless they are controllable). This means
% for each episode of the simulation the predator is moving continuously in
% one direction or another (or remains stationary).
% 
% In what follows, we first demonstrate the utility of belief sharing, in
% resolving uncertainty about latent states, when their consequences can
% only be seen by one agent at a time. We then illustrate the acquisition
% of language (i.e., likelihood mappings) using active learning or
% accumulation of Dirichlet counts, where a child of one parent learns from
% its parents and conspecifics. Finally, we illustrate the emergence of
% language (i.e., precise and shared likelihood mappings) in an ensemble of
% agents with initially imprecise, random likelihood tensors, when they are
% exposed to the same environment. This emergence rests upon structure
% learning in which priors over the Dirichlet counts of likelihood tensors
% are updated using Bayesian model reduction – and the updates are in the
% direction of maximising expected free energy or the mutual information of
% the likelihood mapping. In short, distributed inference emerges from
% minimising (variational and expected) free energy at three levels of
% joint inference, learning and selection.
%_________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%% set up and preliminaries
%==========================================================================
rng(1);

% Latent states
%==========================================================================
% First, we will specify latent causes in terms of a mean field
% approximation; namely, factors and their associated states (called
% names). These can be specified as a structure of text fields as follows.
% Each factor has an associated dynamic, that may or may not be
% controllable called, action. In this example, there are four factors:
% location, proximity, pose and gaze, where location and gaze are equipped
% with different actions or paths. The agent controls the action on gaze,
% while the action associated with the location of an object represents
% movement of an object in the environment.
%--------------------------------------------------------------------------

% lines of sight for each agent
%--------------------------------------------------------------------------
los   = 9;                          % number of lines of sight
cm    = [3 5 7];                    % central lines of sight for each agent

% Latent factors, states and paths or actions
%--------------------------------------------------------------------------
label.factor = {'location','proximity','pose','gaze'};
for i = 1:los
    label.name{1}{i}   = sprintf('at %i',i);
end
label.name{2}   = {'close ','far '};
label.name{3}   = {'friend ','foe '};
label.name{4}   = {'right','centre','left'};

label.action{1} = {'stay','right','left'};
label.action{2} = {'none'};
label.action{3} = {'none'};
label.action{4} = {'right','centre','left'};

% number of levels for each state and path
%--------------------------------------------------------------------------
Nf = numel(label.name);
for f = 1:Nf
    Ns(f) = numel(label.name{f});
    Nu(f) = numel(label.action{f});
end


% Outcomes
%==========================================================================
% We now repeat the exercise for outputs, where output factors are called
% modalities and the accompanying states outcomes. Here, the first modality
% is a central (e.g., foveal) vision. Subsequent modalities are peripheral
% visions registering contrast; proprioceptive input reporting the
% direction of gaze and, crucially, three auditory modalities corresponding
% to the posterior beliefs about latent factors articulated by other
% agents.
%--------------------------------------------------------------------------
label.modality = {...
    'what',...
    'contrast-left',...
    'contrast-centre',...
    'contrast-right',...
    'gaze',...
    'location',...
    'proximity',...
    'pose'};

% Central ('what') outcomes generated by object attributes
%--------------------------------------------------------------------------
label.outcome{1} = {  ...
    'close friend',...
    'close foe',...
    'far person',...
    'nothing'};

% and label contrast energy in peripheral field of vision
%--------------------------------------------------------------------------
for g = 2:4
    label.outcome{g} = {'close','near','none'};
end

% proprioception
%--------------------------------------------------------------------------
label.outcome{5} = label.name{4};

% spoken outcome
%--------------------------------------------------------------------------
label.outcome{6} = label.name{1};
label.outcome{7} = label.name{2};
label.outcome{8} = label.name{3};

% number of levels for each attribute
%--------------------------------------------------------------------------
Ng = numel(label.modality);
for g = 1:Ng
    No(g) = numel(label.outcome{g});
end
Cg = 6:8;


%% Transitions: B
%==========================================================================
% Next, we specify the probabilistic transitions among hidden states for
% each Factor. These transitions correspond to a path or action, resulting
% in a three tensor for each factor. These generally have a relatively
% simple form; either a single identity matrix for states that remain in
% their initial states over and Epoque; or that move in a structured way.
% Here, location can change according to whether the subject of
% surveillance is moving to the right, not moving or moving to the left.
% She can also approach or withdraw, transitioning from near to far with a
% small probability. The final gaze factor is controllable and simply
% causes a transition to the right, the centre or the left.
%--------------------------------------------------------------------------
d     = [0 -1  1];
f     = 1;
for h = 1:Nu(f)
    B{f}(:,:,h) = full(spm_speye(Ns(f),Ns(f),d(h),1));
end
f     = 2;
for h = 1:Nu(f)
    B{f}(:,:,h) = full(spm_speye(Ns(f),Ns(f))) + 1/8;
end
f     = 3;
for h = 1:Nu(f)
    B{f}(:,:,h) = full(spm_speye(Ns(f),Ns(f)));
end
f     = 4;
for h = 1:Nu(f)
    B{f}(:,:,h) = zeros(Ns(f),Ns(f));
    B{f}(h,:,h) = 1;
end


%% outcome probabilities: A
%==========================================================================
% Next, we specify the probabilistic mappings between latent states and
% outcomes with a tensor for each outcome modality, which models the high
% order interactions among the causes of outcomes.
%--------------------------------------------------------------------------
% Attributes:
%  line of sight: 1,2,...Ns(1)
%  proximity:     Close and far
%  pose:          friend or foe
%--------------------------------------------------------------------------
for g = 1:5
    A0{g}     = zeros([No(g),Ns]);        % agent-specific outcomes
end
for f = 1:3
    A0{f + 5} = zeros([Ns(f),Ns(f)]);     % shared (communication) outcomes
end

% for each agent
%--------------------------------------------------------------------------
nm    = numel(cm);                        % number of agents
for m = 1:nm
    A = A0;                               % reset likelihood tensor
    for s1 = 1:Ns(1)                      % location
        for s2 = 1:Ns(2)                  % proximity
            for s3 = 1:Ns(3)              % pose
                for u = 1:Ns(4)           % gaze

                    % agent's line of sight
                    %======================================================
                    c  = cm(m) + u - 2;
                    a  = [s1,s2,s3,c];

                    % object attributes
                    %======================================================

                    % {'what'}: outcomes
                    %------------------------------------------------------
                    o  = spm_what(a,c);
                    A{1}(o,s1,s2,s3,u) = 1;

                    % contrast energy - right
                    %------------------------------------------------------
                    o  = spm_contrast_energy(a,c - 1);
                    A{2}(o,s1,s2,s3,u) = 1;

                    % contrast energy - center                                  
                    %------------------------------------------------------
                    o  = spm_contrast_energy(a,c);
                    A{3}(o,s1,s2,s3,u) = 1;

                    % contrast energy - left 
                    %------------------------------------------------------
                    o  = spm_contrast_energy(a,c + 1);
                    A{4}(o,s1,s2,s3,u) = 1;

                    % proprioception
                    %------------------------------------------------------
                    A{5}(u,s1,s2,s3,u) = 1;

                    % communication (see belief sharing)
                    %------------------------------------------------------
                    % A{6}(s1,s1,s2,s3,u) = 1;

                    % communication
                    %------------------------------------------------------
                    % A{7}(s2,s1,s2,s3,u) = 1;

                    % communication
                    %------------------------------------------------------
                    % A{8}(s3,s1,s2,s3,u) = 1;

                end
            end
        end
    end

    % belief sharing and domains of likelihood mappings (ID)
    %----------------------------------------------------------------------
    for g = Cg
        f       = g - 5;
        A{g}    = eye(Ns(f),Ns(f));
        id.A{g} = f;
    end
    ID{m} = id;                  

    % save this agent's likelihood tensor
    %----------------------------------------------------------------------
    Am{m} = A;

end


%% priors: (utility) C
%==========================================================================
% Finally, we have to specify the prior constraints in terms of log
% probabilities over outcomes. Here, there are no prior constraints which
% means that all active inference rests upon expected information gain
% (i.e., epistemic affordances).
%--------------------------------------------------------------------------
% {'what'}: outcomes
%--------------------------------------------------------------------------
%     'close friend',...
%     'close foe',...
%     'far person',...
%     'nothing'};
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify allowable actions (with
% an action for each controllable state; here, the fourth factor, namely,
% gaze)
%--------------------------------------------------------------------------
U      = zeros(1,Nf);
U(end) = 1;

% Now specify which states are shared and which are agent-specific
%--------------------------------------------------------------------------
m       = [1,1,1,0];               % scene factors are shared
n       = zeros(Ng,1);             % joint outcomes (Cg)
n(Cg,:) = -1;                      % are shared (with n < 0)

% MDP Structure, specifying 8 epochs (i.e., 4 seconds of active vision)
%==========================================================================
mdp.T = 8;                        % numer of moves
mdp.U = U;                        % controllable actions
mdp.A = A;                        % likelihood probabilities
mdp.B = B;                        % transition probabilities
mdp.C = C;                        % prior constraints
mdp.N = 0;                        % policy depth
mdp.m = m;                        % shared outcomes
mdp.n = n;                        % shared states

mdp.label = label;                % labels
mdp.id    = ID{1};                % domains
mdp.eta   = 32;                   % Dirichlet hyperprior [learnability]

% create a cohort of agents, each with their unique likelihood mapping
%--------------------------------------------------------------------------
mdp.s = 1;                        % start in the first location
for m = 1:numel(Am)
    MDP(m,1)    = mdp;
    MDP(m,1).A  = Am{m};
end

% Solve - an example with multiple epochs
%==========================================================================
OPTIONS.N = 1;                    % enable synthetic neural responses
MDP       = spm_MDP_VB_XXX(MDP,OPTIONS);

% illustrate scene construction and perceptual synthesis
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference'); clf
spm_surveillance_percept(MDP,cm)

% now compare behaviour in the absence of language
%--------------------------------------------------------------------------
for m = 1:numel(Am)
    NDP(m,1)   = mdp;
    NDP(m,1).A = Am{m};
end

% reproduce scene and disable language (with an imprecise likelihood)
%----------------------------------------------------------------------
for m = 1:numel(Am)
    NDP(m,1).s = MDP(m,1).s;
    for g = Cg
        NDP(m,1).A{g} = ones(size(MDP(m,1).A{g}));
    end
end

% Solve - an example with no communication
%==========================================================================
NDP  = spm_MDP_VB_XXX(NDP,OPTIONS);


%% compare MDP and NDP
%==========================================================================
spm_figure('GetWin','comparison'); clf
T     =  mdp.T;
for m = 1:numel(MDP)

    % illustrate posterior beliefs and action selection
    %----------------------------------------------------------------------
    subplot(6,2,(m*2) - 1), image(64*(1 - MDP(m).X{1}))
    xlabel('time'), ylabel('location')
    hold on

    % show where the agents are looking
    %----------------------------------------------------------------------
    plot(1:T, MDP(m).s(4,:) + cm(m) - 2,'.c','MarkerSize',16)
    plot(1:T, MDP(m).s(1,:),'.r','MarkerSize',16)

    % repeat for naïve agents
    %----------------------------------------------------------------------
    subplot(6,2,(m*2) - 0), image(64*(1 - NDP(m).X{1}))
    xlabel('time'), ylabel('location')
    hold on

    % show where the agents are looking
    %----------------------------------------------------------------------
    plot(1:T, NDP(m).s(4,:) + cm(m) - 2,'.c','MarkerSize',16)
    plot(1:T, NDP(m).s(1,:),'.r','MarkerSize',16)

    % compare in terms of free energy (negative ELBO)
    %----------------------------------------------------------------------
    subplot(4,1,3)
    plot(MDP(m).F - NDP(m).F), hold on

end

subplot(6,2,1)
title('With communication','FontSize',14)
subplot(6,2,2)
title('and without','FontSize',14)

subplot(4,1,3), 
hold on, plot([1 T],[0 0],'--'), hold off, box off
xlabel('time'), ylabel('natural units')
title('Free energy differences','FontSize',14)
legend({'first','second','third'})

% illustrate belief updating of second agent
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(2),1:4,1:5);

% illustrate physiological responses with and withoout language
%--------------------------------------------------------------------------
t     = 3;                               % at time t
m     = 2;                               % in agent m
[i,j] = max(MDP(m).X{1}(:,t) - NDP(m).X{1}(:,t));
spm_figure('GetWin','Figure 2a with language'); clf
spm_MDP_VB_LFP(MDP(m),[j;t],1);
subplot(3,2,6), d = axis;

spm_figure('GetWin','Figure 2b without language'); clf
spm_MDP_VB_LFP(NDP(m),[j;t],1);

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
subplot(3,2,6), axis(d);


%% Language acquisition and cultural niche construction
%==========================================================================
% Now, we turn to learning the likelihood tensor modelling language
% acquisition. Here, we will pair each agent in turn with a child who sees
% and hears the same things as their parent and accumulates Dirichlet
% parameters so that she can understand what her conspecifics believe.
% After a sufficient number of episodes (i.e., exposures) the child
% replaces her parent, and the next parent is selected. By repeating this
% four times, all the parents are replaced by children. The idea here is to
% see whether the language survives in a transgenerational sense.
%==========================================================================
clear MDP
spm_figure('GetWin','language acquisition'); clf
spm_figure('GetWin','generations'); clf

ndp   = rmfield(mdp,'s');
NG    = 4;                                 % number of generations
NT    = 32;                                % number of episodes
ndp.T = 16;                                % number of epochs per episodes
ndp.a = ndp.A;                             % initialise Dirichlet likelihood
ndp   = rmfield(ndp,'A');                  % remove likelihood array

BMR.g     = Cg;                            % modalities [default: 1]
BMR.T     = 0;                             % Occam's window [default: 0]
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.B = 0;                             % suppress replay
OPTIONS.BMR = BMR;                         % enable BMR

% set up likelihood mappings of parents
%--------------------------------------------------------------------------
for m = 1:3                                % for each agent
    MDP(m,1) = ndp;                        % initialise
    for i = 1:numel(Am{m})
        MDP(m,1).a{i} = Am{m}{i}*exp(16);
    end
end

for gn = 1:NG                              % for each generation

    % add a child of agent 1 + rem(g,3)
    %======================================================================
    parent   = 1 + rem(gn,3);
    MDP(4,1) = MDP(parent,1);

    % initialise learnable language mapping for the child
    %----------------------------------------------------------------------
    for g = Cg
        MDP(4).a{g}  = MDP(4).a{g}*0 + 1;
    end

    % experience-dependent learning over NT exposures
    %======================================================================
    KL    = [];                            % divergence (child vs parent)
    MI    = [];                            % Mutual information of likelihood
    FQ    = [];                            % Free energy (inference)
    FA    = [];                            % free energy (learning)
    for t = 1:NT                           % for each episode


        % active inference and learning
        %------------------------------------------------------------------
        PDP = spm_MDP_VB_XXX(MDP,OPTIONS);

        % update Dirichlet parameters with Bayesian model reduction
        %------------------------------------------------------------------
        MDP = spm_MDP_VB_update(MDP,PDP,OPTIONS);


        % report acquisition 
        %==================================================================

        % inference
        %------------------------------------------------------------------
        if false
            spm_figure('GetWin','inference');
            spm_MDP_VB_trial(PDP(4),1:4,1:8);

            subplot(6,2,5)
            bar(PDP(4).F)
            title('Free energies'),ylabel('agent')

            subplot(6,2,7)
            image(64 + full(spm_cat({PDP.F}')))
            title('Free energies'),ylabel('agent')
        end

        spm_figure('GetWin','language acquisition');
        %------------------------------------------------------------------
        for j = 1:3
            if t ~= 7 && gn == 1
                subplot(12,8,(j - 1)*8 + min(t,8))
                a = PDP(4).a{j + 5};
                a = spm_norm(a);
                image(a*64), axis image
                if t == 1, ylabel(label.factor{j}), end
            end
        end, drawnow
        
        % record divergence and free energies
        %------------------------------------------------------------------
        KL(t)   = spm_KL_cat(MDP(parent).a{6},MDP(4).a{6});
        MI(t)   = spm_MDP_MI(MDP(4).a{6});
        FQ(:,t) = sum(spm_cat({PDP.F}'),2);
        FA(:,t) = sum(spm_cat({PDP.Fa}'),2);

    end

    % illustrate language acquisition (first generation)
    %======================================================================
    if gn == 1
        subplot(4,2,5), plot(1:t,KL), set(gca,'XLim',[1,NT]);
        title({'KL divergence','(likelihood)'}),xlabel('time'),ylabel('nats')
        box off, axis square

        subplot(4,2,6), plot(-FA'), set(gca,'XLim',[1,NT])
        title({'Free energy','(learning)'}),xlabel('time'),ylabel('nats')
        box off, axis square

        subplot(4,2,7), plot(-FQ'), set(gca,'XLim',[1,NT])
        title({'Free energy','(inference)'}),xlabel('time'),ylabel('nats')
        box off, axis square

        subplot(4,2,8), plot(1:t,MI), set(gca,'XLim',[1,NT]);
        title({'MI [negative EFE]','(likelihood)'}),xlabel('time'),ylabel('nats')
        box off, axis square
    end

    % illustrate language conservation over parents
    %----------------------------------------------------------------------
    spm_figure('GetWin','generations');
    for i = 1:3
        for j = 1:3
            subplot(12,8,(j - 1)*8 + i + 5*(gn > 1))
            a = PDP(i).a{j + 5};
            a = spm_norm(a);
            image(a*64), axis image

            if i == 1
                ylabel(label.factor{j})
            end
            if j == 1 && gn == 1
                title(sprintf('agent %i',i))
            end
            if j == 1 && i == 1 && g > 1
                title(sprintf('generation %i',gn))
            end
        end
    end, drawnow

    % replace parent with a child and continue with the next generation
    %======================================================================
    MDP(parent,1) = MDP(4,1);

end


%% finally, see if a language emerges via free energy minimisation
%==========================================================================
% We now repeat the procedure but, in this instance, starting with naïve
% agents whose likelihood mapping is imprecise (with small random additions
% to uninformative likelihood tensor). The idea here is to see if the
% likelihood tensors converge to a common frame of reference or ground;
% thereby enabling belief sharing.
%--------------------------------------------------------------------------
spm_figure('GetWin','language emergence'); clf

% set up naïve population
%--------------------------------------------------------------------------
clear MDP

NT    = 512;                             % number of episodes
ndp.T = 16;                              % number of epochs per episode
for m = 1:3                              % for each agent

    % likelihood mappings
    %----------------------------------------------------------------------
    MDP(m,1) = ndp;                      % initialise
    for i = 1:numel(Am)
        MDP(m,1).a{i} = Am{m}{i}*exp(16);
    end

    % initialise learnable (random) language mapping for all agents
    %----------------------------------------------------------------------
    for g = Cg
        MDP(m,1).a{g} = abs(randn(size(MDP(m).a{g}))) + 1;
        a0{m,g}       = MDP(m,1).a{g};
    end

end

%% experience-dependent learning over NT episodes
%==========================================================================
tp    = [1 NT/4 NT/2 NT];                % times to plot
tt    = 1;

k     = Cg(1);                           % modality to plot
KL    = [];                              % divergence among agents
FQ    = [];                              % Free energy (inference)
FA    = [];                              % Free energy (learning)
MI    = [];                              % Mutual information (likelihood)
CI    = [];                              % complexity (Dirichlet KL)
for t = 1:NT % for each exposure

    % report accumulation of Dirichlet counts
    %======================================================================
    if ismember(t,tp)
        for i = 1:numel(MDP)

            % plot likelihood tensor
            %--------------------------------------------------------------
            subplot(6,4,(i - 1)*4 + tt)
            imagesc(spm_dir_norm(MDP(i).a{k})), axis image
            if t == 1, ylabel(sprintf('agent %i',i)),   end
            if i == 1, title(sprintf('t = %i',tp(tt))), end

            % plot correlations among Dirichlet counts
            %--------------------------------------------------------------
            subplot(6,4,numel(MDP)*4 + tt)
            for m = (i + 1):numel(MDP)
                plot(MDP(i).a{k}(:),MDP(m).a{k}(:),'.'), axis square, hold on
                xlabel('Dirichlet counts')
            end

        end, drawnow, set(gca,'ColorOrderIndex',1);
        tt = tt + 1;
    end

    % exposure and evidence accumulation
    %======================================================================

    % active inference and learning
    %----------------------------------------------------------------------
    PDP = spm_MDP_VB_XXX(MDP,OPTIONS);

    % update Dirichlet parameters
    %----------------------------------------------------------------------
    MDP = spm_MDP_VB_update(MDP,PDP,OPTIONS);


    % report acquisition
    %======================================================================
    subplot(3,4,9), hold off
    for i = 1:numel(MDP)
        for j = (i + 1):numel(MDP)
            KL(t,i,j) = spm_KL_cat(MDP(i).a{k},MDP(j).a{k});
            plot(1:t,KL(:,i,j)), set(gca,'XLim',[1,NT]); hold on
        end
        MI(t,i) = spm_MDP_MI(MDP(i).a{k}(:,:));
        CI(t,i) = spm_KL_dir(MDP(i).a{k}(:,:),a0{i,k});
    end
    title({'KL divergence','(likelihood)'}),xlabel('time'),ylabel('nats')
    box off, axis square

    FQ(:,t) = sum(spm_cat({PDP.F}'),2);
    FA(:,t) = sum(spm_cat({PDP.Fa}'),2);

    subplot(3,4,10)
    plot(-FQ'), axis([1 NT 0 256])
    title({'Free energy','(inference)'}),xlabel('time'),ylabel('nats')
    box off, axis square

    subplot(3,4,11)
    plot(-FA'), axis([1 NT -2 2])
    title({'Free energy','(learning)'}),xlabel('time'),ylabel('nats')
    box off, axis square

    subplot(3,4,12)
    plot(MI), set(gca,'XLim',[1,NT])
    title({'Mutual information','(likelihood)'}),xlabel('time'),ylabel('nats')
    box off, axis square
    drawnow

end


%% Summarise learning in terms of joint free energy & structural complexity
%==========================================================================

spm_figure('GetWin','Generalised synchrony');
%--------------------------------------------------------------------------
subplot(3,2,1)
imagesc(FQ > -128)
title('Joint free energy','FontSize',14)
xlabel('episodes'),ylabel('agent')
axis square

subplot(3,2,2)
plot(CI)
title('Structural complexity','FontSize',14)
xlabel('episodes'),ylabel('nats')
axis square, spm_axis tight

return



% subroutines
%==========================================================================

function o = spm_what(s,c)
% returns an outcome from a list of attributes
% FORMAT o = spm_what(a,c)
% a   - latent states (attributes)
% c   - agent's line of sight
%
% {'what'}: outcomes
%--------------------------------------------------------------------------
%     'close friend',...
%     'close foe',...
%     'near person',...
%     'nothing'};
%--------------------------------------------------------------------------

% if there is an object or natural kind
%--------------------------------------------------------------------------
if s(1) == c

    % animate object:
    %     'proximity'  'pose'
    %----------------------------------------------------------------------
    if     s(2) == 1 && s(3) == 1
        o = 1;
    elseif s(2) == 1 && s(3) == 2
        o = 2;
    elseif s(2) == 2
        o = 3;
    else
        % nothing to see
        %------------------------------------------------------------------
        o = 4;
    end
else
    % nothing to see
    %----------------------------------------------------------------------
    o = 4;
end


function o = spm_contrast_energy(a,c)
% returns an outcome from a list of attributes
% FORMAT o = spm_contrast_energy(a,c)
% a   - latent states (attributes)
% c   - agent's line of sight
%
% {'where'}: outcomes
%--------------------------------------------------------------------------
% close ...1
% near  ...2
% none  ...3
%--------------------------------------------------------------------------

% if there is an object
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


function spm_surveillance_percept(MDP,c)
%% illustrates visual search graphically
%==========================================================================
% This subroutine creates little movies of each episode, for each agent
% illustrating what the agents saw and the — actual and inferred — latent
% states of the scene
%
% {'what'}: outcomes in DEM_scene
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

% {'what'}: outcomes in this demo
%--------------------------------------------------------------------------
%     'close friend',...        3
%     'close foe',...            4
%     'near person',...         9
%     'nothing'};               14

%--------------------------------------------------------------------------

% load images
%--------------------------------------------------------------------------
load DEM_scene
outcomes{1} = outcomes{1}(:,:,[3 4 9 14]);
k           = size(outcomes{1},1);
l           = size(outcomes{1},2);
spoken      = MDP(1).label.outcome;

% loop over agents
%--------------------------------------------------------------------------
Nm    = size(MDP,1);
LOS   = 9;
for m = 1:Nm

    % loop over time
    %----------------------------------------------------------------------
    for t = 1:MDP(m).T

        % what the agent actually sees
        %==================================================================
        % peripheral vision, or magnocellular (contrast energy)
        % foveal (central) vision, or parvocellular
        %------------------------------------------------------------------
        seen  = { ...
            outcomes{2}(:,:,MDP(m).o(2,t)) outcomes{2}(:,:,MDP(m).o(3,t)) outcomes{2}(:,:,MDP(m).o(4,t));
            outcomes{1}(:,:,end)           outcomes{1}(:,:,MDP(m).o(1,t)) outcomes{1}(:,:,end)};
        seen  = spm_cat(seen);

        subplot(4,Nm,m)
        image(spm_cat(seen))
        axis image, axis off

        str = [spoken{7}{MDP(m).o(7,t)} spoken{8}{MDP(m).o(8,t)} spoken{6}{MDP(m).o(6,t)}];

        text(32,32,str,'Color','r','FontSize',10)
        vision{m}(t) = getframe(gca);

        % actual scene
        %==================================================================
        subplot(4,1,2)
        for j = 1:LOS
            o        = spm_what(MDP(m).s(:,t),j);
            scene{j} = outcomes{1}(:,:,o);
        end

        % add direction of gaze and display
        %------------------------------------------------------------------
        image(spm_cat(scene))
        axis image, box off
        actual(t) = getframe(gca);


        % what the agent thinks she sees
        %==================================================================
        subplot(2*Nm,1,Nm + m)

        % find the most likely state of each object
        %------------------------------------------------------------------
        for f = 1:numel(MDP(m).X)
            [q,s] = max(MDP(m).X{f}(:,t));
            a(f)  = s;
        end

        % object in line of sight
        %------------------------------------------------------------------
        for j = 1:LOS
            o    = spm_what(a,j);
            scene{j} = outcomes{1}(:,:,o);
        end

        % add direction of gaze and display
        %------------------------------------------------------------------
        s    = c(m) + MDP(m).s(4,t) - 2;
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



function A  = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A           = bsxfun(@rdivide,A,sum(A,1));
A(isnan(A)) = 1/size(A,1);
