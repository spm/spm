function mdp = DEM_dSprites
% Demo of active inference and structure learning (i.e.,disentanglement)
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
% inference to illustrate structure learning. Structure learning here is
% read as optimising the structure of a generative model that, crucially,
% includes dynamics. This foregrounds the sequential order and temporal
% scheduling of various updates. In this example, we start with a simple
% problem in which one or more objects can be removed around in a
% two-dimensional space. The kind of structure learning considered here can
% be likened to nonparametric Bayes; namely, the addition of a model
% component if licensed in terms of model evidence or marginal likelihood.
% Procedurally, this means that with each new stimulus (sequence of
% observations) various models are compared that entail the addition of a
% new latent state, path or factor. If the ELBO (i.e., negative variational
% free energy) increases the addition is accepted but not otherwise. The
% training sequences are carefully constructed to reflect the ordinal
% structure of observations. In other words, structure learning is
% predicated on both the content and dynamics of the generative process.
% 
% This demonstration calls a belief propagation scheme with factorisation
% of latent states into factors. Furthermore, the likelihood mapping is
% factorised into conditionally independent outcome modalities. This means
% that the size of any requisite tensor for belief updating is upper
% bounded by the factorisation or mean field approximation). This mitigates
% the van Neumann bottleneck; leading to increased efficiency at all three
% levels of optimisation (inference, learning and model selection).
% 
% A key aspect of this demonstration routine is that it deals with discrete
% state space generative models (and observations). This means that belief
% propagation and updating can be implemented using linear (tensor)
% operators; without worrying about nonlinearities of the sort found in
% continuous state space models.
%
%_________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8439 2023-03-27 18:41:45Z guillaume $

%% set up and preliminaries
%==========================================================================
rng(1)

% create three two-dimensional objects
%--------------------------------------------------------------------------
dSprite{1} = [ ...
    0 0 0 0 0 0 0 0;
    0 1 1 1 1 1 1 0;
    0 1 1 1 1 1 1 0;
    0 1 1 1 1 1 1 0;
    0 1 1 1 1 1 1 0;
    0 1 1 1 1 1 1 0;
    0 1 1 1 1 1 1 0;
    0 0 0 0 0 0 0 0];

dSprite{2} = [ ...
    0 0 0 0 0 0 0 0;
    0 0 1 1 0 1 1 0;
    0 1 1 1 1 1 1 1;
    0 1 1 1 1 1 1 1;
    0 0 1 1 1 1 1 0;
    0 0 0 1 1 1 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0];

dSprite{3} = [ ...
    0 0 0 0 0 0 0 0;
    0 0 1 1 1 1 0 0;
    0 1 1 1 1 1 1 0;
    1 1 1 1 1 1 1 1;
    1 1 1 1 1 1 1 1;
    0 1 1 1 1 1 1 0;
    0 0 1 1 1 1 0 0;
    0 0 0 0 0 0 0 0];


% Latent states
%--------------------------------------------------------------------------
% First, specify latent states
%--------------------------------------------------------------------------
% This sets up a label structure with the names and cardinality of each
% factor and the accompanying latent states. These are further equipped
% with a number of actions that require a probability transition matrix for
% each allowable action. Here, the actions correspond to moving things in
% each dimension
%--------------------------------------------------------------------------
label.factor  = {'x','y','kind'};

Nx    = 8;                          % number of pixels
Nk    = numel(dSprite);             % number of dSprite kinds

for i = 1:Nx, label.name{1}{i} = sprintf('x%i',i); end
for i = 1:Nx, label.name{2}{i} = sprintf('y%i',i); end
for i = 1:Nk, label.name{3}{i} = sprintf('k%i',i); end

label.action{1} = {'stay','right','left'};
label.action{2} = {'stay','down','up'};
label.action{3} = {'none'};

% visual output modalities
%--------------------------------------------------------------------------
g     = 1;
for i = 1:Nx
    for j = 1:Nx
        label.modality{g} = sprintf('v(%i,%i)',i,j);
        label.outcome{g} = {'on','off'};
        g = g + 1;
    end
end

% and reward
%--------------------------------------------------------------------------
label.modality{g} = 'reward';
label.outcome{g}  = {'on','off'};


%% Transitions: B
%==========================================================================

% Transitions among these states are characterised by movement
%--------------------------------------------------------------------------
disp('specifying generative process'), disp(' ')
WRAP        = 1;
nx          = numel(label.name{1});
B{1}(:,:,1) = full(spm_speye(nx,nx,-0,WRAP));
B{1}(:,:,2) = full(spm_speye(nx,nx,-1,WRAP));
B{1}(:,:,3) = full(spm_speye(nx,nx,+1,WRAP));

nx          = numel(label.name{2});
B{2}(:,:,1) = full(spm_speye(nx,nx,-0,WRAP));
B{2}(:,:,2) = full(spm_speye(nx,nx,-1,WRAP));
B{2}(:,:,3) = full(spm_speye(nx,nx,+1,WRAP));

nx          = numel(label.name{3});
B{3}(:,:,1) = full(spm_speye(nx,nx,0,2));

for f = 1:numel(B)
    B{f}(:) = bsxfun(@rdivide,B{f}(:,:),sum(B{f}(:,:),1));
    Nf(f)   = size(B{f},1);
    Nu(f)   = size(B{f},3);
end



%% outcome probabilities: A
%==========================================================================
% Next, we specify the probabilistic mappings between latent states and
% outcomes, with a tensor for each outcome modality, which could model the
% high order interactions among the causes of outcomes (e.g., occlusion).
%--------------------------------------------------------------------------
for g = 1:numel(label.modality)

    A{g} = single(zeros([numel(label.outcome{g}),Nf]));

    % default outcomes (absent)
    %----------------------------------------------------------------------
    A{g}(2,:,:,:) = 1;
end

%  spatial targets for different kinds  of objects
%--------------------------------------------------------------------------
target{1} = [3,1]*Nx/4;
target{2} = [3,3]*Nx/4;
target{3} = [1,2]*Nx/4;
precision = eye(2,2)/(Nx/4)^2;

% for every combination  of latent states, specify an outcome
%--------------------------------------------------------------------------
wrap  = @(i,x)mod( i - 1,x) + 1;             % torus grid
for x = 1:Nf(1)
    for y = 1:Nf(2)
        for k = 1:Nf(3)

            % reward outcome
            %==============================================================
            d = [x,y] - target{k};           % distance from target
            r = exp(-d*precision*d'/2);      % probability of reward

            A{end}(:,x,y,k) = double([r; 1 - r]);

            %  visual ouutcomes
            %==============================================================

            % location of stimuli (subscript)
            %--------------------------------------------------------------
            [i,j] = find(dSprite{k});
            i     = i + x - 4;
            j     = j + y - 5;

            % location of stimuli (indices)
            %--------------------------------------------------------------
            g     = sub2ind(Nf(1:2),wrap(i,Nf(1)),wrap(j,Nf(2)));

            % place in likelihood tensor
            %--------------------------------------------------------------
            for o = g'
                A{o}(:,x,y,k) = [1; 0];
            end

        end
    end
end



%% priors: (cost) C
%--------------------------------------------------------------------------
% Finally, specify the prior preferences in terms of log probabilities over
% outcomes:

% uninformative preferences over vision
%--------------------------------------------------------------------------
for g = 1:numel(A)
    C{g} = zeros(1,size(A{g},1));
end

% with a preference for rewarding outcomes
%--------------------------------------------------------------------------
C{end}   = [128,0];

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify (uninformative) prior
% beliefs about initial states (D) and paths through those states (D)
%--------------------------------------------------------------------------
for f = 1:numel(B)
    D{f} = ones(Nf(f),1);
    E{f} = ones(Nu(f),1);
end

% specify controllable factors; here, the first two factors
%--------------------------------------------------------------------------
U     = [1 1 0];                  % controlable factors

% MDP Structure, specifying 64 epochs (i.e., 16 seconds of active vision)
%==========================================================================
mdp.T = 16;                       % numer of moves
mdp.U = U;                        % controllable factors
mdp.A = A;                        % likelihood probabilities
mdp.B = B;                        % transition probabilities
mdp.C = C;                        % prior preferences
mdp.D = D;                        % prior over initial states
mdp.E = E;                        % prior over initial paths

mdp.label = label;

% Solve an example with known (veridical) structure and parameters to
% establish – and illustrate – optimal performance
%==========================================================================
disp('inverting generative model (c.f., active inference)'), disp(' ')

% specify a trial, with initial conditions, s
%--------------------------------------------------------------------------
OPTIONS.N = 1;                    % simulate neuronal responses
MDP   = mdp;
MDP.s = [1,1,2];
MDP   = spm_MDP_VB_XXX(MDP,OPTIONS);


% illustrate behavioural responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1));

% illustrate physiological responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP(1),[],3);

% illustrate action and perceptual inference
%--------------------------------------------------------------------------
spm_figure('GetWin','optimal behaviour'); clf
spm_report(MDP)


% Illustrate structure learning via assimilation of epochs of observations
%==========================================================================
mdp = spm_MDP_structure_learning(mdp);


% Illustrate mmotor babbling (motor learning)
%==========================================================================
mdp = spm_MDP_motor_learning(mdp);

% active learning
%==========================================================================
% So far, only a small number of likelihood parameters have been learned.
% We can now call upon the active learning aspect of active inference
% (i.e., novelty) to sample inputs efficiently in a way that reduces
% uncertainty about the likelihood mapping (i.e., fill in the gaps). This
% can be done efficiently because there are now precise beliefs about
% latent states and their dynamics. In other words, provided the agent
% knows the initial states and paths it can associate any combination of
% states with the most likely outcome; in the spirit of one shot learning
%--------------------------------------------------------------------------
spm_figure('GetWin','learning'); clf
disp('inverting generative model (c.f., active inference)'), disp(' ')

OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.G = 1;                             % suppress graphics

MDP   = mdp;
MDP.T = 512;
MDP.N = 0;
for d = 1:numel(dSprite)

    % initial conditions for this object
    %----------------------------------------------------------------------
    MDP.s = [1,1,d];

    % active inference and learning
    %----------------------------------------------------------------------
    PDP   = spm_MDP_VB_XXX(MDP,OPTIONS);   % invert
    MDP.a = PDP.a;                         % update likelihood parameters


    % report graphically
    %==================================================================
    if OPTIONS.G

        spm_figure('GetWin','inference'); clf
        spm_MDP_VB_trial(PDP);

        spm_figure('GetWin','behaviour'); clf
        spm_report(PDP)

        % illustrate learning over trials
        %----------------------------------------------------------------------
        spm_figure('GetWin','learning');

        subplot(4,1,1), plot(PDP.F)             % free energy (states)
        xlabel('trials'), ylabel('nats')
        title('Variational free energy (ELBO: states)','FontSize',14)
        set(gca,'XLim',[1 (MDP.T - 16)]), hold on

        subplot(4,1,2), plot(PDP.v)             % expected free energy
        xlabel('trials'), ylabel('nats')
        title('Expected free energy','FontSize',14)
        set(gca,'XLim',[1 (MDP.T - 16)]), hold on

        subplot(4,1,3), plot(PDP.w)             % precision over policies
        xlabel('trials'), ylabel('nats')
        title('Precision (policies)','FontSize',14)
        set(gca,'XLim',[1 (MDP.T - 16)]), hold on


        % reward outcomes
        %----------------------------------------------------------------------
        subplot(4,1,4), plot(conv(PDP.o(end,:) == 1, ones(32,1)/32,'same'))
        xlabel('trials'), ylabel('number')
        title('Performance (rewards)','FontSize',14)
        set(gca,'XLim',[1 (MDP.T - 16)]), hold on
        drawnow

    end

end

%--------------------------------------------------------------------------
% The graphics illustrate the emergent transition from epistemic
% (exploratory) behaviour to preference seeking (exploitative) behaviour.
% This as evinced by a reduction in expected free energy as expected
% information gain (about parameters) is resolved through active learning;
% revealing the preference seeking – that itself depends upon learning the
% likelihood mapping from latent states to the preferred reward channel.
%__________________________________________________________________________


return


% subroutines
%==========================================================================

function spm_report(MDP)

%% illustrates active inference graphically
%--------------------------------------------------------------------------
for f = 1:numel(MDP(1).B)
    Nf(f) = size(MDP.B{f},1);
end

% show (initial) stimulus and percept in outcome space
%==========================================================================
for g = 1:(numel(MDP.A) - 1)
    o(g)     = MDP.o(g,1);
    if isfield(MDP,'a')
        A{g} = spm_dir_norm(MDP.a{g});
    else
        A{g} = MDP.A{g};
    end
    for f = 1:numel(Nf)
        Q{f} = MDP.X{f}(:,1);
    end
    QO   = spm_dot(A{g},Q);
    O(g) = QO(2);
end

subplot(3,2,1), imagesc(reshape(o,Nf(1),Nf(2)))
axis image, title('Stimulus','FontSize',14)
subplot(3,2,2), imagesc(reshape(O,Nf(1),Nf(2)))
axis image, title('Percept' ,'FontSize',14)
subplot(3,2,1), hold on
plot(MDP.s(2,:),MDP.s(1,:),'r.','MarkerSize',8), hold off

% show trajectory
%--------------------------------------------------------------------------
for t = 1:MDP.T

    % what the agent perceives
    %======================================================================
    for g = 1:(numel(MDP.A) - 1)
        O(g) = MDP.Y{g,t}(2,:);
    end

    subplot(3,2,2), hold off
    imagesc(reshape(O,Nf(1),Nf(2))), axis image
    vision(t) = getframe(gca);

end

% assign movies to graph object
%--------------------------------------------------------------------------
subplot(3,2,2)
set(gca,'Userdata',{vision,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Percept','FontSize',14)

% show likelihood mapping to reward
%--------------------------------------------------------------------------
for f = 1:Nf(3)
    if isfield(MDP,'a')
        R = spm_dir_norm(MDP.a{end}(:,:,:,f));
    else
        R = MDP.A{end}(:,:,:,f);
    end
    R = squeeze(R(1,:,:));
    subplot(6,3,f + 2*3)
    imagesc(R), axis image, title(sprintf('Reward mapping (%i)',f))
    
end

% show Dirichlet counts
%--------------------------------------------------------------------------
for f = 1:Nf(3)
    if isfield(MDP,'a')
        R = sum(MDP.a{end}(:,:,:,f),1);
    else
        R = sum(MDP.A{end}(:,:,:,f),1);
    end
    R = squeeze(R(1,:,:));
    subplot(6,3,f + 3*3)
    image(64*R/(MDP.T/prod(Nf(1:2))))
    axis image, title(sprintf('Dirichlet counts (%i)',f))
end

% place fields - first factor (i.e., marginal likelihoods)
%--------------------------------------------------------------------------
Np    = 4;
for i = 1:Np
    for g = 1:numel(A)
        P(g) = sum(A{g}(1,i,:,:),[3,4]);
    end
    subplot(6,Np,i + Np*4)
    imagesc(reshape(P,Nf(1),Nf(2))), axis image
    title(sprintf('Place fields (%i,%i)',1,i))

end

% place fields - second factor
%--------------------------------------------------------------------------
Np    = 4;
for i = 1:Np
    for g = 1:numel(A)
        P(g) = sum(A{g}(1,:,i,:),[2,4]);
    end
    subplot(6,Np,i + Np*5)
    imagesc(reshape(P,Nf(1),Nf(2))), axis image
    title(sprintf('Place fields (%i,%i)',2,i))
end
drawnow

return


%% NOTES: on factorisation of Dirichlet likelihood tensors
%--------------------------------------------------------------------------
clear K L
for i = 1:256
    for j = 1:size(A{i},1)
        L{i}(j,:) = A{i}(1,:);
        m    = spm_marginal(A{i}(1,:,:,:));
        K{i}(j,:) = spm_vec(spm_cross(m(2:end)))';
    end
end
L = spm_cat(L');
K = spm_cat(K');


%% routines that call XXX
%--------------------------------------------------------------------------
DEMO_MDP_maze_X
DEM_surveillance
DEM_dSprites
DEM_syntax





