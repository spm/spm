function mdp = DEM_MNIST
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
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $

%% set up and preliminaries
%==========================================================================
rng(1)

%% Load MNIST dataset
%--------------------------------------------------------------------------
load('mnist.mat');

% initialise learnable language mapping for all agents
%--------------------------------------------------------------------------
mdp.p = exp(-0);%%%%
Ns    = [32,10];
Nf    = numel(Ns);
Ng    = numel(spm_MNIST2o(training,1));

for g = 1:Ng
    mdp.a{g} = zeros([2,Ns]) + mdp.p;
end
for f = 1:Nf
    mdp.b{f} = eye(Ns(f),Ns(f));
end

% illustrate behavioural responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf

for n = 1:Ns(2)

    j     = find(training.labels == (n - 1));
    j     = j(randperm(numel(j)));
    for u = 1:Ns(1)
        O     = spm_MNIST2o(training,j(u));
        for g = 1:Ng
            mdp.a{g}(:,u,n) = O{g} + mdp.p;
        end
        for g = 1:Ng
            Y{g} = spm_dir_norm(mdp.a{g}(:,u,n));
        end

        subplot(2,3,1)
        imagesc(spm_o2MNIST(O))
        axis square, drawnow

    end
end

% default priors
%--------------------------------------------------------------------------
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 1;                             % suppress explicit action
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.O = 1;                             % probabilistic outcomes
OPTIONS.G = 1;                             % suppress graphics
OPTIONS.P = 1;                             % suppress graphics


[Nf,Ns,Nu] = spm_MDP_size(mdp);
for f = 1:Nf
    mdp.D{f} = ones(Ns(f),1);              % initial state
    mdp.E{f} = ones(Nu(f),1);              % initial control
end

%  train
%--------------------------------------------------------------------------
mdp.T = 1;
mdp.U = zeros(1,Nf);
for j = 1:10000

    [O,D]    = spm_MNIST2o(training,j);

    pdp      = mdp;
    pdp.O    = O;
    pdp.D{2} = D;
    pdp      = spm_MDP_VB_XXX(pdp,OPTIONS);

    subplot(3,4,5)
    imagesc(spm_o2MNIST(O))
    [m,d] = max(D); axis square
    title(sprintf('Content: %i',d - 1))

    subplot(3,4,6)
    imagesc(spm_o2MNIST(pdp.Y))
    [m,d] = max(pdp.X{2}); axis square
    title(sprintf('Content: %i',d - 1))

    % Bayesian model reduction
    %----------------------------------------------------------------------
    for g = 1:numel(mdp.a)
%        mdp.a{g} = spm_MDP_VB_prune(pdp.a{g},mdp.a{g},0,-2,0,'MI');
         mdp.a{g} = pdp.a{g};
    end

    % mutual information of likelihood
    %----------------------------------------------------------------------
    I     = 0;
    for g = 1:numel(mdp.a)
        I = I + spm_MDP_MI(mdp.a{g});
    end
    MI(j) = I;
    FI(j) = pdp.F;

    subplot(3,2,2)
    plot(MI), axis square
    subplot(3,2,4)
    plot(FI), axis square, drawnow

    subplot(3,2,6)
    ma    = spm_marginal(mdp.a{1});
    for i = 1:Nf
         bar(ma{i + 1}), hold on
    end

    save mdp mdp

end

OPTIONS.B = 0;                             % suppress explicit action

%     % Bayesian model reduction
%     %----------------------------------------------------------------------
%     for g = 1:numel(mdp.a)
%         mdp.a{g} = spm_MDP_VB_prune(pdp.a{g},mdp.a{g},0,0,0,'SIMPLE');
%     end

%  train
%--------------------------------------------------------------------------
for j = 1:100

    [O,D] = spm_MNIST2o(test,j);

    pdp   = mdp;
    pdp.O = O;
    pdp   = spm_MDP_VB_XXX(pdp,OPTIONS);


    subplot(3,4,5)
    imagesc(spm_o2MNIST(O)), axis image
    [m,d] = max(D);
    title(sprintf('Content: %i',d - 1))

    subplot(3,4,6)
    imagesc(spm_o2MNIST(pdp.Y)), axis image
    [m,p] = max(pdp.X{2});
    title(sprintf('Prediction: %i',p - 1))
    drawnow

    C(j) = d == p;
    disp(100*mean(C))

end

return

%  illustrate
%--------------------------------------------------------------------------
n     = 7;
for u = 1:Ns(1)
    for g = 1:Ng
        Y{g} = spm_dir_norm(mdp.a{g}(:,u,n));
    end
    YY{u} = spm_o2MNIST(Y);
end
subplot(4,1,1)
imagesc(spm_cat(reshape(YY,2,16)))
axis image, axis off
drawnow


% NOTES: 
%==========================================================================


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



% subroutines
%==========================================================================

function [O,D] = spm_MNIST2o(mnist,id)
% FORMAT [O,D] = spm_MNIST2o(mnist,id)
% mnist - struct
% id     - content
%
% o     - observation
% D     - label (prior)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% create array of observations: each column is a word (i.e., letter sequence)
%--------------------------------------------------------------------------
for i = 1:mnist.height
    for j = 1:mnist.width
        o      = mnist.images(i,j,id);
        O{i,j}(1,1) = o;
        O{i,j}(2,1) = 1 - o;
    end
end

% return a cell for this image
%--------------------------------------------------------------------------
O = O(:);
D = sparse(mnist.labels(id) + 1,1,1,10,1);

return


function [I] = spm_o2MNIST(Y)
% FORMAT [I] = spm_o2MNIST(Y)
% Y      - posterior prediction

% I      - image
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% create array of observations: each column is a word (i.e., letter sequence)
%--------------------------------------------------------------------------
for i = 1:numel(Y)
    I(i) = Y{i}(1);
end
I     = reshape(I,28,28);

return



