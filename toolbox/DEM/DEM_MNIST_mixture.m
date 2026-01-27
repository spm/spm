function [RDP,RGB] = DEM_MNIST_mixture
% Structure learning from pixels
%__________________________________________________________________________
%
% This example uses image classification to showcase the ability of
% recursive generative models to compress images in a way that leverages
% the composition of small image components into larger arrangements. Uses
% the MNIST digit classification problem. In brief, images are discretised
% using singular value decomposition of local image patches. These are then
% recursively assembled using a fast structure learning routine. Subsequent
% training images and then assimilated using parametric learning. The
% ensuing generative model is then illustrated in terms of digit
% classification and the underlying architecture.
%
% This routine is used to demonstrate the difficulties in using variational
% free energy for classification under mixtures of experts (i.e.,
% renormalizable generative models) with different parameterisations.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% This subroutine expects to find MNIST.mat in its working directory. The
% requisite MATLAB file can be downloaded from the following website
%--------------------------------------------------------------------------
% https://lucidar.me/en/matlab/load-mnist-database-of-handwritten-digits-in-matlab/


% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%--------------------------------------------------------------------------
rng(1), close all

% Load MNIST
%==========================================================================
cd 'C:\Users\karl\Dropbox\matlab'
load MNIST

% Read into working memory and show the movie
%==========================================================================
spm_figure('GetWin','Images'); clf

% Load images into RGB format suitable for videos
%--------------------------------------------------------------------------
N     = 4000;                                     % size of training set
T     = 100;                                      % size of test set
S     = 24;                                       % histogram width
for t = 1:(N + T)
    temp    = imresize(training.images(:,:,t),[32,32]);
    [~,i]   = sort(temp(:));                      % histogram equalisation
    n       = numel(i);                           % number of pixels
    temp(i) = exp(-((1:n) - n).^2/(S*n));

    temp    = min(max(temp,0),1);
    temp    = temp*255;
    
    I(:,:,1,t) = temp;
    I(:,:,2,t) = temp;
    I(:,:,3,t) = temp;
    
    % true class
    %----------------------------------------------------------------------
    D{t}       = full(sparse(training.labels(t) + 1,1,1,10,1));
    d(t)       = training.labels(t);
end

% And show an illustrative image
%--------------------------------------------------------------------------
subplot(2,2,1)
spm_imshow(I(:,:,:,1))
title('MNIST image')

% Map from image to discrete state space (c.f., Amortisation) 
%==========================================================================
RGB.nd  = 8;                     % Diameter of tiles in pixels
RGB.nb  = 7;                     % Number of discrete singular variates 
RGB.mm  = 16;                    % Maximum number of singular modes
RGB.su  = 8;                     % Variance threshold
RGB.R   = 1;                     % temporal resampling
[O,L,RGB] = spm_rgb2O(I,RGB);

% And show the images generated from a discrete representation
%--------------------------------------------------------------------------
subplot(2,2,2)
spm_imshow(spm_O2rgb(O(:,1),RGB))
drawnow

% Use a small number (128) for (RG) structure learning
%--------------------------------------------------------------------------
No    = 8;                                 % Number of compressed images
Ns    = 10;                                % Number of classes
for n = 1:Ns

    % Train: Structure learning based upon initla samples
    %----------------------------------------------------------------------
    train    = find(ismember(d,n - 1));
    mdp      = spm_MB_structure_learning(O(:,train(1:No)),L,1);

    % Train: Parametric learning based upon a training dataset
    %----------------------------------------------------------------------
    spm_figure('GetWin','Train'); clf
    MDP(n,1) = spm_MNIST_train(mdp,RGB,O(:,train));

end

% Test: on unseen data
%==========================================================================

% Posterior predictive classification accuracy using a test set of data
%--------------------------------------------------------------------------
spm_figure('GetWin','Confusion'); clf
test   = (1:T) + N;
[C,F]  = spm_MNIST_test(MDP,RGB,O(:,test),D(test));

% classifiable test data
%--------------------------------------------------------------------------
[F,i]  = sort(F);
C      = C(i);
Nc     = fix(numel(F)/2);
fprintf('Accuracy (total) %.2f \n',100*mean(C));
fprintf('Accuracy (upper) %.2f \n',100*mean(C((1:Nc) + Nc)));

% Classification accuracy is a function of free energy
%==========================================================================
spm_figure('GetWin','Classification'); clf

subplot(2,2,1)
histogram(F(~~C),32), hold on
histogram(F(~C),32)
xlabel('ELBO'), ylabel('frequency')
title('Correct and incorrect classification'), axis square

f     = linspace(min(F),max(F),64);
for i = 1:numel(f)
    c(i) = mean(C(F > f(i)));
end
subplot(2,2,2)
plot(f,c*100)
xlabel('ELBO threshold'), ylabel('classification accuracy')
title('Classification and confidence'), axis square


return


% subroutines
%=========================================================================

function mdp = spm_MNIST_train(MDP,RGB,O)
% FORMAT RDP = spm_MNIST_train(MDP,RGB,O)
% MDP      - model struct (cell array): prior
% RDP      - model struct (nested): posterior
% O        - train data
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% default priors
%--------------------------------------------------------------------------
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 0;                             % suppress backwards pass
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.G = 1;                             % graphics
OPTIONS.P = 4;                             % graphics

FIX.A     = 0;                             % enable likelihood learning
FIX.B     = 1;                             % but not transitions

% upper bound on mutual information (for plotting)
%--------------------------------------------------------------------------
U     = zeros(1,numel(MDP));
for l = 1:numel(MDP)
    for g = 1:numel(MDP{l}.A)
        U(l) = U(l) + log(size(MDP{l}.A{g},2))/2;
    end
end

% enable active learning (with minimal forgetting)
%--------------------------------------------------------------------------
for m = 1:numel(MDP)
    MDP{m}.beta = 512;
    MDP{m}.eta  = 512;
end

% train: with small concentration parameters (1/16)
%--------------------------------------------------------------------------
MDP   = spm_mdp2rdp(MDP,1/16,0,1,FIX);
mdp   = MDP;
mdp.T = 1;
for j = 1:size(O,2)

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    pdp   = spm_RDP_O(mdp,O(:,j),1);

    % active inference and learning
    %----------------------------------------------------------------------
    pdp   = spm_MDP_VB_XXX(pdp,OPTIONS);

    % update Dirichlet parameters
    %----------------------------------------------------------------------
    mdp   = spm_RDP_update(mdp,pdp);

    if OPTIONS.G

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(3,4,5)
        spm_imshow(spm_O2rgb(pdp.Q.O{1},RGB))
        axis square
        title(sprintf('Content: %i',j - 1))

        subplot(3,4,6)
        spm_imshow(spm_O2rgb(pdp.Q.Y{1},RGB))
        axis square
        title(sprintf('Content: %i',j - 1))

    end

    % mutual information of likelihood mapping
    %----------------------------------------------------------------------
    for l = 1:numel(pdp.Q.a)
        I(l,j) = spm_MDP_MI(pdp.Q.a{l});
    end
    I(l + 1,j) = spm_MDP_MI(pdp.a);
    F(j)       = pdp.Q.F + sum(pdp.F);

    subplot(3,2,2), hold off, plot(1:j,I), hold on
    plot(1:j,U(1)*ones(1,j),'--')
    xlabel('number of exemplars'), ylabel('nats')
    title('Mutual information: first level'), axis square
    subplot(3,2,4), hold off, plot(1:j,I(end,:)), hold on
    plot(1:j,U(end)*ones(1,j),'--')
    xlabel('number of exemplars'), ylabel('nats')
    title('Mutual information: final level'), axis square
    subplot(3,2,6), plot(1:j,F)
    xlabel('number of exemplars'), ylabel('nats')
    title('ELBO'), axis square
    drawnow

    % disable graphics after intial training iterations
    %----------------------------------------------------------------------
    if j > (2*size(mdp.B{1},1))
        OPTIONS.G = 0;
        OPTIONS.P = 0;
    end
    clc, disp(j), disp(I(:,j))
end


return

function [C,G,pdp] = spm_MNIST_test(RDP,RGB,O,D)
% FORMAT [C,F,PDP] = spm_MNIST_test(RDP,RGB,O,D)
% RDP   - model struct
% O     - test data
% D     - class labels for accuracy
%
% C     - correct classification vector
% F     - ELBO (log marginal likelihood)
% PDP   - posterior RDP
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging


%  train
%--------------------------------------------------------------------------
C     = 0;
for j = 1:size(O,2)

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    rdp   = spm_RDP_O(RDP,O(:,[j j]),1);

    % active inference and learning
    %----------------------------------------------------------------------
    for m = 1:numel(rdp)
        pdp(m) = spm_MDP_VB_XXX(rdp(m));
        F(j,m) = pdp(m).F + pdp(m).Q.F;
    end

    % Bayesian model selection
    %----------------------------------------------------------------------
    [~,d] = max(D{j});
    [f,m] = max(F(j,:));

    % classification accuracy (and ELBO)
    %----------------------------------------------------------------------
    C(j)  = d == m;
    G(j)  = f;
    clc, disp(100*mean(C)), disp(j)

     if ~C(j)

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(3,4,5)
        spm_imshow(spm_O2rgb(pdp(m).Q.O{1}(:,end),RGB))
        title(sprintf('Content: %i',d - 1))

        subplot(3,4,6)
        spm_imshow(spm_O2rgb(pdp(m).Q.Y{1}(:,end),RGB))
        title(sprintf('Content: %i',m - 1))
        drawnow

    end

end

return


