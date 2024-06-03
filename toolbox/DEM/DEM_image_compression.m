function [C,F,G,RDP,RGB] = DEM_image_compression(OPTIONS)
% Structure learning from pixels
%__________________________________________________________________________
%
% This routine demonstrates the mapping from pixels to (discrete)
% representations of episodes. It illustrates the use of a certain kind of
% deep generative model based upon reduction and grouping operators (RG)
% found in the renormalisation group. In brief, this allows an efficient
% encoding (and generation) of high dimensional content by maintaining
% conditionally independent partitions of hidden states in a recursive
% fashion.
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
rng(1)

try N = OPTIONS.N; catch, N = 4000; end           % Number of training data
try T = OPTIONS.T; catch, T = 1000; end           % Number of test data
try q = OPTIONS.q; catch, q = 64;   end           % concentration parameter  
try s = OPTIONS.s; catch, s = 2;    end           % smoothing in pixels  
try S = OPTIONS.S; catch, S = 24;   end           % histogram width

try Ne = OPTIONS.Ne; catch, Ne = 8;  end          % exemplars per class
try Ns = OPTIONS.Ns; catch, Ns = 10;  end         % Number of classes

try RGB.nd = OPTIONS.nd; catch, RGB.nd  = 5;  end % Diameter of tiles in pixels
try RGB.nb = OPTIONS.nb; catch, RGB.nb  = 9;  end % Number of discrete singular variates 
try RGB.mm = OPTIONS.mm; catch, RGB.mm  = 16; end % Maximum number of singular modes

RGB.su  = 16;                                     % Variance threshold
RGB.R   = 1;                                      % temporal resampling

% Load MNIST
%==========================================================================
cd 'C:\Users\karl\Dropbox\matlab'
load MNIST

% Read into working memory and show the movie
%==========================================================================
spm_figure('GetWin','Images'); clf

% Load images into RGB format suitable for videos
%--------------------------------------------------------------------------
for t = 1:(N + T)
    if t > N

        % test images
        %------------------------------------------------------------------
        j    = t - N;
        temp = imresize(test.images(:,:,j),[32,32]);

        % true class
        %------------------------------------------------------------------
        D{t} = full(sparse(test.labels(j) + 1,1,1,10,1));

    else

        % training images
        %------------------------------------------------------------------
        temp = imresize(training.images(:,:,t),[32,32]);

        % true class
        %------------------------------------------------------------------
        D{t} = full(sparse(training.labels(t) + 1,1,1,10,1));
    end

    temp    = spm_conv(temp,s,s);                 % convolve
    [~,i]   = sort(temp(:));                      % histogram equalisation
    n       = numel(i);                           % number of pixels
    temp(i) = exp(-((1:n) - n).^2/(S*n));

    temp    = min(max(temp,0),1);
    temp    = temp*255;
    
    I(:,:,1,t) = temp;
    I(:,:,2,t) = temp;
    I(:,:,3,t) = temp;
    
end

labels = training.labels;
clear training test

% marginal likelihood
%----------------------------------------------------------------------
j     = [];
for n = 1:Ns
    i = find(ismember(labels(1:N),n - 1));
    k = randperm(numel(i),Ne);
    j = [j;i(k)];
end

% And show an illustrative image
%--------------------------------------------------------------------------
subplot(2,2,1)
for i = 1:numel(j)
    spm_imshow(I(:,:,:,j(i)))
    axis on, title('MNIST image'), drawnow
end

% Map from image to discrete state space (c.f., Amortisation) 
%==========================================================================
[O,L,RGB] = spm_rgb2O(I,RGB);

% And show the images generated from a discrete representation
%--------------------------------------------------------------------------
subplot(2,2,2)
spm_imshow(spm_O2rgb(O(:,1),RGB))
axis on, title('Quantised image')

% exemplar images for (RG) structure learning
%--------------------------------------------------------------------------
MDP    = spm_MB_structure_learning(O(:,j),L,1);

% Add class contraints (i.e., number priors) as a final level
%--------------------------------------------------------------------------
No    = size(MDP{end}.B{1},1);             % Number of compressed images 
A     = zeros(No,Ns);
for n = 1:Ns
    i        = ismember(labels(j),n - 1);
    A(i,n)   = 1;
end

mdp.A{1}     = A;                          % likelihood (state)
mdp.A{2}     = ones(1 ,Ns);                % likelihood (paths)
mdp.B{1}     = ones(Ns,Ns);
mdp.id.A     = {1,1};
mdp.id.D     = {1};
mdp.id.E     = {2};
mdp.RG{1}    = 1;
mdp.LG       = mean(MDP{end}.LG);
mdp.L        = numel(MDP) + 1;
MDP{end + 1} = mdp;


% Show deep model structure in terms of (concatenated) likelihoods
%==========================================================================
spm_figure('GetWin','Paramters'); clf

% Show concatenated likelihoods at each level
%--------------------------------------------------------------------------
Nm    = numel(MDP);
for m = 2:Nm
    subplot(4,Nm,m - 1)
    Ns    = numel(MDP{m}.RG);
    A     = cell(Ns,Ns);
    for s = 1:Ns
        i      = MDP{m}.RG{s};
        A{s,s} = spm_cat(MDP{m}.A(i));
    end
    spm_spy(spm_cat(A)',4), axis square
end

% Show the same images but above the lower level they constrain
%--------------------------------------------------------------------------
for m = 1:Nm
    subplot(4,Nm,m + Nm)
    Ns    = numel(MDP{m}.RG);
    A     = cell(Ns,Ns);
    for s = 1:Ns
        i      = MDP{m}.RG{s};
        A{s,s} = spm_cat(MDP{m}.A(i));
    end
    spm_spy(spm_cat(A),4), axis square
    title(sprintf('Level %i',m))
end

% Illustrate compositionality in terms of changes in outputs with states
%==========================================================================
spm_figure('GetWin','Composition'); clf

% Create recursive model
%--------------------------------------------------------------------------
NDP   = spm_mdp2rdp(MDP);
for m = 1:Nm

    % size at this level
    %----------------------------------------------------------------------
    [Nf,Ns,Nu] = spm_MDP_size(NDP);
    NDP.T      = 1;

    % Get the generated image (O) under the first states at this level
    %----------------------------------------------------------------------
    for f = 1:Nf
        NDP.D{f} = sparse(1,1,1,Ns(f),1);
        NDP.E{f} = sparse(1,1,1,Nu(f),1);
    end
    PDP   = spm_MDP_VB_XXX(NDP);
    if m < Nm
        I = double(spm_O2rgb(PDP.Q.O{1}(:,1),RGB));
    else
        I = double(spm_O2rgb(PDP.O(:,1),RGB));
    end

    % Evaluate the change in the output for subsequent states
    %----------------------------------------------------------------------
    for s = 1:min(8,Ns(1))
        subplot(2*Nm,8,s + (m - 1)*8)
        NDP.D{1} = sparse(s,1,1,Ns(1),1);
        PDP      = spm_MDP_VB_XXX(NDP);
        if m == 1
            d = double(spm_O2rgb(PDP.Q.O{1}(:,1),RGB));
        elseif m < Nm
            d = double(spm_O2rgb(PDP.Q.O{1}(:,1),RGB)) - I;
        else
            d = double(spm_O2rgb(PDP.O(:,1),RGB)) - I;
        end
        imagesc(d(:,:,1)); axis image
        if s == 1, title(sprintf('Scale %i',m)), end
        drawnow
    end

    % Move one level down
    %----------------------------------------------------------------------
    try NDP = NDP.MDP; end
end

% Train: Parametric learning based upon a training dataset
%==========================================================================
% Learn by adding latent states if warranted by active selection to the
% latent factor style, under a supervised or precise prior over the latent
% factor (Digit) class
%--------------------------------------------------------------------------
spm_figure('GetWin','Train'); clf
train     = 1:N;
[RDP,G]   = spm_MNIST_train(MDP,RGB,O(:,train),D(train),1/q);

% Test: on unseen data
%==========================================================================

% Posterior predictive classification accuracy using a test set of data
%--------------------------------------------------------------------------
spm_figure('GetWin','Confusion'); clf
test      = (1:T) + N;
[C,F,PDP] = spm_MNIST_test(RDP,RGB,O(:,test),D(test));

% quantile analysis
%--------------------------------------------------------------------------
Fm     = median(F);
i      = F > median(F);
fprintf('Median F %.2f \n',Fm);
fprintf('Accuracy (total) %.2f \n',100*mean(C));
fprintf('Accuracy (upper median) %.2f \n',100*mean(C(i)));

% graphics for the last classification
%--------------------------------------------------------------------------
spm_figure('GetWin','Test'); clf
spm_show_RGB(PDP,RGB,1,false);

% Classification accuracy is a function of free energy
%==========================================================================
spm_figure('GetWin','Classification'); clf

subplot(2,2,1)
histogram(F(~~C),32), hold on
histogram(F(~C),32)
plot([Fm Fm],get(gca,'YLim'),'--k')
xlabel('ELBO'), ylabel('frequency')
title('Correct and incorrect classification'), axis square

f     = linspace(min(F),max(F),64);
for i = 1:numel(f)
    c(i) = mean(C(F > f(i)));
end
subplot(2,2,2)
plot(f,c*100), hold on
plot([Fm Fm],get(gca,'YLim'),'--k')

xlabel('ELBO threshold'), ylabel('classification accuracy')
title('Classification and confidence'), axis square


% Generate from recursive generative model
%==========================================================================
spm_figure('GetWin','Generation'); clf

% Create deep recursive model
%--------------------------------------------------------------------------
NDP       = RDP;
[~,Ns,Nu] = spm_MDP_size(NDP);

% Illustrate trained model in generative mode
%--------------------------------------------------------------------------
NDP.T = 1;
for i = 1:10
    subplot(8,10,i)
    NDP.D{1}  = sparse(i,1,1,Ns(1),1);
    NDP.E{1}  = sparse(1,1,1,Nu(1),1);
    PDP       = spm_MDP_VB_XXX(NDP);
    spm_imshow(spm_O2rgb(PDP.Q.O{1}(:,1),RGB));
end

fprintf('total samples %i\n',N)
try
    fprintf('effective samples %i\n',fix(sum(RDP.a{1},'all')))
end

% summary outputs
%--------------------------------------------------------------------------
C = mean(100*C);
F = mean(F);

return


% subroutines
%==========================================================================

function [RDP,F] = spm_MNIST_train(MDP,RGB,O,D,q)
% FORMAT [RDP,F] = spm_MNIST_train(MDP,RGB,O,D,q)
% MDP      - model struct (cell array): prior
% RGB      - image structure
% O        - train data
% D        - labels
% q        - concentration parameter
%
% RGB      - image structure
% RDP      - model struct (nested): posterior
% F        - ELBO
%
% This subroutine demonstrates learning of a recursive generative model. In
% the special case of static models, it is sufficient to place each
% training exemplar in the lowest level of the model — in the field MDP.O —
% and learn in the usual way by accumulating Dirichlet parameters (if the
% expected free energy improves). For dynamic models, sequences of outcomes
% or stimuli have to be specified — in the field MDP.S.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% default priors
%--------------------------------------------------------------------------
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 0;                             % suppress backwards pass
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.G = 1;                             % graphics
OPTIONS.P = 4;                             % graphics

if nargin < 5, q = 1/16; end               % concentration parameter

FIX.A  = 0;                                % enable likelihood learning
FIX.B  = 1;                                % but not transitions

% enable active learning (with minimal forgetting)
%--------------------------------------------------------------------------
N     = size(O,2);
Nm    = numel(MDP);
for m = 1:Nm
    MDP{m}.beta = 512;
    MDP{m}.eta  = N;
end

% train: with small concentration parameters (1/16)
%--------------------------------------------------------------------------
RDP      = spm_mdp2rdp(MDP,q,0,1,FIX);
RDP.A    = MDP{end}.A;
RDP      = rmfield(RDP,'a');

RDP.T = 1;
for j = 1:N

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    PDP   = spm_RDP_O(RDP,O(:,j),1);
    PDP.D = D(j);

    % active inference and learning
    %----------------------------------------------------------------------
    PDP   = spm_MDP_VB_XXX(PDP,OPTIONS);

    % update Dirichlet parameters
    %----------------------------------------------------------------------
    RDP   = spm_RDP_update(RDP,PDP);

    if OPTIONS.G

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(3,4,5)
        spm_imshow(spm_O2rgb(PDP.Q.O{1},RGB))
        [~,d] = max(PDP.D{1}); axis square
        title(sprintf('Content: %i',d - 1))

        subplot(3,4,6)
        spm_imshow(spm_O2rgb(PDP.Q.Y{1},RGB))
        [~,d] = max(PDP.X{1}); axis square
        title(sprintf('Content: %i',d - 1))

    end

    % mutual information of likelihood mapping
    %----------------------------------------------------------------------
    for l = 1:numel(PDP.Q.a)
        I(l,j) = spm_MDP_MI(PDP.Q.a{l});
    end
    F(j)       = PDP.Q.F + sum(PDP.F);

    subplot(3,2,2), hold off, plot(1:j,I(1:(end - 1),:))
    xlabel('number of exemplars'), ylabel('nats')
    title('Mutual information: first level[s]'), axis square
    subplot(3,2,4), hold off, plot(1:j,I(end,:))
    xlabel('number of exemplars'), ylabel('nats')
    title('Mutual information: final level'), axis square
    subplot(3,2,6), plot(1:j,F)
    xlabel('number of exemplars'), ylabel('nats')
    title('ELBO'), axis square
    try
        subplot(3,2,5), imagesc(RDP.a{1})
    catch
        subplot(3,2,5), imagesc(RDP.A{1})
    end
    title('Class likelihood'), axis square
    drawnow

    % disable graphics after first training iteration
    %----------------------------------------------------------------------
    if j > 2
        OPTIONS.G = 0;
        OPTIONS.P = 0;
    end
    % clc, disp(j), disp(I(:,j))
end


return

function [C,F,PDP] = spm_MNIST_test(RDP,RGB,O,D)
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

% uniform priors
%--------------------------------------------------------------------------
Ns       = size(RDP.B{1},1);
RDP.D{1} = ones(Ns,1);

% test
%--------------------------------------------------------------------------
N     = size(O,2);
RDP.T = 1;
for j = 1:N

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    PDP   = spm_RDP_O(RDP,O(:,[j j]),1);

    % active inference and learning
    %----------------------------------------------------------------------
    PDP   = spm_MDP_VB_XXX(PDP);

    % Bayesian model averaging
    %----------------------------------------------------------------------
    [~,d] = max(D{j});
    [~,p] = max(PDP.X{1}(:,end));

    % classification accuracy (and ELBO)
    %----------------------------------------------------------------------
    C(j)  = d == p;
    F(j)  = PDP.Q.F + sum(PDP.F);

    if N > 32
        % clc, disp(100*mean(C)), disp(j)
    end

    if ~C(j)

        spm_MDP_VB_trial(PDP)

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(3,4,5)
        spm_imshow(spm_O2rgb(PDP.Q.O{1}(:,end),RGB))
        title(sprintf('Content: %i',d - 1))

        subplot(3,4,6)
        spm_imshow(spm_O2rgb(PDP.Q.Y{1}(:,end),RGB))
        title(sprintf('Content: %i',p - 1))
        drawnow

    end

end

return
