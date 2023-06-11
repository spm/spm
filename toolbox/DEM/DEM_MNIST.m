function mdp = DEM_MNIST
% Demo of structure learning (i.e.,disentanglement) of MNIST digits
%__________________________________________________________________________
%
% This routine provides a simple illustration of disentanglement or
% structure learning using a static (image) recognition problem; namely the
% MNIST digits. It casts the structure learning problem as a data
% assimilation problem by accumulating latent states as discrete levels of
% a factor in a hidden Markov model with two factors — the second factor
% being the digit class. The key procedure demonstrated here is the
% addition of a new latent (style) state using Bayesian model comparison
% under the prior that additional states will increase expected free
% energy. In the absence of any prior costs or constraints; this reduces to
% accepting a new training exemplar if adding Dirichlet counts to a new
% state (i.e., Column of the likelihood tensor) decreases variational and
% expected free energy. In other words, increases the marginal likelihood
% of the new observation under the prior constraint that the resulting
% likelihood mapping has a greater mutual information.
% 
% In this example, an upper bound of 128 styles is implemented as an
% additional hyperprior and the model is learned during exposure to 1024
% training exemplars for each digit class. The ensuing classification
% accuracy is established using a test dataset of 10,000. Performance is
% between 95 and 100%; depending upon whether one accepts all the test data
% or just those with a sufficiently high marginal likelihood (as scored by
% the variational free energy).
% 
% The code below contains a number of subroutines and alternative
% formulations of structure learning. It also examines the information
% geometry inherent in the discovered structures. This routine can take
% some time to ingest the training data (and will require downloading of a
% MATLAB file as described in the preparation subroutine).
% 
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $

%% set up and preliminaries
%==========================================================================

% Load (pre-processed) MNIST dataset
%--------------------------------------------------------------------------
spm_MNIST_prepare;
load spm_MNIST

% train
%==========================================================================
% Learn by adding latent states if warranted by active selection to the
% latent factor style, under a supervised or precise prior over the latent
% factor (Digit) class
%--------------------------------------------------------------------------
NS    = 128;                                % number of styles
NT    = 1024;                               % number of training digits
mdp   = spm_MNIST_learn(training,NS,NT);    % supervised structure learning

% Illustrate learned latent states as images
%--------------------------------------------------------------------------
spm_figure('GetWin','Disentangled');
spm_show_a(mdp.a);


% test (e.g., accuracy 95.1 %)
%==========================================================================
% Posterior predictive classification accuracy using a test set of data
%--------------------------------------------------------------------------
[C,F] = spm_MNIST_test(mdp,test,10000);


% graphics: mutual information, free energy and classification accuracy
%==========================================================================
spm_figure('GetWin','Test'); clf

subplot(2,2,1)
histogram(F(C),32), hold on
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

% illustrate images with a high and low marginal likelihood
%--------------------------------------------------------------------------
for n = 1:10
    j    = find(test.labels == n - 1);
    j    = j(C(j));
    j    = j(F(j) == max(F(j)));
    I{n} = spm_o2MNIST(spm_MNIST2o(test,j));
end
subplot(8,1,5), imagesc(spm_cat(I)), axis image off
title('Platonic digits','FontSize',14)
for n = 1:10
    j    = find(test.labels == n - 1);
    j    = j(~C(j));
    j    = j(F(j) == min(F(j)));
    I{n} = spm_o2MNIST(spm_MNIST2o(test,j));
end
subplot(8,1,6), imagesc(spm_cat(I)), axis image off
title('Unlikely digits','FontSize',14)

% Bayesian model reduction
%==========================================================================
% Illustrate Bayesian model reduction by simply setting tensor elements
% with small values zero
%--------------------------------------------------------------------------
% Number of parameters before and after reduction
%    1253376
%    1194895
% ELBO before and after reduction (per modality)
%   -5.3413e+03
%   -5.3255e+03
%    15.8223
%--------------------------------------------------------------------------
rdp   = mdp;
Ng    = numel(mdp.a);
for g = 1:Ng
    rdp.a{g} = spm_MDP_VB_prune(mdp.a{g},mdp.p,0,0,0,'SIMPLE');
end

disp('number parameters before and after reduction')
disp(sum(~~spm_vec(mdp.a)))
disp(sum(~~spm_vec(rdp.a)))

% Evaluate the predictive accuracy of the reduced model
%--------------------------------------------------------------------------
[Cr,Fr] = spm_MNIST_test(rdp,test,10000);

disp('ELBO before and after reduction (per modality)')
disp(sum(F)/Ng)
disp(sum(Fr)/Ng)
disp(sum(Fr)/Ng - sum(F)/Ng)

% graphics: increasing classification accuracy with marginal likelihood
%--------------------------------------------------------------------------
f     = linspace(min(Fr),max(Fr),64);
for i = 1:numel(f)
    cr(i) = mean(Cr(Fr > f(i)));
end
subplot(2,2,2), hold on
plot(f,cr*100,'c'), hold off


% illustrate information geometry
%==========================================================================
spm_dir_disentangle(mdp.a);

% Uncomment the following code to look at particular test digits (j)
%----------------------------------------------------------------------
% [O,D]    = spm_MNIST2o(test,j);
% pdp      = mdp;
% pdp.O    = O;
% pdp.D{2} = D;
% pdp.T    = 1;
% 
% % active inference
% %----------------------------------------------------------------------
% OPTIONS.A = 0;                             % suppress explicit action
% OPTIONS.B = 0;                             % suppress explicit action
% OPTIONS.N = 0;                             % suppress neuronal responses
% OPTIONS.O = 1;                             % probabilistic outcomes
% OPTIONS.P = 1;                             % inference graphics
% pdp = spm_MDP_VB_XXX(pdp,OPTIONS);
% subplot(6,4,9),  imagesc(spm_o2MNIST(O)),     axis image off
% subplot(6,4,10), imagesc(spm_o2MNIST(pdp.Y)), axis image off

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subroutines
%==========================================================================

function mdp = spm_MNIST_learn(training,NS,NT)
% FORMAT mdp = spm_MNIST_learn(training,NS,T)
% training - training data
% NS       - number of styles
% NT       - number of training data
% mdp      - model struct
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimagin

% structure_learning
%--------------------------------------------------------------------------
spm_figure('GetWin','Structure_learning'); clf;
spm_figure('GetWin','Disentangled'); clf;

OPTIONS.G  = 0;    % suppress graphics
OPTIONS.NF = 1;    % maxmium number of factors
OPTIONS.NS = NS;   % maxmium number of states
OPTIONS.NU = 1;    % maxmium number of paths
OPTIONS.UB = 1;    % upper bound prior
OPTIONS.UG = 1;    % expected free energy prior

Ng    = numel(spm_MNIST2o(training,1));
for g = 1:Ng
    a{g} = zeros(2,NS,10);
end
Nd    = zeros(10000,10);
for n = 1:10
  nstr{n} = sprintf('%i', n - 1);
end

% Structure learning
%==========================================================================
mdp.p    = 1/32;                           % initial Dirichlet prior
mdp.zeta = 8;                              % precision of selection priors
for n = 1:10                               % For each digit class

    % exemplars
    %----------------------------------------------------------------------
    spm_figure('GetWin','Structure_learning');

    j     = find(training.labels == (n - 1));
    pdp   = mdp;
    for u = 1:NT

        % structure_learning
        %------------------------------------------------------------------
        O       = spm_MNIST2o(training,j(u));
        pdp     = spm_MDP_structure_learning(pdp,{O},OPTIONS);

        % graphics
        %------------------------------------------------------------------
        Nd(u,n) = size(pdp.a{1},2);
        [nn,jj] = sort(sum(pdp.a{1}),'descend');
        i       = 1:numel(nn);
        Nn(i,n) = nn;

        % mutual information of likelihood
        %------------------------------------------------------------------
        MI(u,n) = spm_MDP_MI(pdp.a);

        subplot(2,2,3)
        plot(MI(1:u,:)), xlabel('number of exemplars'), ylabel('nats')
        title('Mutual information'), axis square

        subplot(2,2,1)
        plot(Nd(1:u,:)), xlabel('number of exemplars'), ylabel('number of styles')
        title('Style learning'), axis square

        subplot(2,2,2)
        plot(Nn), xlabel('styles'), ylabel('occurences')
        title('Frequency'), axis square
        if u > 1, legend(nstr(1:n)), end, drawnow

    end
    
    % Illustrate learned styles
    %----------------------------------------------------------------------
    spm_figure('GetWin','Disentangled');
    j     = 1:numel(jj);
    for g = 1:numel(pdp.a)
        a{g}(:,j,n) = pdp.a{g}(:,jj);
    end
    spm_show_a(a);

end

% save likelihood tensor in MDP structure
%----------------------------------------------------------------------
mdp.a = a;


return


function [C,G] = spm_MNIST_test(mdp,test,NT)
% FORMAT [C,F] = spm_MNIST_test(mdp,test,NT)
% mdp   - model struct
% test  - test data
% NT    - number of test images
%
% C     - correct classification vector
% F     - ELBO (log marginal likelihood)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% number of test images
%--------------------------------------------------------------------------
if nargin < 3, NT = 10000; end
if NT < 128
    OPTIONS.P = get(gcf,'Number');
else
    OPTIONS.P = 0;
end

% test options
%--------------------------------------------------------------------------
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 0;                             % suppress explicit action
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.O = 1;                             % probabilistic outcomes

mdp.A = mdp.a;                             % likelihood tensors
mdp.T = 1;                                 % no transitions
mdp   = rmfield(mdp,'a');                  % suppress learning
for j = 1:NT
    
    % infer (i.e., classify)
    %----------------------------------------------------------------------
    for i = 1 %%%%:25                      % saccadic serach
        [O,D]  = spm_MNIST2o(test,j);
        pdp    = mdp;
        pdp.O  = O;
        pdp    = spm_MDP_VB_XXX(pdp,OPTIONS);
        Q(:,i) = pdp.X{2};
        F(i)   = pdp.F;
    end

    % Bayesian model averaging
    %----------------------------------------------------------------------
    Q     = Q*spm_softmax(F(:));
    [m,d] = max(D);
    [m,p] = max(Q);

    % graphics
    %----------------------------------------------------------------------
    if OPTIONS.P
        subplot(3,4,5)
        imagesc(spm_o2MNIST(O)), axis image
        title(sprintf('Content: %i',d - 1))

        subplot(3,4,6)
        imagesc(spm_o2MNIST(pdp.Y)), axis image
        title(sprintf('Prediction: %i',p - 1))
        drawnow
    end

    % classification accuracy (and ELBO)
    %----------------------------------------------------------------------
    C(j) = d == p;
    G(j) = F*spm_softmax(F(:));
    clc, disp(100*mean(C)),disp(j)

end

return


function spm_show_a(a)
% display learned images in likelihood tensors
% FORMAT spm_show_a(a)
% a{g} - likelihood tensors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% illustrate images
%--------------------------------------------------------------------------
for n = 1:size(a{1},3)
    for u = 1:size(a{1},2)
        for g = 1:numel(a)
            Y{g}   = spm_dir_norm(a{g}(:,u,n));      % normalised
            % Y{g} = a{g}(:,u,n);                    % unnormalised
        end
        YY{1,u}  = spm_o2MNIST(Y);
    end
    subplot(10,1,n)
    if numel(YY) > 31
        imagesc(spm_cat(reshape(YY(1:32),2,16)))
    else
        imagesc(spm_cat(YY))
    end
    axis image, axis off
    drawnow
end

return


function [O,D] = spm_MNIST2o(mnist,id,d)
% convert image into a probabilistic output
% FORMAT [O,D] = spm_MNIST2o(mnist,id,[d])
% mnist - struct
% id    - content
% d     - location (1 to 25) [optional for mnist.mat images]
%
% O     - cell array of outcome probabilities
% D     - label (prior)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% create array of outcomes (one modality for each pixel)
%--------------------------------------------------------------------------
if isfield(mnist,'I')
    I = mnist.I(:,id);
    D = mnist.D(:,id);
    o = [I, 1 - I];
    O = num2cell(o,2);
    return
end


% If we are dealing with raw MNIST data, pre-process on the fly
%==========================================================================
h    = mnist.height;
w    = mnist.width;
n    = h*w;
I    = mnist.images(:,:,id);

% move image around, if requested
%--------------------------------------------------------------------------
if nargin > 2
    c    = zeros(5,5);
    c(d) = 1;
    I    = conv2(I,c,'same');
end

% smooth and histogram equalisation
%--------------------------------------------------------------------------
I    = spm_conv(I,2,2);                       % convolve
[m,i] = sort(I(:));                           % histogram equalisation
I(i) = exp(-((1:n) - n).^2/(32*n));
for i = 1:h
    for j = 1:w
        o           = I(i,j);
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
% Converts a probabilistic output into an image
% FORMAT [I] = spm_o2MNIST(Y)
% Y      - posterior prediction
% I      - image
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% create array from (output) structures
%--------------------------------------------------------------------------
for i = 1:numel(Y)
    D(i) = Y{i}(1);
end

% assign pixel values
%--------------------------------------------------------------------------
[n,i] = sort(D(:));
n     = numel(D);
D(i)  = exp(-((1:n) - n).^2/(24*n));

% re-arrange into image
%--------------------------------------------------------------------------
if numel(D) ~= (28^2)
    load spm_MNIST pixels;
    I         = zeros(28,28);
    I(pixels) = D;
else
    I         = reshape(D,28,28);
end

return


function [training,test,pixels] = spm_MNIST_prepare
% Prepares and saves spm_MNIST.mat
% FORMAT [training,test,pixels] = spm_MNIST_prepare
% training - struct
% test     - struct
% pixels   - vector of indices corresponding to  observation  modalities
%
% Preprocessing is implemented by a subroutine that says the requisite
% pre-processed images in a MATLAB file. This preprocessing entails
% smoothing with gaussian convolution kernel with a full width at half
% maximum of two pixels. This is followed by a simple histogram
% equalisation so that the distribution of pixel values are the same for
% each image. Finally, the sum of the accumulated counts is used to select
% the 512 most informative pixels that constitute the 512 factors or
% observation modalities. The format of these inputs is probabilistic. In
% other words, each pixel is a categorical distribution encoding the
% probability of this pixel being black or white.
%
% The requisite structures (and pixel indices) are saved in spm_MNIST (in
% the current working directory)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% This subroutine expects to find mnist.mat in its working directory. The
% requisite MATLAB file can be downloaded from the following website
%--------------------------------------------------------------------------
% https://lucidar.me/en/matlab/load-mnist-database-of-handwritten-digits-in-matlab/

% load training and test data
%--------------------------------------------------------------------------
load('mnist.mat');

% create array of outcomes (one modality for each pixel)
%--------------------------------------------------------------------------
m     = mean(training.images,3);                  % average intensity
[m,d] = sort(m(:),'descend');                     % sort
N     = 512;                                      % number of modalities
s     = 2;                                        % smoothing in pixels  
S     = 32;                                       % histogram width
d     = d(1:N);                                   % indices of pixels                                 
n     = numel(training.images(:,:,1));            % number of pixels

% smooth and histogram equalisation (training data)
%--------------------------------------------------------------------------
for id = 1:20000

    I     = training.images(:,:,id);
    I     = spm_conv(I,s,s);                      % convolve
    [m,i] = sort(I(:));                           % histogram equalisation
    I(i)  = exp(-((1:n) - n).^2/(S*n));

    % return a cell for this image
    %----------------------------------------------------------------------
    training.I(:,id) = I(d);
    training.D(:,id) = sparse(training.labels(id) + 1,1,1,10,1);

end

% smooth and histogram equalisation  (test data)
%--------------------------------------------------------------------------
for id = 1:10000

    I     = test.images(:,:,id);
    I     = spm_conv(I,1,1);                      % convolve
    [m,i] = sort(I(:));                           % histogram equalisation
    I(i)  = exp(-((1:n) - n).^2/(S*n));

    % return a cell for this image
    %----------------------------------------------------------------------
    test.I(:,id) = I(d);
    test.D(:,id) = sparse(test.labels(id) + 1,1,1,10,1);

end

% prepare structures for saving
%--------------------------------------------------------------------------
test.count     = 10000;
training.count = 20000;
test           = rmfield(test,'images');
training       = rmfield(training,'images');
pixels         = d;

% save in spm_MNIST
%--------------------------------------------------------------------------
save spm_MNIST training test pixels

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below are some additional subroutines that explore alternative
% formulations of active learning and selection
%--------------------------------------------------------------------------



function mdp = spm_initialise(training,NS)
% FORMAT mdp = spm_initialise(training,NS)
% training - training data
% NS       - number of styles
%
% mdp      - model struct
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging


% initial Dirichlet prior
%--------------------------------------------------------------------------
mdp.p = 1/32;

% initialise learnable language mapping for all agents
%--------------------------------------------------------------------------
Ns    = [NS,10];
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
for n = 1:Ns(2)
    j     = find(training.labels == (n - 1));
    j     = j(randperm(numel(j)));
    for u = 1:Ns(1)
        O     = spm_MNIST2o(training,j(u));
        for g = 1:Ng
            mdp.a{g}(:,u,n) = O{g} + mdp.p;
        end
    end
end

return




function mdp = spm_initialise_resolve(training,NS)
% FORMAT mdp = spm_initialise_resolve(training,NS)
% training - training data
% NS       - number of styles
%
% mdp      - model struct
%
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% resolve 
%--------------------------------------------------------------------------
spm_figure('GetWin','resolve'); clf

Ng    = numel(spm_MNIST2o(training,1));
for g = 1:Ng
    A{g} = zeros(2,NS,10);
end
Nd    = zeros(10000,10);
for n = 1:10
  nstr{n} = sprintf('%i', n - 1);
end

% initial Dirichlet prior
%--------------------------------------------------------------------------
mdp.p = 1/32;
mdp.T = 1;
for n = 1:10

    % exemplars
    %----------------------------------------------------------------------
    for g = 1:Ng
        mdp.a{g} = zeros(2,0);
    end
    j     = find(training.labels == (n - 1));
    for u = 1:numel(j)

        % structure_learning
        %------------------------------------------------------------------
        O     = spm_MNIST2o(training,j(u));
        for g = 1:Ng
            mdp.a{g}(:,end + 1) = O{g} + mdp.p;
        end
        
        mdp = spm_resolve(mdp);

        % graphics
        %------------------------------------------------------------------
        spm_figure('GetWin','resolve');

        Nd(u,n) = size(mdp.a{1},2);
        [nn jj] = sort(sum(mdp.a{1}),'descend');
        i       = 1:numel(nn);
        Nn(i,n) = nn;

        % mutual information of likelihood
        %------------------------------------------------------------------
        MI(u,n) = spm_MDP_MI(mdp.a);

        subplot(2,2,3),   plot(MI(1:u,:))
        xlabel('number of exemplars'), ylabel('nats')
        title('Mutual information'), axis square

        subplot(2,2,1),   plot(Nd(1:u,:))
        xlabel('number of exemplars'), ylabel('number of styles')
        title('Style learning'), axis square

        subplot(2,2,2),   plot(Nn)
        xlabel('styles'), ylabel('occurences')
        title('Frequency'), axis square, drawnow
        if u > 1, legend(nstr(1:n)), legend('Location','SouthEast'), end

        if size(mdp.a{1},2) == NS, break, end

    end
    
    % disentangle
    %----------------------------------------------------------------------
    spm_figure('GetWin','Disentangled');
    j     = 1:numel(jj);
    for g = 1:numel(mdp.a)
        A{g}(:,j,n) = mdp.a{g}(:,jj);
    end
    spm_show_a(mdp.a);

end

% save in mdp
%----------------------------------------------------------------------
mdp.a = A;

return


function mdp = spm_train(mdp,training,NT)
% FORMAT mdp = spm_train(mdp,training,NT)
% mdp      - model struct
% training - train data
% NT       - number of training samples
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% illustrate behavioural responses
%--------------------------------------------------------------------------
spm_figure('GetWin','train'); clf;

% default priors
%--------------------------------------------------------------------------
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 0;                             % suppress explicit action
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.O = 1;                             % probabilistic outcomes
OPTIONS.G = 1;                             % suppress graphics
OPTIONS.P = 0;                             % suppress graphics

[Nf,Ns,Nu]   = spm_MDP_size(mdp);
for f = 1:Nf
    mdp.D{f} = ones(Ns(f),1);              % initial state
    mdp.E{f} = ones(Nu(f),1);              % initial control
end

% initial prior (Uncomment first line to invoke Bayesian model reduction)
%--------------------------------------------------------------------------
% OPTIONS.BMR.T = 0;                       % Occam's window [default: 0]
try mdp.p;    catch, mdp.p    = 1/32; end
try mdp.beta; catch, mdp.beta = 1;    end

%  train
%--------------------------------------------------------------------------
mdp.T = 1;
mdp.U = zeros(1,Nf);
for j = 1:NT

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    [O,D]    = spm_MNIST2o(training,training.count - j);
    pdp      = mdp;
    pdp.O    = O;
    pdp.D{2} = D;

    % active inference and learning
    %----------------------------------------------------------------------
    pdp = spm_MDP_VB_XXX(pdp,OPTIONS);

    % update Dirichlet parameters [with/out Bayesian model reduction]
    %----------------------------------------------------------------------
    mdp = spm_MDP_VB_update(mdp,pdp,OPTIONS);

    if OPTIONS.G

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(3,4,5)
        imagesc(spm_o2MNIST(O))
        [m,d] = max(D); axis square
        title(sprintf('Content: %i',d - 1))

        subplot(3,4,6)
        imagesc(spm_o2MNIST(pdp.Y))
        [m,d] = max(pdp.X{2}); axis square
        title(sprintf('Content: %i',d - 1))

        % mutual information of likelihood mapping
        %------------------------------------------------------------------
        I     = spm_zeros(spm_o2MNIST(O));
        for g = 1:numel(mdp.a)
            I(g) = I(g) + spm_MDP_MI(mdp.a{g});
        end
        MI(j) = sum(I,'all');
        FI(j) = sum(pdp.F,'all');

        subplot(3,2,2), plot(MI)
        xlabel('number of exemplars'), ylabel('nats')
        title('Mutual information'), axis square
        subplot(3,2,4), plot(FI)
        xlabel('number of exemplars'), ylabel('nats')
        title('ELBO'), axis square, drawnow
    end

end

return


