function mdp = DEM_MNIST_XXX
% Active sampling of MNIST digits
%__________________________________________________________________________
%
% This routine demonstrates smart or sparse sampling of outcomes to
% accumulate evidence for posterior beliefs in an optimal way. Effectively,
% it scores the expected free energy (i.e., expected information gain) of
% different policies. Crucially, these policies reflect covert or internal
% action; namely, subscribing to or querying different subsets of output
% modalities. This can be regarded as one aspect of proactive message
% passing; in the sense, it rests upon expected free energy to send
% requests to nodes that evaluate outcome probabilities. The saving here is
% in terms of compute time and message passing. Interestingly, this can
% lead to characteristic overconfidence in classification. In other words,
% by just sampling a limited subset of outcome modalities, the routine can
% find a plausible explanation very quickly; however, not necessarily the
% right explanation.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $


%% set up and preliminaries
%==========================================================================

% Load (pre-processed) MNIST dataset (uncomment to save mnist_mdp)
%--------------------------------------------------------------------------
% load mnist
% NS  = 8;                                   % number of styles
% mdp = spm_initialise_resolve(training,NS);
% save mnist_mdp mdp
%--------------------------------------------------------------------------
cd 'C:\Users\karl\Dropbox\matlab'
load mnist
load mnist_mdp

% Indices of outcome nodes
%==========================================================================
d     = 9;
id    = {};
for i = 0:d:(28 - d)
    for j = 0:d:(28 - d)
        n  = sparse((1:d) + i,1,1,28,1)*sparse((1:d) + j,1,1,28,1)';
        id{end + 1} = find(n(:));
    end
end

% illustrate Sparse sampling of outcomes
%--------------------------------------------------------------------------
spm_figure('GetWin','train'); clf;

% default priors
%--------------------------------------------------------------------------
% test options
%--------------------------------------------------------------------------
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 0;                             % suppress replay
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.G = 1;                             % suppress neuronal responses

mdp.a = spm_vecfun(mdp.a,@times,128);      % likelihood tensors

[Nf,Ns,Nu]   = spm_MDP_size(mdp);
for f = 1:Nf
    mdp.D{f} = ones(Ns(f),1);              % initial state
    mdp.E{f} = ones(Nu(f),1);              % initial control
end

mdp.beta  = 0;                             % active selection
mdp.eta   = exp(16);                       % active forgetting

%  train
%--------------------------------------------------------------------------
mdp.id.ig = 1;
mdp.id.g  = id;
mdp.T = 6;
mdp.U = zeros(1,Nf);
for j = 1:4

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    [O,n]    = spm_MNIST2o(training,training.count - j);
    pdp      = mdp;
    pdp.O    = repmat(O,1,mdp.T);
    % pdp.D{2} = D;                       % Cmt., precise veridical priors

    % active inference and learning
    %----------------------------------------------------------------------
    pdp = spm_MDP_VB_XXX(pdp,OPTIONS);

    % Graphics
    %======================================================================
    if OPTIONS.G

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(4,3,10)
        I     = spm_o2MNIST(O);
        imagesc(I)
        [m,d] = max(n); axis square
        title(sprintf('Content: %i',d - 1))

        for t = 1:min(mdp.T,6)

            % Observations
            %--------------------------------------------------------------
            subplot(8,4,4*(t - 1) + 1)
            i     = NaN(size(I));
            ig    = pdp.id.g{pdp.id.ig(t)};
            i(ig) = I(ig);
            imagesc(i)
            axis square, title(sprintf('Seen: %i',d - 1))

            % Predictions
            %--------------------------------------------------------------
            subplot(8,4,4*(t - 1) + 2)
            imagesc(spm_o2MNIST(pdp.Y(:,t)))
            [m,s] = max(pdp.X{2}(:,end)); axis square
            title(sprintf('Prediction: %i',s - 1))

            % Posteriors
            %--------------------------------------------------------------
            subplot(8,4,4*(t - 1) + 3)
            bar(pdp.X{1}(:,t))
            title('Posteriors')

            % Posteriors
            %--------------------------------------------------------------
            subplot(8,4,4*(t - 1) + 4)
            bar(pdp.X{2}(:,t))
            title('Posteriors')

            % Dirichlet counts
            %--------------------------------------------------------------
            ax     = 0;
            for g  = 1:numel(pdp.a)
                ax = ax + sum(pdp.a{g},1);
            end
            subplot(4,3,11),
            imagesc(reshape(squeeze(ax),Ns))
            title('Dirichlet counts')

            % Dirichlet counts
            %--------------------------------------------------------------
            xa     = I;
            for g  = 1:numel(pdp.a)
                xa(g) = sum(pdp.a{g},'all');
            end
            subplot(4,3,12),
            imagesc(xa)
            title('Dirichlet counts')

            drawnow

        end

        % update Dirichlet parameters [with/out Bayesian model reduction]
        %----------------------------------------------------------------------
        mdp = spm_MDP_VB_update(mdp,pdp,OPTIONS);

    end

end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subroutines
%==========================================================================

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
    o = [I,1 - I]';
    O = num2cell(o,1);
    O = O(:);
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
I     = spm_conv(I,2,2);                       % convolve
[m,i] = sort(I(:));                            % histogram equalisation
I(i)  = exp(-((1:n) - n).^2/(32*n));
for i = 1:h
    for j = 1:w
        o           = I(i,j);
        O{i,j}(1,1) = o;
        O{i,j}(2,1) = 1 - o;
    end
end

% return a cell array for this image
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
cd 'C:\Users\karl\Dropbox\matlab'
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
Ntrain = 30000;
for id = 1:Ntrain

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
Ntest  = 10000;
for id = 1:Ntest

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
training.count  = Ntrain;
training.labels = training.labels(1:Ntrain);
training        = rmfield(training,'images');

test.count      = Ntest;
test.labels     = test.labels(1:Ntest);
test            = rmfield(test,'images');

pixels          = d;

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
        [nn,jj] = sort(sum(mdp.a{1}),'descend');
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
        title('Frequency'), axis square
        if u > 1, legend(nstr(1:n)), legend('Location','SouthEast'), end
        drawnow
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
OPTIONS.G = 1;                             % suppress graphics
OPTIONS.P = 0;                             % suppress graphics

[Nf,Ns,Nu]   = spm_MDP_size(mdp);
for f = 1:Nf
    mdp.D{f} = ones(Ns(f),1);              % initial state
    mdp.E{f} = ones(Nu(f),1);              % initial control
end

% Uncomment first line to invoke Bayesian model reduction
%--------------------------------------------------------------------------
% OPTIONS.BMR.T = 0;                       % Occam's window [default: 0]


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


%% NOTES

%==========================================================================
n = 8;
Y = randn(n,n);
subplot(4,4,1), imagesc(Y), axis image
subplot(4,4,2), imagesc(spm_speye(n,n,1,2,2)'*Y), axis image
subplot(4,4,3), imagesc(spm_speye(n,n,1,2,2)'*Y*spm_speye(n,n,-1,2)), axis image
subplot(4,4,4), imagesc(reshape(kron(spm_speye(n,n,-1,2),spm_speye(n,n,1,2,2))'*Y(:),n,n)), axis image





