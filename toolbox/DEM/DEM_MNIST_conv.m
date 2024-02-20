function mdp = DEM_MNIST_conv
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
% In addition, it illustrates active sampling under a simple form of weight
% sharing implemented by introducing a state dependent mapping between
% modalities and outcomes. This allows for inference over various
% transformations between the modalities that are generated and the final
% outcomes. In this demonstration, the transformations are restricted to X
% and Y translations of an image. The range of (a fine) transformations can
% be specified by editing the code below.
%
% See also: spm_DCM_conv.m
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $


%% set up and preliminaries
%==========================================================================
rng(1);

% Load  MNIST dataset to create a generative model of number images
%--------------------------------------------------------------------------
cd 'C:\Users\karl\Dropbox\matlab'
load mnist

% Specify the number of styles (ONS) and create a generative model
%--------------------------------------------------------------------------
NS    = 4;
pdp   = spm_initialise(training,NS);

% likelihood tensors
%--------------------------------------------------------------------------
Ng    = numel(pdp.a);
for g = 1:Ng
    mdp.A{g} = pdp.a{g} + 1;               % add ambiguity
end
Nf    = ndims(mdp.A{1}) - 1;

% Partition of outcomes for selective (active) sampling
%==========================================================================
d     = 6;                                 % size of partition (pixels)
id.g  = {};
for i = 0:d:(28 - d)
    for j = 0:d:(28 - d)
        n = sparse((1:d) + i,1,1,28,1)*sparse((1:d) + j,1,1,28,1)';
        id.g{end + 1} = find(n(:));
    end
end
id.ig = 8;                                 % Initial selection

% State dependent mapping from modalities to outcomes
%==========================================================================
d     = 1;                                 % size of translation
id.gg = spm_MDP_conv([28,28],1);

% Supplement hidden states with domain factors
%--------------------------------------------------------------------------
Nff   = ndims(id.gg{1});
id.ff = Nf + (1:Nff);
for f = 1:Nf
    Ns(f) = size(mdp.A{1},f + 1);
    Nu(f) = 1;
end
for f = 1:Nff
    Ns(f + Nf) = size(id.gg{1},f);
    Nu(f + Nf) = 1;
end

% latent states
%--------------------------------------------------------------------------
for f = 1:Nf
    mdp.B{f} = eye(Ns(f),Ns(f));           % transitions (Class and Style)
    mdp.D{f} = ones(Ns(f),1);              % initial state
    mdp.E{f} = ones(Nu(f),1);              % initial control
end

% domain states
%--------------------------------------------------------------------------
for f = Nf + (1:Nff)
    mdp.B{f} = eye(Ns(f),Ns(f));           % transitions
    mdp.D{f} = sparse(2,1,8,Ns(f),1) + 1;  % initial state
    mdp.E{f} = ones(Nu(f),1);              % initial control
end

Nf  = Nf + Nff;                            % number of factors

% Show likelihood mappings to images
%--------------------------------------------------------------------------
spm_figure('GetWin','Numbers and styles');
spm_show_a(mdp.A); drawnow

% default priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Inference'); clf;

% test options
%--------------------------------------------------------------------------
OPTIONS.A = 0;                             % suppress explicit action
OPTIONS.B = 0;                             % suppress replay
OPTIONS.N = 0;                             % suppress neuronal responses
OPTIONS.G = 1;                             % enable graphics
OPTIONS.P = get(gcf,'Number');             % enable graphics

mdp.id    = id;                            % domains and co-domains
mdp.T     = 8;                             % number of samples
mdp.U     = zeros(1,Nf);                   % no controlled paths
mdp.N     = 0;                             % depth of planning

% Illustrate classification (inference) among 10 x NS (translated) images
%==========================================================================
gdp   = mdp;
gdp.T = 1;
for j = 1:4

    % Generate image
    %----------------------------------------------------------------------
    pdp = spm_MDP_VB_XXX(gdp);
    O   = pdp.O(:,1);
    n   = pdp.s(2,1);

    % get observation and place in MDP structure
    %----------------------------------------------------------------------
    pdp   = mdp;
    pdp.O = repmat(O,1,mdp.T);

    % active inference and sampling
    %----------------------------------------------------------------------
    pdp   = spm_MDP_VB_XXX(pdp,OPTIONS);

    % Observed image
    %------------------------------------------------------------------
    subplot(3,4,5)
    I     = spm_o2MNIST(O);
    imagesc(I), axis square
    title(sprintf('Content: %i',n - 1))

    % Final prediction
    %--------------------------------------------------------------
    subplot(3,4,6)
    imagesc(spm_o2MNIST(pdp.Y(:,end)))
    [~,s] = max(pdp.X{2}(:,end)); axis square
    title(sprintf('Prediction: %i',s - 1))


    % Graphics: Belief updating over active sampling
    %======================================================================
    if OPTIONS.G

        % illustrate sparse sampling of outcomes
        %------------------------------------------------------------------
        spm_figure('GetWin','Belief updating'); clf;

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(6,1,6)
        I     = spm_o2MNIST(O);
        imagesc(I), axis square
        title(sprintf('Content: %i',n - 1))

        for t = 1:min(mdp.T,6)

            % Observations
            %--------------------------------------------------------------
            subplot(8,5,5*(t - 1) + 1)
            i     = NaN(size(I));
            try
                ig = pdp.id.g{pdp.id.ig(t)};
            catch
                ig = pdp.id.g{1};
            end
            i(ig) = I(ig);
            imagesc(i)
            axis square, title(sprintf('Seen: %i',n - 1))

            % Predictions
            %--------------------------------------------------------------
            subplot(8,5,5*(t - 1) + 2)
            imagesc(spm_o2MNIST(pdp.Y(:,t)))
            [~,s] = max(pdp.X{2}(:,end)); axis square
            title(sprintf('Prediction: %i',s - 1))

            % Posteriors
            %--------------------------------------------------------------
            subplot(8,5,5*(t - 1) + 3)
            imagesc(pdp.X{1}(:,t)*pdp.X{2}(:,t)')
            title('Posteriors (Class)')

            % Posteriors
            %--------------------------------------------------------------
            if Nf > 3
                subplot(8,5,5*(t - 1) + 4)
                imagesc(pdp.X{3}(:,t)*pdp.X{4}(:,t)')
                title('and translation')
            end

            % Posteriors
            %--------------------------------------------------------------
            if Nf > 4
                subplot(8,5,5*(t - 1) + 5)
                imagesc(pdp.X{5}(:,t)*pdp.X{6}(:,t)')
                title('and rotation')
            end
            drawnow
        end
    end
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subroutines
%==========================================================================

function mdp = spm_initialise(training,NS)
% FORMAT mdp = spm_initialise(training,NS)
% This simply loads in some randomly selected numbers into a likelihood
% mapping and creates an MDP structure training - training data
%
% NS       - number of styles
% mdp      - model struct
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging


% initialise likelihood mapping
%--------------------------------------------------------------------------
Ns    = [NS,10];
Nf    = numel(Ns);
Ng    = numel(spm_MNIST2o(training,1));
for g = 1:Ng
    mdp.a{g} = zeros([2,Ns]);
end
for f = 1:Nf
    mdp.b{f} = eye(Ns(f),Ns(f));
end

% populate likelihood mapping with exemplar numbers
%--------------------------------------------------------------------------
for n = 1:Ns(2)
    j     = find(training.labels == (n - 1));
    j     = j(randperm(numel(j)));
    for u = 1:Ns(1)
        O     = spm_MNIST2o(training,j(u));
        for g = 1:Ng
            mdp.a{g}(:,u,n) = O{g};
        end
    end
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
[~,i] = sort(I(:));                            % histogram equalisation
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
[~,i] = sort(D(:));
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


%% NOTES (Discrete transformations)
%==========================================================================
n = 8;
Y = randn(n,n);
subplot(4,4,1), imagesc(Y), axis image
subplot(4,4,2), imagesc(spm_speye(n,n,1,2,2)'*Y), axis image
subplot(4,4,3), imagesc(spm_speye(n,n,1,2,2)'*Y*spm_speye(n,n,-1,2)), axis image
subplot(4,4,4), imagesc(reshape(kron(spm_speye(n,n,-1,2),spm_speye(n,n,1,2,2))'*Y(:),n,n)), axis image





