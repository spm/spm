function [RDP,RGB,F,C] = DEM_MNIST_compression
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

% Load MNIST
%==========================================================================
cd 'C:\Users\karl\Dropbox\matlab'
load MNIST

% Read into working memory and show the movie
%==========================================================================
spm_figure('GetWin','Images'); clf

% Load images into RGB format suitable for videos
%--------------------------------------------------------------------------
N     = 10000;                                    % Number of training data
T     = 1000;                                     % Number of test data 
s     = 2;                                        % smoothing in pixels  
S     = 24;                                       % histogram width
for t = 1:(N + T)
    temp    = imresize(training.images(:,:,t),[32,32]);
    temp    = spm_conv(temp,s,s);                 % convolve
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
end

% And show an illustrative image
%--------------------------------------------------------------------------
subplot(2,2,1)
spm_imshow(I(:,:,:,1))
title('MNIST image')

% Map from image to discrete state space (c.f., Amortisation) 
%==========================================================================
RGB.nd  = 4;                     % Diameter of tiles in pixels
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

labels = training.labels;
clear training test

% get training images for stucture learning
%--------------------------------------------------------------------------
Ns    = 10;                                % Number of classes
Ne    = 13;
j     = [];
for n = 1:10
    i = find(ismember(labels(1:N),n - 1));
    j = [j;i(1:Ne)];
end


% Use a small number (128) for (RG) structure learning
%--------------------------------------------------------------------------
MDP   = spm_MB_structure_learning(O(:,j),L,1);

% Add class contraints (i.e., number priors) as a final level
%--------------------------------------------------------------------------
No    = size(MDP{end}.B{1},1);             % Number of compressed images 
A     = zeros(No,Ns);
for n = 1:Ns
    i        = ismember(labels(j),n - 1);
    A(i,n)   = 1;
end

mdp.A{1}     = spm_dir_norm(A);            % likelihood (state)
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
    for s = 1:8
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
    try, NDP = NDP.MDP; end
end


% Train: Parametric learning based upon a training dataset
%==========================================================================
% Learn by adding latent states if warranted by active selection to the
% latent factor style, under a supervised or precise prior over the latent
% factor (Digit) class
%--------------------------------------------------------------------------
spm_figure('GetWin','Train'); clf
train = 1:N;
RDP   = spm_MNIST_train(MDP,RGB,O(:,train),D(train));

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
%=========================================================================

function mdp = spm_MNIST_train(MDP,RGB,O,D)
% FORMAT RDP = spm_MNIST_train(MDP,RGB,O,D)
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

% upper (half) bound on mutual information (for plotting)
%--------------------------------------------------------------------------
U     = zeros(1,numel(MDP));
for l = 1:numel(MDP)
    for g = 1:numel(MDP{l}.A)
        U(l) = U(l) + log(min(size(MDP{l}.A{g})));
    end
end

% enable active learning (with minimal forgetting)
%--------------------------------------------------------------------------
N     = size(O,2);
for m = 1:numel(MDP)
    MDP{m}.beta = 512;
    MDP{m}.eta  = 512;
end

% train: with small concentration parameters (1/16)
%--------------------------------------------------------------------------
mdp   = spm_mdp2rdp(MDP,1/16,0,1,FIX);
% mdp.A = MDP{end}.A;
% mdp   = rmfield(mdp,'a');
mdp.T = 1;
for j = 1:N

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    pdp   = spm_RDP_O(mdp,O(:,j),1);
    pdp.D = D(:,j);

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
        [~,d] = max(pdp.D{1}); axis square
        title(sprintf('Content: %i',d - 1))

        subplot(3,4,6)
        spm_imshow(spm_O2rgb(pdp.Q.Y{1},RGB))
        [~,d] = max(pdp.X{1}); axis square
        title(sprintf('Content: %i',d - 1))

    end

    % mutual information of likelihood mapping
    %----------------------------------------------------------------------
    for l = 1:numel(pdp.Q.a)
        I(l,j) = spm_MDP_MI(pdp.Q.a{l});
    end
    try
        I(l + 1,j) = spm_MDP_MI(pdp.a);
    catch
        I(l + 1,j) = spm_MDP_MI(pdp.A);
    end
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

    % disable graphics after first training iteration
    %----------------------------------------------------------------------
    if j > 2
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

% test
%--------------------------------------------------------------------------
for j = 1:size(O,2)

    % get training observation and place in MDP structure
    %----------------------------------------------------------------------
    pdp   = spm_RDP_O(RDP,O(:,[j j]),1);

    % active inference and learning
    %----------------------------------------------------------------------
    pdp   = spm_MDP_VB_XXX(pdp);

    % Bayesian model averaging
    %----------------------------------------------------------------------
    [~,d] = max(D{j});
    [~,p] = max(pdp.X{1}(:,end));

    % classification accuracy (and ELBO)
    %----------------------------------------------------------------------
    C(j)  = d == p;
    G(j)  = pdp.Q.F + sum(pdp.F);
    clc, disp(100*mean(C)), disp(j)

     if ~C(j)

        spm_MDP_VB_trial(pdp)

        % Predicted and observed images
        %------------------------------------------------------------------
        subplot(3,4,5)
        spm_imshow(spm_O2rgb(pdp.Q.O{1}(:,end),RGB))
        title(sprintf('Label: %i',d - 1),'FontSize',12)

        subplot(3,4,6)
        spm_imshow(spm_O2rgb(pdp.Q.Y{1}(:,end),RGB))
        title(sprintf('Class: %i',p - 1),'FontSize',12)
        drawnow

    end

end

function RDP = mdp2rdp(MDP,p,q,T,FIX)
% Converts a cell array of MDPs into a recursive MDP
% FORMAT RDP = mdp2rdp(MDP,p,q,T,FIX)
% MDP{n} - Cell array of MDPs
%  MDP{n}.id.D - cell array of parents of D in supraordinate outcomes
%  MDP{n}.id.E - cell array of parents of E in supraordinate outcomes
%
% p      - likelihood concentration; e.g., p = 1/32; [default: p = 0] 
% q      - prior (transition) decay; e.g., q = 1/4;  [default: q = 0] 
% T      - path lengths [default: T = 2]
% 
% FIX.A = 0 for learning likelihoods
% FIX.B = 0 for learning transitions
%
% RDP    - likelihood and transition tensors for generating this sequence
%  RDP.L - level or depth
%  RDP.T - time steps
%
% This auxiliary routine takes a cell array of hierarchically arranged
% MDP’s and creates a single MDP where each subordinate MDP is a field of
% the superordinate MDP. The vertical dependencies are encoded in the cell
% arrays of parents at each level of the hierarchy, where the outcomes of
% one level are the parents of the initial states and paths of the lower
% level (encoded probabilistically in the vectors D and E of the lower
% level).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% Assume time scaling with a scale doubling
%--------------------------------------------------------------------------
if nargin < 2, p = 0; end
if nargin < 3, q = 0; end
if nargin < 4, T = 2; end

if nargin < 5
    FIX.A = 1;
    FIX.B = 1;
end

Nm    = numel(MDP);
try p = p(1:Nm); catch, p = repmat(p(1),1,Nm); end
try q = q(1:Nm); catch, q = repmat(q(1),1,Nm); end

% Check for single level models
%--------------------------------------------------------------------------
if numel(MDP) < 2
    RDP   = MDP{1};
    RDP.L = 1;
    RDP.T = T;
    return
end

% prior concentration parameters
%--------------------------------------------------------------------------
for n = 1:Nm

    % normalised likelihoods
    %----------------------------------------------------------------------
    for g = 1:numel(MDP{n}.A)
        a           = spm_dir_norm(MDP{n}.A{g});
        MDP{n}.A{g} = spm_dir_norm(a + p(n));
    end

    % normalised transitions
    %----------------------------------------------------------------------
    for f = 1:numel(MDP{n}.B)
        b           = spm_dir_norm(MDP{n}.B{f});
        MDP{n}.B{f} = spm_dir_norm(b + q(n));
    end


    % remove fields
    %----------------------------------------------------------------------
    if FIX.A
        if isfield(MDP{n},'a')
            MDP{n} = rmfield(MDP{n},'a');
        end
    else
        MDP{n}.a = MDP{n}.A;
        MDP{n}   = rmfield(MDP{n},'A');
    end
    if FIX.B
        if isfield(MDP{n},'b')
            MDP{n} = rmfield(MDP{n},'b');
        end
    else
        MDP{n}.b = MDP{n}.B;
        MDP{n}   = rmfield(MDP{n},'B');
    end

end

% Recursively place MDP in MDP.MDP and specify path lengths
%--------------------------------------------------------------------------
SDP   = MDP{1};
SDP.T = T;
SDP.L = 1;
for n = 2:Nm
    RDP     = MDP{n};
    RDP.MDP = SDP;
    SDP     = RDP;
    SDP.T   = T;
    SDP.L   = n;
end
RDP.L = n;

return

return


