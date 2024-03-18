function MDP = DEM_compression
% Structure learning from pixels
%__________________________________________________________________________
%
% This routine demonstrates the mapping from pixels to (discrete)
% representations of episodes. It illustrates the use of a certain kind of
% deep generative model based upon reduction and grouping operators (RG)
% found in the renormalisation group. In brief, this allows an efficient
% encoding (and generation) of high dimensional content by maintaining
% conditionally independent partitions of hidden states in a recursive
% fashion. This example takes the view of efficient inference as,
% effectively, compression and efficient coding of videos. However, the
% focus of this compressed encoding is on the trajectory or path as visual
% objects move.
% 
% A key architectural aspect of these deep structures is an appeal to the
% renormalisation group; in the sense that there is a recursive application
% of grouping and reduction operators, as we move from one hierarchical
% level to the next. These operate over both space and time: in the sense
% that groups of pixels are jointly inferred over a short period of time;
% such that the state any given level of the model generates the
% instantaneous state and path (i.e., velocity) of a group of latent states
% at the lower level. In these examples, the group operator effectively
% tiles a lattice of latent states (and instantaneous paths) into little
% contiguous squares. This necessary induces a separation of temporal
% scales such as the states at the highest level see all the pixels at the
% lowest level and, effectively, encode a trajectory over to the 2^n
% timesteps, where n is the depth of the model.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Copyright free videos obtained from the following website:
% https://pixabay.com/videos/search/birds/

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%--------------------------------------------------------------------------
close all
clc, rng(1)

% Load a video clip
%==========================================================================
vid   = VideoReader('C:\Users\Karl\Dropbox\matlab\dove.mp4');

% Read into working memory and show the movie
%==========================================================================
spm_figure('GetWin','video'); clf

% Crop and resize to 128 x128 pixels
%--------------------------------------------------------------------------
i     = (1:256) - 128 + round(vid.Width/2);
j     = (1:256) - 128 + round(vid.Height/2);
k     = 1:3:vid.NumFrames;
for t = 1:numel(k)
    temp       = read(vid,k(t));
    I(:,:,:,t) = imresize(temp(j,i,:),[128,128]);
end

subplot(2,2,1)
spm_imshow(I)

% Map from image to discrete state space (c.f., Amortisation) 
%--------------------------------------------------------------------------
[O,L,RGB] = spm_rgb2O(I);

% And show the images generated from a discrete representation
%--------------------------------------------------------------------------
subplot(2,2,2)
T     = size(O,2);
for t = 1:T
    I = spm_O2rgb(O(:,t),RGB);
    spm_imshow(I), drawnow
end

% Use the ensuing sequence for (RG) structure learning
%--------------------------------------------------------------------------
MDP   = spm_MB_structure_learning(O,L);

% Show model structure (at the deepest level)
%--------------------------------------------------------------------------
Nm    = numel(MDP);
spm_figure('GetWin',sprintf('Paramters: level %i',Nm)); clf
spm_MDP_params(MDP{Nm})


% Generate movie by sampling from the resulting deep generative model
%==========================================================================
spm_figure('GetWin','Generative AI'); clf

% sample outcomes
%--------------------------------------------------------------------------
pdp   = MDP{end};
pdp.T = 128/(2^(Nm - 1));
pdp.s = 1;
pdp.u = 1;
pdp   = spm_MDP_VB_XXX(pdp);
Q     = cell(Nm,1);
Q{Nm} = pdp.O;
for n = Nm:-1:2

    % set empirical priors over states and paths
    %----------------------------------------------------------------------
    for t = 1:size(Q{n},2)
        pdp   = MDP{n - 1};
        pdp.T = 2;
        for g = 1:numel(pdp.id.D)
            pdp.D{g} = Q{n}{pdp.id.D{g},t};
            pdp.E{g} = Q{n}{pdp.id.E{g},t};
        end
        pdp      = spm_MDP_VB_XXX(pdp);
        Q{n - 1} = [Q{n - 1} pdp.O];
    end

end

% show generated image
%--------------------------------------------------------------------------
spm_show_outcomes(Q,RGB)

% Generate play from recursive generative model
%==========================================================================
spm_figure('GetWin','Active Inference'); clf

% Create deep recursive model
%--------------------------------------------------------------------------
RDP       = spm_mdp2rdp(MDP);
[~,Ns,Nu] = spm_MDP_size(RDP);

RDP.T     = 128/(2^(Nm - 1));
RDP.D{1}  = sparse(1,1,1,Ns(1),1);
RDP.E{1}  = sparse(1,1,1,Nu(1),1);
PDP       = spm_MDP_VB_XXX(RDP);

% Illustrate the model in generative mode
%--------------------------------------------------------------------------
spm_show_MDP(PDP,RGB)


% Biomimetic characterisation in terms of evoked responses
%==========================================================================
spm_figure('GetWin','Evoked responses'); clf

% Create partial stimulus (c.f., Moving bar or grating stimuli)
%--------------------------------------------------------------------------
S     = O(:,1:end);
for g = 1:size(S,1)
    for t = 1:size(S,2)
        if ~ismember(g,MDP{1}.RG{3})
            S{g,t} = spm_dir_norm(ones(size(O{g,t})));
        end
    end
end

%%% S     = [O(:,1:16) O(:,1:end)];

% Invert deep (RG) model
%--------------------------------------------------------------------------
MDP{1}.S  = S;
RDP       = spm_mdp2rdp(MDP,1/32);

RDP.T     = fix(size(S,2)/(2^(Nm - 1)));
RDP.D{1}  = ones(Ns(1),1);
RDP.E{1}  = ones(Nu(1),1);
PDP       = spm_MDP_VB_XXX(RDP);

% Illustrate neuronal responses
%--------------------------------------------------------------------------
spm_show_MDP(PDP,RGB)

return

% subroutines
%==========================================================================


function spm_imshow(I)
% rrendering of an RGB image sequence
% FORMAT spm_imshow(I)
%
% Converts discrete (probabilistic) image into an RGB format for image
% display
%__________________________________________________________________________

% image
%--------------------------------------------------------------------------
for t = 1:size(I,4)
   imshow(uint8(I(:,:,:,t)))
   drawnow
end

return


function [O,L,RGB] = spm_rgb2O(I)
% RGB image format to O
% FORMAT [O,L,RGB] = spm_rgb2O(I);
% 
% I    - RGB image (pixels x pixels x 3 colours): uint8
%
% O    - Cell array of probabilistic outputs for pixel partition
% L    - Location of each probabilistic output modality
% RBG  - structure containing image size, pixel locations, singular vectors
%        and bins for singular variants
%
% Converts an RGB image into discrete probabilistic format using singular
% value decomposition of a partition into local pixels.
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
nd    = 32;               % Diameter of tiles in pixels 
nb    = 8 + 1;            % Number of discrete singular variatesaal levels  
mm    = 32;               % Maximum number of singular modes
su    = 1/16;             % Singular value (normalised) threshold
R     = 2;                % temprial resampling

% locations for spatial mapping
%--------------------------------------------------------------------------
T     = R*fix(size(I,4)/R);
I     = I(:,:,:,1:T);
I     = permute(I,[3,4,1,2]);
N     = size(I,[1,2,3,4]);
I     = reshape(I,3*R,[],N(3),N(4));

I     = permute(I,[2,1,3,4]);
N     = size(I,[2,3,4]);
T     = size(I,1);
for n = 1:prod(N)
    ind    = spm_index(N,n);
    L(n,:) = [ind(2),ind(3)];
end

[G,M] = spm_tile(L,nd);
RGB.M = M;
RGB.N = N;
RGB.G = G;

O     = {1,T};
L     = [1,1];
o     = 1;
for g = 1:numel(G)

    % singular value decomposition for this group
    %----------------------------------------------------------------------
    [u,s,v]  = spm_svd(double(I(:,G{g})),su);
    Nm       = min(length(s),mm);
    if Nm
        RGB.V{g} = v(:,1:Nm)*s(1:Nm,1:Nm);
    end

    % generate (probability over discrete) outcomes
    %----------------------------------------------------------------------
    for m = 1:Nm
        
        % dicretise singular variates
        %------------------------------------------------------------------
        d     = max(abs(u(:,m)));
        a     = linspace(-d,d,nb);
        for t = 1:T
            [~,U]  = min(abs(u(t,m) - a));
            O{o,t} = sparse(U,1,1,nb,1);
        end

        % record locations and group for this outcome
        %------------------------------------------------------------------
        L(o,:)     = RGB.M(g,:);
        RGB.O(o)   = g;
        RGB.A{o}   = a;
        o          = o + 1;
    end
end

return


function I = spm_O2rgb(O,RGB)
% O to RGB image format
% FORMAT I = spm_O2rgb(O,RGB)
%
% Converts discrete (probabilistic) image into an RGB format for image
% display
%__________________________________________________________________________

% Generate image
%--------------------------------------------------------------------------
I     = zeros(RGB.N);
for g = 1:numel(RGB.G)

    % recover singular vectors
    %----------------------------------------------------------------------
    o     = find(ismember(RGB.O,g));
    Nm    = numel(o);
    u     = zeros(Nm,1);
    for m = 1:Nm
        u(m) = RGB.A{o(m)}*O{o(m)};
    end
    
    % place in image, unless zeros
    %----------------------------------------------------------------------
    if numel(u)
        I(RGB.G{g}) = RGB.V{g}*u;
    end

end

% permute colour to trailing dimension
%--------------------------------------------------------------------------
I = reshape(I,3,[],RGB.N(2),RGB.N(3));
I = uint8(permute(I,[3,4,1,2]));

return


function Q = spm_get_O(RDP)
% hierarchical outcomes from a recursive MDP
% FORMAT Q = spm_get_O(RDP)
%
% Extracts hierarchal outcomes from a recursive MDP that have been
% accumulated during recursive inversion
%__________________________________________________________________________

% recover outcomes
%--------------------------------------------------------------------------
try
    Q = RDP.Q.O;
catch
    Q = {};
end
Q{end + 1} = RDP.O;

return


function spm_show_MDP(MDP,RGB)
% Graphical illustration of active inference
% FORMAT spm_show_MDP(MDP,RGB)
% MDP   - Hierarchical (RG) generative model (inverted)
% RGB   - Structure for mapping discrete outcomes to RGB pixels
%
%--------------------------------------------------------------------------
% This auxiliary routine plots the hierarchal outcomes from a deep
% generative model. The separation of temporal scales is illustrated by
% showing the outcomes as a movie saved as the user data in the current
% figure. This routine separates predicted states and paths and shows the
% first few (predictions and observations) in pixel space.
%__________________________________________________________________________

% Get observations and predictions
%--------------------------------------------------------------------------
O     = spm_get_O(MDP);                      % observations
try
    Y = MDP.Q.Y;
catch
    Y = {};
end
Y{end + 1} = MDP.Y;                          % predictions

% Number of timesteps and levels
%==========================================================================
T     = size(O{1},2);
Nn    = numel(O);

% preditive posteriors over hidden states and paths: raster format
%--------------------------------------------------------------------------
mdp = MDP;
for n = 1:Nn
    
    L   = Nn - n + 1;                        % hierarchical level
    id  = spm_cat(mdp.id.A);                 % parents 
    is  = ismember(rem(id,2),1);             % odd: states
    iu  = ismember(rem(id,2),0);             % even: paths

    % predicted states (and entropy)
    %----------------------------------------------------------------------
    subplot(Nn + 2,2,(n - 1)*2 + 1)

    R   = spm_cat(Y{L}(is,:));
    hold off, image((1 - spm_cat(R))*64), axis xy
    title(sprintf('Predictive posteriors at level %i',L),'FontSize',12)
    
    R   = spm_dir_norm(R);
    H   = sum(-R.*spm_log(R));
    KL  = sum(R(:,2:end).*(spm_log(R(:,2:end)) - spm_log(R(:,1:end - 1))));
    hold on, plot(KL,'r')
    hold on, plot(10*H, 'g')

    % predicted paths (and entropy)
    %----------------------------------------------------------------------
    subplot(Nn + 2,2,(n - 1)*2 + 2)
    R   = spm_cat(Y{L}(iu,:));
    hold off, image((1 - spm_cat(R))*64), axis xy
    title(sprintf('Predictive posteriors at level %i',L),'FontSize',12)

    R   = spm_dir_norm(R);
    H   = sum(-R.*spm_log(R));
    KL  = sum(R(:,2:end).*(spm_log(R(:,2:end)) - spm_log(R(:,1:end - 1))));
    hold on, plot(KL,'r')
    hold on, plot(10*H, 'g')

    % next level
    %----------------------------------------------------------------------
    if n < Nn
        mdp = mdp.MDP;
    end
end

% At the lowest (pixel) level
%--------------------------------------------------------------------------
Nt    = 4;                             % Number of timesteps to show images 
for t = 1:T

    % observed outcomes
    %----------------------------------------------------------------------
    X  = spm_O2rgb(O{1}(:,t),RGB);

    % predicted outcomes
    %----------------------------------------------------------------------
    P  = spm_O2rgb(Y{1}(:,t),RGB);

    % plot first level outcomes as a coloured image
    %----------------------------------------------------------------------
    subplot(Nn + 2,Nt,Nt*Nn + min(Nt,t))
    spm_imshow(P)
    title(sprintf('Predictons at time %i',t),'FontSize',12)

    % plot first level predictions as a coloured image
    %----------------------------------------------------------------------
    subplot(Nn + 2,Nt,Nt*(Nn + 1) + min(Nt,t))
    spm_imshow(X)
    title(sprintf('Stimulus at time %i',t),'FontSize',12)

    % save movie
    %----------------------------------------------------------------------
    drawnow
    I(t) = getframe(gcf);

end

% Place movie in graphic subject
%--------------------------------------------------------------------------
set(gcf,'Userdata',[])
set(gcf,'Userdata',{I,8})
set(gcf,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

return


function spm_show_outcomes(O,RGB,Y)
% Plots the inferred sequence of moves
% FORMAT spm_show_outcomes(O,RGB)
% MDP   - Hierarchical (RG) generative model (inverted)
% RGB   - Structure for mapping discrete outcomes to RGB pixel
% Y     - Cell of image predictions
%--------------------------------------------------------------------------
% This auxiliary routine plots the hierarchal outcomes from a deep
% generative model. The separation of temporal scales is illustrated by
% showing the outcomes as a movie saved as the user data in the current
% figure.
%__________________________________________________________________________


% Get observations and predictions
%==========================================================================
T     = size(O{1},2);
Nn    = numel(O);

% Show predictive posteriors over hierarchical outcomes
%--------------------------------------------------------------------------
for n = 1:Nn
    subplot(Nn + 1,1,n + 1)
    imagesc(1 - spm_cat(O{n})), axis xy
    title(sprintf('Predictive posteriors at level %i',n),'FontSize',12)
end

% And that show final predictions in pixel space
%--------------------------------------------------------------------------
for t = 1:T

    % Predicted and observed (RGB) image
    %----------------------------------------------------------------------
    X  = spm_O2rgb(O{1}(:,t),RGB);

    if nargin > 2

        % predicted outcomes
        %------------------------------------------------------------------
        P  = spm_O2rgb(Y{1}(:,t),RGB);

        % plot first level outcomes as a coloured image
        %------------------------------------------------------------------
        subplot(Nn + 1,2,1)
        spm_imshow(P)
        title(sprintf('Outcomes at time %i',t),'FontSize',12)

        % plot first level predictions as a coloured image
        %------------------------------------------------------------------
        subplot(Nn + 1,2,2)
        spm_imshow(X)
        title(sprintf('Predictions at time %i',t),'FontSize',12)

    else

        % plot first level predictions as a coloured image
        %------------------------------------------------------------------
        subplot(Nn + 1,1,1)
        spm_imshow(X)
        title(sprintf('Predictions at time %i',t),'FontSize',12)

    end

    % save movie
    %----------------------------------------------------------------------
    drawnow
    I(t) = getframe(gcf);

end

% Place movie in graphic subject
%--------------------------------------------------------------------------
set(gcf,'Userdata',[])
set(gcf,'Userdata',{I,8})
set(gcf,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

return

function [G,M] = spm_tile(L,d)
% Grouping into a partition of non-overlapping outcome tiles
% FORMAT [G,M] = spm_tile(L,d)
%--------------------------------------------------------------------------
% L  - locations
% d  - Number of rows per tile [default: 2 or 3] 
%
% G  - Cell array of outcome indices
% M  - and their locations
%
% Effectively, this leverages the conditional independencies that inherit
% from local interactions; of the kind found in metric spaces that preclude
% action at a distance.
%--------------------------------------------------------------------------

% locations
%--------------------------------------------------------------------------
u     = unique(L,'rows','stable');
m     = size(L,1)/size(u,1);
Nr    = numel(unique(L(:,1)));
Nc    = numel(unique(L(:,2)));

% defaults
%--------------------------------------------------------------------------
if nargin < 2

    % use 3 x 3 tiles (or smaller)
    %----------------------------------------------------------------------
    if ~rem(Nr,3), dr = 3; else, dr = 2; end
    if ~rem(Nc,3), dc = 3; else, dc = 2; end

else
    dr = d;
    dc = d;
end

% deal with single row (or column) cases
%--------------------------------------------------------------------------
dr    = min(dr,Nr);
dc    = min(dc,Nc);

% Decimate rows and columns
%--------------------------------------------------------------------------
r     = 0:dr:(Nr - dr);
c     = 0:dc:(Nc - dc);
for i = 1:numel(r)
    for j = 1:numel(c)
        n = sparse((1:dr*m) + r(i)*m,1,1,Nr*m,1)*sparse((1:dc) + c(j),1,1,Nc,1)';
        g{i,j} = find(n(:));
    end
end
G      = g(:);

% mean location of groups
%--------------------------------------------------------------------------
for g = 1:numel(G)
    M(g,:) = mean(L(G{g},:));
end

return