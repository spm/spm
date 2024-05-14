function MDP = DEM_video_compression
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
rng(1)

% Load a video clip
%==========================================================================
vid   = VideoReader('C:\Users\Karl\Dropbox\matlab\robin.mp4');

% Read into working memory and show the movie
%==========================================================================
spm_figure('GetWin','video'); clf

% Crop and resize to 128 x 128 pixels
%--------------------------------------------------------------------------
k     = fix(linspace(1,vid.NumFrames/2,128));
k     = [k k(1:32)];
for t = 1:numel(k)
    temp       = read(vid,k(t));
    I(:,:,:,t) = imresize(temp(:,1:360,:),[128,128]);
end

subplot(2,2,1)
spm_imshow(I)

% Map from image to discrete state space (c.f., Amortisation) 
%--------------------------------------------------------------------------
RGB.nd    = 16;                    % Diameter of tiles in pixels
RGB.nb    = 16;                    % Number of discrete singular variates 
RGB.mm    = 32;                    % Maximum number of singular modes
RGB.su    = 32;                    % Variance threshold
RGB.R     = 2;                     % temporal resampling
[O,L,RGB] = spm_rgb2O(I,RGB);

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
spm_show_RGB(PDP,RGB,8,1)


% Biomimetic characterisation in terms of evoked responses
%==========================================================================
spm_figure('GetWin','Evoked responses'); clf

% Create partial stimulus (c.f., Moving bar or grating stimuli)
%--------------------------------------------------------------------------
S     = O(:,1:end);
for g = 1:size(S,1)
    for t = 1:size(S,2)
        if ismember(t,[3 4 7])
            S{g,t} = spm_dir_norm(ones(size(O{g,t})));
        end
    end
end

% Invert deep (RG) model
%--------------------------------------------------------------------------
MDP{1}.S  = S;
RDP       = spm_mdp2rdp(MDP,1/256);

RDP.T     = fix(size(S,2)/(2^(Nm - 1)));
RDP.D{1}  = ones(Ns(1),1);
RDP.E{1}  = ones(Nu(1),1);
PDP       = spm_MDP_VB_XXX(RDP);

% Illustrate neuronal responses
%--------------------------------------------------------------------------
spm_show_RGB(PDP,RGB,8,1)

return

% subroutines
%==========================================================================

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
