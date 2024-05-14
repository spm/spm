function spm_show_WAV(MDP,WAV)
% Graphical illustration of active inference
% FORMAT spm_show_MDP(MDP,WAV)
% MDP   - Hierarchical (RG) generative model (inverted)
% WAV   - Structure for mapping discrete outcomes to WAV pixels
%
%--------------------------------------------------------------------------
% This auxiliary routine plots the hierarchal outcomes (posterior
% predictive distributions) from a deep generative model. The separation of
% temporal scales is illustrated by showing the outcomes as a sonogram.
%
% This version is specialised for sound files, rendering the posterior
% predictions in alp come space as a solar gram. The image graphic is
% equipped with a callback function that will play the corresponding sound
% following an inverse wavelet transform.
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
Nm    = numel(O);

% highest states and transitions
%--------------------------------------------------------------------------
subplot(Nm + 3,2,1)
image((1 - spm_cat(MDP.X))*64)
title(sprintf('Posterior (states) level %i',Nm),'FontSize',12)

subplot(Nm + 3,2,2)
imagesc(1 - MDP.B{1}), axis square
title('Transitions','FontSize',12)

% preditive posteriors over hidden states and paths: raster format
%--------------------------------------------------------------------------
mdp = MDP;
for n = 1:Nm

    % hierarchical level and parents
    %----------------------------------------------------------------------
    L   = Nm - n + 1;                        
    if L > 1
        iD  = spm_cat(mdp.MDP.id.D);         % states
        iE  = spm_cat(mdp.MDP.id.E);         % paths
    else
        iD  = 1:size(Y{L},1);
        iE  = [];
    end

    % predicted states (and entropy)
    %----------------------------------------------------------------------
    subplot(Nm + 3,2,(n - 1)*2 + 3)
    image((1 - spm_cat(Y{L}(iD,:)))*64)
    title(sprintf('Predictive posterior (states) level %i',L),'FontSize',12)

    % predicted paths (and entropy)
    %----------------------------------------------------------------------
    if L > 1
        subplot(Nm + 3,2,(n - 1)*2 + 4)
        image((1 - spm_cat(Y{L}(iE,:)))*64)
        title(sprintf('Predictive posterior (paths) level %i',L),'FontSize',12)
    end

    % next level
    %----------------------------------------------------------------------
    if n < Nm
        mdp = mdp.MDP;
    end
end

% At the lowest (voxel) level
%==========================================================================

% observed outcomes
%--------------------------------------------------------------------------
X  = spm_O2wav(O{1},WAV);

% predicted outcomes
%--------------------------------------------------------------------------
P  = spm_O2wav(Y{1},WAV);

% plot first level outcomes as a sonogram
%--------------------------------------------------------------------------
subplot(Nm + 3,1,Nm + 2)
imagesc(log(abs(P) + exp(-8)))
title(sprintf('Predictions'),'FontSize',11)

% Place sound in graphic object
%--------------------------------------------------------------------------
s      = spm_wavshow(P,WAV);
h      = get(gca,'Children');
set(h(1),'Userdata',{s,WAV.Fs})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

% plot stimuli, if provided
%--------------------------------------------------------------------------
if isfield(mdp,'S')
    subplot(Nm + 3,1,Nm + 3)
    imagesc(log(abs(X) + exp(-8)))
    title(sprintf('Stimulus'),'FontSize',11)

    % Place sound in graphic
    %----------------------------------------------------------------------
    s      = spm_wavshow(X,WAV);
    h      = get(gca,'Children');
    set(h(1),'Userdata',{s,WAV.Fs})
    set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')
end

return


% subroutines
%==========================================================================

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
