function MDP = DEM_sound_compression
% Structure learning from sound files
%__________________________________________________________________________
%
% This routine illustrates the use of renormalizable generative models to
% recognise, compress and generate sound files or auditory streams. This
% application uses a blocking or grouping operator, which operates on the
% neighbouring frequencies of a sonogram (i.e., time frequency
% representation of a sound file). This means that local frequencies are
% grouped together over time and assimilated using fast structure learning.
%
% From the perspective of the generative model, the model generates
% sequences of sequences that can be regarded as successively higher
% abstractions or compositions of frequency glides. The ensuing generative
% models can then generate sequences of sequences in a way that emulates
% previous simulations of birdsong with hierarchically composed non-linear
% dynamical (chaotic) systems but are generated in discrete state spaces
% via the explicit representation of paths or trajectories.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Copyright free files obtained from the following website:
% https://sound-effects.bbcrewind.co.uk

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%--------------------------------------------------------------------------
rng(1)

% Load a sound file
%==========================================================================
START  = 50000;
END    = 500000;
% info = audioinfo('C:\Users\Karl\Dropbox\matlab\Crossbill.wav');
[s,Fs] = audioread('C:\Users\Karl\Dropbox\matlab\Crossbill.wav',[START,END]);
s      = s(:,1);

N      = fix(Fs/8000);
s      = downsample(s,N);
Fs     = Fs/N;

% ensure length can be recursively bipartitioned
%--------------------------------------------------------------------------
t      = 2.^(1:16);
ni     = find(t <= numel(s),1,'last');
s      = s(1:t(ni));

% Read into working memory and play
%==========================================================================
spm_figure('GetWin','Sound'); clf
spm_wav2I = @(s,WAV)real(spm_wft(s,WAV.k,WAV.n));

% Crop and resize to 128 x128 pixels
%--------------------------------------------------------------------------
WAV.Fs = Fs;                                   % sample rate
WAV.n  = fix(Fs*32/1000);                      % window (ms)
WAV.k  = fix(linspace(40,4000,64)/(Fs/WAV.n)); % frequencies (Hz)

% Map from image to discrete state space (c.f., Amortisation) 
%--------------------------------------------------------------------------
R      = 2.^(ni - 6);       % length of time bin
WAV.nd = 4;                 % Size of tiles in bins        [default: 4]
WAV.nb = 5;                 % Number of discrete singular variates  [9]  
WAV.mm = 16;                % Maximum number of singular modes      [32]
WAV.su = 16;                % Singular value (normalised) threshold [16]
WAV.R  = R;                 % temporal resampling                   [128]

I         = spm_wav2I(s,WAV);
[O,L,WAV] = spm_wav2O(I,WAV);

% And show the CWT images generated from a discrete representation
%--------------------------------------------------------------------------
I = spm_O2wav(O,WAV);
spm_wavshow(I,WAV);

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

RDP.T     = 16;
RDP.D{1}  = sparse(1,1,1,Ns(1),1);
RDP.E{1}  = sparse(1,1,1,Nu(1),1);
PDP       = spm_MDP_VB_XXX(RDP);

% Illustrate the model in generative mode
%--------------------------------------------------------------------------
spm_show_WAV(PDP,WAV)

% Biomimetic characterisation in terms of evoked responses
%==========================================================================
spm_figure('GetWin','Evoked responses'); clf

% Create partial stimulus (by briefly presenting the second call)
%--------------------------------------------------------------------------
S     = cell(size(O,1),128);
for g = 1:size(S,1)
    for t = 1:size(S,2)
        if t < 16
            S{g,t} = O{g,t + 32};
        else
            S{g,t} = spm_dir_norm(ones(size(O{g,1})));
        end
    end
end

% Invert deep (RG) model
%--------------------------------------------------------------------------
MDP{1}.S  = S;
RDP       = spm_mdp2rdp(MDP);
RDP.T     = 16;
RDP.D{1}  = ones(Ns(1),1);
RDP.E{1}  = ones(Nu(1),1);
PDP       = spm_MDP_VB_XXX(RDP);

% Illustrate responses
%--------------------------------------------------------------------------
spm_show_WAV(PDP,WAV)

return





