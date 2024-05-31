function MDP = DEM_music_compression
% Structure learning from sound files with active learning
%__________________________________________________________________________
%
% This routine illustrates the use of renormalizable generative models to
% recognise, compress and generate sound files or auditory streams. This
% application uses a blocking or grouping operator, which operates on the
% neighbouring frequencies of a sonogram (i.e., time frequency
% representation of a sound file). This means that local frequencies are
% grouped together over time and assimilated using fast structure learning.
%
% This particular demonstration ingests a sound file of music discretising
% it into a few hundred time-frequency voxels, which are subsequently
% re-normalised to a succession of (high-level) episodes corresponding
% roughly to a bar of music. This demonstration includes active learning,
% after structure learning to illustrate how transitions from the final
% episode to a subsequent episode can be learnt in an experience -dependent
% fashion.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Copyright free file obtained from the following website:
% https://sound-effects.bbcrewind.co.uk/

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%--------------------------------------------------------------------------
rng(1)

% Load a sound file (about 16 bars of Jazz piano)
%==========================================================================
% info = audioinfo('C:\Users\Karl\Dropbox\matlab\piano-112623.mp3');
[s,Fs] = audioread('C:\Users\Karl\Dropbox\matlab\piano-112623.mp3');

% ensure length can be recursively bipartitioned
%--------------------------------------------------------------------------
s      = s((1:1580000) + 90000,1);
N      = fix(Fs/8000);
s      = downsample(s,N);
Fs     = Fs/N;

% Read and discretise
%==========================================================================
spm_figure('GetWin','Sound'); clf
spm_wav2I = @(s,WAV)real(spm_wft(s,WAV.k,WAV.n));

% Crop and resize to 128 x128 pixels
%--------------------------------------------------------------------------
Nk     = 32;                                    % frequency bins
WAV.Fs = Fs;                                    % sample rate
WAV.n  = fix(Fs*8/1000);                        % window (8 ms)
WAV.k  = fix(linspace(1,4000,Nk)/(Fs/WAV.n));   % frequencies (Hz)

% Map from CWT image to discrete state space (c.f., Amortisation) 
%--------------------------------------------------------------------------
No     = 64;                  % number of outcomes in time series
WAV.nd = 4;                   % Size of tiles in bins        [default: 4]
WAV.nb = 9;                   % Number of discrete singular variates  [9]  
WAV.mm = 32;                  % Maximum number of singular modes      [32]
WAV.su = 16;                  % Singular value (normalised) threshold [16]
WAV.R  = fix(numel(s)/No);    % temporal resampling                   [512]

I         = spm_wav2I(s,WAV);
[O,L,WAV] = spm_wav2O(I,WAV);

% And show CWT images generated from a discrete representation
%--------------------------------------------------------------------------
I   = spm_O2wav(O,WAV);
spm_wavshow(I,WAV);

% Use the ensuing sequence for (RG) structure learning
%--------------------------------------------------------------------------
MDP = spm_MB_structure_learning(O,L);

% learn from subsequent episodes
%==========================================================================

% active learning (with minimal forgetting)
%--------------------------------------------------------------------------
FIX.A = 1;                             % disable likelihood learning
FIX.B = 0;                             %  enable transition learning
for m = 1:numel(MDP)
    MDP{m}.beta = 4;
    MDP{m}.eta  = 512;
end

% training sequences
%--------------------------------------------------------------------------
S     = cell(1,2);
seg   = No/4;
for i = 1:numel(S)
    t    = [((1:seg) + 3*seg) ((1:seg) + 2*(i - 1)*seg)];
    S{i} = O(:,t);
end

% active learning
%--------------------------------------------------------------------------
Nm    = numel(MDP);
RDP   = spm_mdp2rdp(MDP,0,0,2,FIX);
RDP.T = fix(size(S{1},2)/(2^(Nm - 1)));
for i = 1:numel(S)
    RDP  = spm_RDP_O(RDP,S{i});
    PDP  = spm_MDP_VB_XXX(RDP);
    RDP  = spm_RDP_update(RDP,PDP);
end

% simple selecton (BMR)
%--------------------------------------------------------------------------
PDP   = spm_MDP_VB_XXX(RDP);
RDP   = spm_RDP_update(RDP,PDP,'SIMPLE');

% Generate music from ensuing generative model
%==========================================================================

% remove stimulus
%--------------------------------------------------------------------------
T     = fix(size(O,2)/(2^(Nm - 1)));
RDP.T = 2*T;
RDP   = spm_RDP_O(RDP,[]);
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate the model in generative mode
%--------------------------------------------------------------------------
spm_figure('GetWin','Generative AI'); clf
spm_show_WAV(PDP,WAV)


% generate (i.e., predict) while listening to music
%==========================================================================

% provide orignal music
%--------------------------------------------------------------------------
RDP.T = T;
RDP   = spm_RDP_O(RDP,O);
PDP   = spm_MDP_VB_XXX(RDP);

% Illustrate the model in accompaniment mode
%--------------------------------------------------------------------------
spm_figure('GetWin','with accompaniment'); clf
spm_show_WAV(PDP,WAV)

T     = No/(2^(Nm - 1));   % number of episodes
bars  = 2.25;              % duration of a bar in seconds
fprintf('duration of music %0.2f secs\n', numel(s)/Fs)
fprintf('duration of music %0.2f bars\n', numel(s)/Fs/bars)
fprintf('duration of voxel %0.2f secs\n', numel(s)/Fs/No)
fprintf('duration of voxel %0.2f bars\n', numel(s)/Fs/No/bars)
fprintf('duration of event %0.2f secs\n', numel(s)/Fs/T)
fprintf('duration of event %0.2f bars\n', numel(s)/Fs/T/bars)

return






