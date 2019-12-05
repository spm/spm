function [SEG,W,P,R] = spm_voice_read(wfile,P)
% Read and translate a sound file or audio source
% FORMAT [SEG,W,P,R] = spm_voice_read(wfile,[P])
%
% wfile  - .wav file or audio object or (double) timeseries
% P      - prior likelihood of lexical content or
%        - number of words to read (N or size(P,2))
%
% requires the following in the global variable VOX:
% LEX    - lexical structure array
% PRO    - prodidy structure array
% WHO    - speaker structure array
%
% for each (s-th) word:
%
% SEG(s).str - lexical class
% SEG(s).p   - prior
% SEG(s).L   - posterior
% SEG(s).W   - lexical class
% SEG(s).P   - prosody class
% SEG(s).R   - speaker class
% SEG(s).I0  - first index
% SEG(s).IT  - final index
%
% This routine takes a sound file or audio stream as an input and infers the lexical
% content and prosody. In then articulates the phrase or
% sequence of word segments (SEG). If called with no output arguments it
% generates graphics detailing the segmentation. This routine assumes that
% all the variables in the VOX structure are set appropriately;
% especially, the fundamental and first formant frequencies (F0 and F1)
% appropriate for speaker identity. If called with no inputs, it will
% create an audio recorder object and record dictation for a few seconds.
%
% see also: spm_voice_speak.m and spm_voice_segmentation.m
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_read.m 7750 2019-12-05 17:54:29Z spm $


%% setup
%==========================================================================
global VOX
if ~isfield(VOX,'LEX')
    try
        load VOX
        VOX.analysis = 0;
        VOX.graphics = 0;
    catch
        error('please create VOX.mat and place in path')
    end
end
if ~isfield(VOX,'disp '), VOX.disp  = 1; end
if ~isfield(VOX,'depth'), VOX.depth = 0; end

% if audio source is not provided, assume a new recording
%--------------------------------------------------------------------------
if ~nargin
    try
        wfile     = VOX.audio;
    catch
        wfile     = audiorecorder(22050,16,1);
        VOX.audio = wfile;
    end
end


% record timeseries from audio recorder, for 8 seconds
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    stop(wfile);
    record(wfile,8);
end

%% priors, if specified (uninformative priors otherwise)
%--------------------------------------------------------------------------
FS  = spm_voice_FS(wfile);
nw  = numel(VOX.LEX);
if nargin < 2
    ns = 8;
    P  = ones(nw,ns)/nw;
elseif isscalar(P)
    ns = P;
    P  = ones(nw,ns)/nw;
else
    ns = size(P,2);
end


%% run through sound file and evaluate likelihoods
%==========================================================================
VOX.IT = 1;                                     % final index

for s  = 1:ns
    
    % (deep) search
    %----------------------------------------------------------------------
    IT        = VOX.IT;
    j         = s:min(size(P,2),s + VOX.depth);
    [L,I,J,F] = spm_voice_get_word(wfile,P(:,j));
 
    % break if EOF
    %----------------------------------------------------------------------
    if isempty(L)
        break
    end
       
    % move to end of word
    %----------------------------------------------------------------------
    VOX.IT = J(2);
        
    % most likely word and prosody
    %----------------------------------------------------------------------
    [d,w]  = max(L{1});                        % most likely word
    [d,q]  = max(L{2});                        % most likely prosody
    [d,r]  = max(L{3});                        % most likely identity
    
    % string and intervals
    %----------------------------------------------------------------------
    SEG(s).str = VOX.LEX(w).word;              % lexical string
    SEG(s).I0  = J(1);                         % first
    SEG(s).IT  = J(2);                         % final
    SEG(s).I   = I;                            % extrema
    SEG(s).p   = P(:,s);                       % prior
    SEG(s).L   = L;                            % posteriors
    SEG(s).W   = w(:);                         % lexical class
    SEG(s).P   = q(:);                         % prosody classes
    SEG(s).R   = r(:);                         % speaker class
    
    % display string
    %----------------------------------------------------------------------
    if VOX.disp, disp({SEG.str}), end
    
end


% stop recording audiorecorder object and return if silence
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    stop(wfile);
end
if exist('SEG','var')
    W   = full(spm_cat({SEG.W}));
    P   = full(spm_cat({SEG.P}));
    R   = full(spm_cat({SEG.R}));
else
    SEG = [];
    W   = [];
    P   = [];
    R   = [];
    return
end

%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
if ~nargout || ~VOX.mute
    spm_voice_speak(W,P,R);
end

%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
if ~nargout
    spm_figure('GetWin','Segmentation'); clf;
    spm_voice_segmentation(wfile,SEG);
end
