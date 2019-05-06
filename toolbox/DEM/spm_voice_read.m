function [SEG,W,P,R] = spm_voice_read(wfile,L,N)
% Reads and translates a sound file or audio source
% FORMAT [SEG,W,P,R] = spm_voice_read(wfile,[L],[N])
%
% wfile  - .wav file or audio object or (double) timeseries
% L      - prior likelihood of lexical content
% N      - number of words to read
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
% This routine takes a sound file has an input and infers the lexical
% content, prosody and speaker. In then articulates the phrase or
% sequence of word segments (SEG). If called with no output arguments it
% generates graphics detailing the segmentation.
%
% see also: spm_voice_speak.m and spm_voice_segmentation.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_read.m 7588 2019-05-06 21:26:32Z karl $


%% setup
%==========================================================================
global VOX
if ~isfield(VOX,'LEX')
    try
        load VOX
        VOX.analysis = 1;
        VOX.graphics = 0;
        VOX.mute     = 0;
    catch
        error('please create VOX.mat and place in path')
    end
end


% if audio sources is not provided, assume a new recording
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

%% priors, assuming at most eight words
%--------------------------------------------------------------------------
try ns = N; catch, ns = 8; end

nw    = numel(VOX.LEX);
for s = 1:ns
    try
        p(:,s) = L(:,s);
    catch
        p(:,s) = ones(nw,1)/nw;
    end
end

%% run through sound file and evaluate likelihoods
%==========================================================================
VOX.I0 = 1;                                    % first index
VOX.IT = 1;                                    % final index
for s  = 1:ns
    
    % find next word
    %----------------------------------------------------------------------
    L   = spm_voice_get_word(wfile,p(:,s));
        
    % break if EOF
    %----------------------------------------------------------------------
    if isempty(L), break, end
    
    % identify the most likely word and prosody
    %----------------------------------------------------------------------
    [d,w]  = max(L{1});                        % most likely word
    [d,c]  = max(L{2});                        % most likely prosody
    [d,r]  = max(L{3});                        % most likely identity
    
    % string
    %----------------------------------------------------------------------
    SEG(s).str = VOX.LEX(w,1).word;            % lexical string
    SEG(s).I0  = VOX.I0;                       % first
    SEG(s).IT  = VOX.IT;                       % final
    SEG(s).p   = p(:,s);                       % prior
    SEG(s).L   = L;                            % posteriors
    SEG(s).W   = w(:);                         % lexical class
    SEG(s).P   = c(:);                         % prosody classes
    SEG(s).R   = r(:);                         % speaker class

    disp({SEG.str})
    
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
if ~VOX.mute
    spm_voice_speak(W,P,R);
end

%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
if ~nargout
    spm_figure('GetWin','Segmentation'); clf;
    spm_voice_segmentation(wfile,SEG);
end



