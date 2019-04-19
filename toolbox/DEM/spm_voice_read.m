function SEG = spm_voice_read(wfile,L)
% Reads and translates a sound file 
% FORMAT spm_voice_read(wfile,[L])
%
% wfile  - .wav file or audio object or (double) timeseries
% L      - prior likelihood of lexical content
%
% requires the following in the global variable VOX:
% LEX    - lexical structure array
% PRO    - prodidy structure array
% WHO    - speaker structure array
%
%
% SEG(s).str - lexical class
% SEG(s).p   - prior
% SEG(s).L   - posterior
% SEG(s).P   - prosody class
% SEG(s).R   - speaker class
% SEG(s).I0  - first index
% SEG(s).IT  - final index
%     
% This routine takes a sound file has an input and infers the lexical
% content, prosody and speaker. In then articulates the phrase or
% sequence of word segments (SEG)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_read.m 7574 2019-04-19 20:38:15Z karl $

% get timeseries from audio recorder(or from a path
%--------------------------------------------------------------------------


%% get data features from a wav file or audiorecorder object
%==========================================================================
global VOX
if ~nargin
    wfile     = audiorecorder(22050,16,1);
    VOX.audio = wfile;
end

if isa(wfile,'audiorecorder')
    stop(wfile);
    record(wfile,8);
end

%% priors
%--------------------------------------------------------------------------
ns    = 8;
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
W      = [];                                   % posterior over words
for s  = 1:8
    
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
    W(1,s) = w(:);                             % lexical class
    P(:,s) = c(:);                             % prosody classes
    R(:,s) = r(:);                             % speaker classes
    
    % string
    %----------------------------------------------------------------------
    SEG(s).str = VOX.LEX(w,1).word;            % lexical string
    SEG(s).p   = p(:,s);                       % prior
    SEG(s).L   = L;                            % posteriors
    SEG(s).P   = P(:,s);                       % prosody
    SEG(s).R   = R(:,s);                       % speaker
    SEG(s).I0  = VOX.I0;                       % first
    SEG(s).IT  = VOX.IT;                       % final
    
    disp({SEG.str})
    
end

% stop recording audiorecorder object and return if silence
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    stop(wfile);
end
if isempty(W)
    return
end

%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
if ~VOX.mute
    spm_voice_speak(W,P,R);
end


