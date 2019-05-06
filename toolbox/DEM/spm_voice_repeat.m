function spm_voice_repeat
% Reads and translates a sound file or audio source
% FORMAT spm_voice_read(wfile,[L],[N])
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
% $Id: spm_voice_repeat.m 7587 2019-05-06 16:47:53Z karl $


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

% Prompt for sentence
%==========================================================================
% SAY: "Square Square Square Square"
%--------------------------------------------------------------------------
VOX.audio = audiorecorder(22050,16,1);



VOX.analysis = 1;
VOX.graphics = 1;
VOX.formant  = 1;
VOX.onsets   = 0;
VOX.mute     = 1;

VOX.FS       = 22050;
VOX.IT       = 1;


sound(Y,VOX.FS)
str       = {'yes'};
w         = spm_voice_i(str);
[i,P]     = spm_voice_i(str);
P         = spm_softmax(log(P));
[L,F0,F1] = spm_voice_identity(Y,P);

VOX.F1       = F1;
VOX.F0       = F0;
VOX.mute     = 0;
VOX.onsets   = 0;
VOX.formant  = 0;
spm_voice_read(Y);

VOX = rmfield(VOX,{'F0','F1','FS'});
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
W      = [];                                   % posterior over words
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

%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
if ~nargout
    spm_figure('GetWin','Segmentation'); clf;
    spm_voice_segmentation(wfile,SEG);
end


%% auxiliary code to articulate with and without prosody
%==========================================================================

%% articulate: prosody without lexical content
%--------------------------------------------------------------------------
spm_voice_speak([],P,1,LEX,PRO,WHO);


%% articulate: with no prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,[],1,LEX,PRO,WHO);


%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,P,1,LEX,PRO,WHO);



