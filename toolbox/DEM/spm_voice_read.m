function spm_voice_read(wfile)
% Reads and translates a sound file 
% FORMAT spm_voice_read(wfile)
%
% wfile  - .wav file or audio object

% requires the following in the global variable VOX:
% LEX    - lexical structure array
% PRO    - prodidy structure array
% WHO    - speaker structure array
%
%  This routine takes a sound file has an input and infers the lexical
%  content, prosody and speaker. In then particulates the phrase or
%  sequence of words
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_read.m 7546 2019-03-18 11:02:22Z karl $

% get timeseries from audio recorder(or from a path
%--------------------------------------------------------------------------


%% get data features from a wav file or audiorecorder object
%==========================================================================
if isa(wfile,'audiorecorder')
    stop(wfile);
    record(wfile,8);
end

%% run through sound file and evaluate likelihoods
%==========================================================================
global VOX
VOX.I0 = 1;
VOX.IT = [1,1];
for s  = 1:8
    
    % find next word
    %----------------------------------------------------------------------
    L   = spm_voice_get_word(wfile);
        
    % break if EOF
    %----------------------------------------------------------------------
    if isempty(L), break, end
    
    % identify the most likely word and prosody
    %----------------------------------------------------------------------
    [d,w]  = max(L{1});                        % most likely word
    [d,p]  = max(L{2});                        % most likely prosody
    [d,r]  = max(L{3});                        % most likely identity
    W(1,s) = w(:);                             % lexical class
    P(:,s) = p(:);                             % prosody classes
    R(:,s) = r(:);                             % speaker classes
    
    % string
    %----------------------------------------------------------------------
    SEG(s).str = VOX.LEX(w,1).word;            % lexical string
    SEG(s).I0  = VOX.I0;                       % centre
    SEG(s).IT  = VOX.IT;                       % range
    
   disp({SEG.str})
    
end

% stop recording audiorecorder object
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    stop(wfile);
end

%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,P,R);


%% segmentation graphics
%--------------------------------------------------------------------------
spm_voice_segmentation(wfile,SEG)









