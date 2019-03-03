function [wfile] = spm_voice_listen
% starts recording for subsequent voice recognition 
% FORMAT [wfile] = spm_voice_listen
%
% wfile  - audio object
%
%  This routine creates an audio recorder object and starts recording for
%  eight seconds. This is the use with spm_voice_get_wordthat can retrieve
%  successive words following an index or pointer in the global variable
%  called voice_options (voice_options.I0)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_listen.m 7536 2019-03-03 21:38:19Z karl $

% get timeseries from audio recorder(or from a path
%--------------------------------------------------------------------------


%% get data features from a wav file or audiorecorder object
%==========================================================================
if isa(wfile,'audiorecorder')
    stop(wfile);
    record(wfile,8);
    pause(1);
end


%% run through sound file and evaluate likelihoods
%==========================================================================
global voice_options
voice_options.I0 = 1;
str   = {};
for s = 1:5
    
    % find next word
    %----------------------------------------------------------------------
    L   = spm_voice_get_word(wfile);
        
    % break if EOF
    %----------------------------------------------------------------------
    if isempty(L)
        return
    end
    
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
    str{s} = voice_options.LEX(w,1).word       % lexical string
    
end

% stop recording audiorecorder object
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    stop(wfile);
end

%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,P,R);


