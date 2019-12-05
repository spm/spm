function spm_voice_repeat
% Illustrates voice recognition
% FORMAT spm_voice_repeat
%
% When invoked, this routine takes an audio input to estimate the
% fundamental and formant frequencies of the speaker. It will then plot the
% estimates and segment a short sentence. The sentence can be replayed
% after being recognised, with and without lexical content and prosody.
% this routinely uses dialogue boxes to step through the various
% demonstrations.
%
% See also: spm_voice_speak.m and spm_voice_segmentation.m
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_repeat.m 7750 2019-12-05 17:54:29Z spm $


%% setup
%==========================================================================
global VOX
if ~isfield(VOX,'LEX')
    try
        load VOX
        VOX.analysis = 0;
        VOX.graphics = 0;
        VOX.mute     = 0;
    catch
        error('please create VOX.mat and place in path')
    end
end

%% Prompt for sentence
%==========================================================================
VOX.audio    = audiorecorder(22050,16,1);
VOX.mute     = 1;
VOX.formant  = 1;
VOX.disp     = 0;
VOX.RAND     = 1/2;
try
    VOX      = rmfield(VOX,'F0');
end

% get fundamental and formant frequencies
%--------------------------------------------------------------------------
prompt  = 'please say "is there a square above" after the beep';
uiwait(msgbox(prompt,'modal'))
beep, pause(1/2), 
stop(VOX.audio)

VOX.F0  = 200;
str     = {'is','there','a','square','above'};
[W,L]   = spm_voice_i(str);
[F0,F1] = spm_voice_identity(VOX.audio,L);
VOX.F1  = F1;
VOX.F0  = F0;

% read audio file and segment
%--------------------------------------------------------------------------
[SEG,W,P,R] = spm_voice_read(getaudiodata(VOX.audio),L);
uiwait(msgbox('show segmentation','modal'))
spm_figure('GetWin','Segmentation'); clf;
spm_voice_segmentation(VOX.audio,SEG);


%% articulate with and without prosody
%==========================================================================

% articulate: prosody without lexical content
%--------------------------------------------------------------------------
uiwait(msgbox('play without lexical content','modal'))
spm_voice_speak(W(1),P,R);  pause(1)

% articulate: with no prosody
%--------------------------------------------------------------------------
uiwait(msgbox('play without prosody','modal'))
spm_voice_speak(W,[],R(:,1)); pause(1)

% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
uiwait(msgbox('play with prosody','modal'))
spm_voice_speak(W,P,R);   pause(1)

% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
uiwait(msgbox('change lexical content','modal'))
W  = spm_voice_i({'is','there','a','triangle','below'});
spm_voice_speak(W(1:size(R,2)),P,R); pause(1)

% articulate: without priors
%--------------------------------------------------------------------------
uiwait(msgbox('play without priors','modal'))
VOX.mute = 0;
spm_voice_read(getaudiodata(VOX.audio));
