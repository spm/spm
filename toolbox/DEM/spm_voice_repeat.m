function spm_voice_repeat
% illustrates voice recognition
% FORMAT spm_voice_repeat
%
% auditory input: SAY: "Is ... there a square/triangle above/below/there ...
%
% When invoked, this routine takes an audio input to estimate the
% fundamental and formant frequencies of the speaker. It will then plot the
% estimates and segment a short sentence. The sentence can be replayed
% after being recognised, with and without lexical content and prosody (to
% hear this,cut and paste the text on the bottom of this routine into the
% command window).
%
% see also: spm_voice_speak.m and spm_voice_segmentation.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_repeat.m 7600 2019-06-01 09:30:30Z karl $


%% setup
%==========================================================================
global VOX
try
    load VOX
catch
    error('please create VOX.mat and place in path')
end


% Prompt for sentence
%==========================================================================
% SAY: "Is there a ...
%--------------------------------------------------------------------------
VOX.audio    = audiorecorder(22050,16,1);
VOX.mute     = 1;
VOX.formant  = 1;

% get fundamental and formant frequencies
%--------------------------------------------------------------------------
str{1}  = {'is'};
str{2}  = {'there'};
str{3}  = {'a'};
str{4}  = {'triangle','square'};
str{5}  = {'below','above'};
[~,P]   = spm_voice_i(str);
[F0,F1] = spm_voice_identity(VOX.audio,P);
VOX.F1  = F1;
VOX.F0  = F0;

% read audio file and segment
%--------------------------------------------------------------------------
[SEG,W,P,R] = spm_voice_read(getaudiodata(VOX.audio),P);

spm_figure('GetWin','Segmentation'); clf;
spm_voice_segmentation(VOX.audio,SEG);


%% articulate with and without prosody
%==========================================================================

% articulate: prosody without lexical content
%--------------------------------------------------------------------------
spm_voice_speak(1,P,R);

% articulate: with no prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,[],R);

spm_voice_speak(1,P,[]);

% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,P,R);



