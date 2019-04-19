function spm_voice_P300
% Creates lexical and prosody cell arrays from sound file exemplars
% FORMAT spm_voice_P300
%
% PATH         -  directory containing sound files of exemplar words
%
% saves VOX.mat
%
% VOX.LEX(w,k) -  structure array for k variants of w words
% VOX.PRO(p)   -  structure array for p aspects of prosody
% VOX.WHO(w)   -  structure array for w aspects of idenity
%
%  This routine creates structure arrays used to infer the lexical and
%  prosody of a word. It uses a library of sound files, each containing 32
%  words spoken with varying prosody. The name of the sound file labels the
%  word in question. These exemplars are then transformed (using a series
%  of discrete cosine and Hilbert transforms) into a set of parameters,
%  which summarise the lexical content and prosody. The inverse transform
%  generates  timeseries that can be played to articulate a word. The
%  transform operates on a word structure xY to create lexical and prosody
%  parameters (Q and P respectively). The accuracy of  lexical inference
%  (i.e., voice the word recognition) is assessed  using the exemplar
%  (training) set and a narrative sound file called '../test.wav' (and
%  associated '../test.txt'). The operation of each subroutine can be
%  examined using graphical outputs by selecting the appropriate options in
%  a voice recognition specific global variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_P300.m 7574 2019-04-19 20:38:15Z karl $



%% setup options and files
%==========================================================================

% directory of sound files if necessary
%--------------------------------------------------------------------------
if ~nargin
    PATH = 'C:\Users\karl\Dropbox\Papers\Voice recognition\Sound files';
end

global VOX
VOX.graphics = 0;
VOX.mute     = 1;
VOX.onsets   = 0;


%% get lexical and prosody arrays in sound file directory
%--------------------------------------------------------------------------
load VOX


%% read the first few words of a test file
%--------------------------------------------------------------------------

% prior words
%--------------------------------------------------------------------------
clear str
str{1} = {'is'};
str{2} = {'there'};
str{3} = {'a'};
str{4} = {'triangle','square'};
str{5} = {'below','above','there'};
str{6} = {'no','yes'};
str{7} = {'is','there'};
[i,P]  = spm_voice_i(str);

%% record and repeat some dictation
%==========================================================================

% record a sentence and save
%--------------------------------------------------------------------------
spm_voice_read
Y = getaudiodata(VOX.audio);
save sentence Y

% segment with priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation: with priors'); clf; 
SEG1 = spm_voice_read(Y,P);
EEG1 = spm_voice_segmentation(Y,SEG1);

% segment without priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation: no priors'); clf; 
SEG0 = spm_voice_read(Y);
EEG0 = spm_voice_segmentation(Y,SEG0);

% repeat 
%==========================================================================

% record a sentence and save
%--------------------------------------------------------------------------
spm_voice_read
Y = getaudiodata(VOX.audio);
save sentence Y

%% segment with priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation: with priors'); clf;

% VOX.noise = 1.4;

SEG1  = spm_voice_read(Y);
[i,P] = spm_voice_i({SEG1.str}');
s     = 5;
SEG1  = spm_voice_read(Y,P);
EEG1  = spm_voice_segmentation(Y,SEG1);

% segment without priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation: no priors'); clf; 
SEG0  = spm_voice_read(Y,spm_softmax(log(P)/8));
EEG0  = spm_voice_segmentation(Y,SEG0);

% add evoked responses with priors to show more exuberant ERPs
%--------------------------------------------------------------------------
subplot(4,1,4), hold on
plot(VOX.PST,EEG1,'-.');

% segment with incongruent priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation: bad priors'); clf;
Q     = P; Q(:,s) = spm_softmax(-log(P(:,s))/2);
SEG2  = spm_voice_read(Y,Q);
EEG2  = spm_voice_segmentation(Y,SEG2);

% with priors
%--------------------------------------------------------------------------
try, VOX = rmfield(VOX,'noise'); end


% P300 
%==========================================================================

% simulations of P300
%--------------------------------------------------------------------------
spm_figure('GetWin','P300'); clf;
i   = spm_voice_i(SEG1(s).str);
E0  = EEG0(:,i);
E1  = EEG1(:,i);
E2  = EEG2(:,i);

subplot(3,2,1), plot(VOX.PST,E0,'-.',VOX.PST,E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Waveforms','FontSize',16)
set(gca,'YLim',[-2 2]), axis square, spm_axis tight, box off

subplot(3,2,3), plot(VOX.PST,E0 - E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Differences','FontSize',16)
set(gca,'YLim',[-1 1]), axis square, spm_axis tight, box off

subplot(3,2,2), plot(VOX.PST,E2,'-.',VOX.PST,E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Waveforms','FontSize',16)
set(gca,'YLim',[-2 2]), axis square, spm_axis tight, box off

subplot(3,2,4), plot(VOX.PST,E2 - E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Differences','FontSize',16)
set(gca,'YLim',[-1 1]), axis square, spm_axis tight, box off

%%
%--------------------------------------------------------------------------
for s = 1:numel(SEG2)
    q     = SEG2(s).L{1};
    p     = SEG2(s).p;
    KL(s) = q'*(log(q) - log(p));
end

spm_figure('GetWin','Segmentation: bad priors')
RMS  = std(EEG2,[],2);
RMS  = max(KL)*RMS/max(RMS);
subplot(4,1,2), hold off, plot(VOX.PST,RMS), hold on
xlabel('time (sec)'), ylabel('RMS/KL (nats)')
title('Belief updating','FontSize',16)
spm_axis tight

for s = 1:numel(SEG2)
    t = SEG2(s).IT/VOX.FS;
    plot(t,KL(s),'.','MarkerSize',32)
end
