function spm_voice_P300
% Illustrate voice recognition with lexical priors
% FORMAT spm_voice_P300
%
% loads the global variable VOX.mat
%
% VOX.LEX(w,k) -  structure array for k variants of w words
% VOX.PRO(p)   -  structure array for p aspects of prosody
% VOX.WHO(w)   -  structure array for w aspects of idenity
%
% This routine demonstrates the basic functionality of voice recognition or
% active listening with a special focus on segmentation and the simulated
% neurophysiological correlates of belief updating. It starts by
% demonstrating segmentation; either in response to some spoken sentences
% (read from prompts in the script or by loading exemplar sentences). It
% then moves on to demonstrating the effect of changing the precision of
% prior beliefs about lexical content and how this is expressed in terms of
% simulated belief updating via the minimisation of variational free
% energy.
%
% This routine assumes the necessary files are located in a particular
% (Sound files) directory; that can be specified by editing the script
% below.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_P300.m 7750 2019-12-05 17:54:29Z spm $


%% demo mode loads sentence (.mat) files
%==========================================================================
DEMO = 1;                      

% get lexical and prosody arrays in sound file directory
%--------------------------------------------------------------------------
clear global VOX
global VOX

NAME  = which('VOX.mat');
PATH  = fileparts(NAME);
load(NAME)

VOX.analysis = 0;
VOX.graphics = 0;
VOX.interval = 0;
VOX.mute     = 0;
VOX.onsets   = 0;
VOX.depth    = 0;

%% set up priors a succession of 'triangle' or 'square'
%==========================================================================
nw    = numel(VOX.LEX);                                   % number of words
[i,P] = spm_voice_i(repmat({{'triangle','square'}},1,6)); % prior words   

% record corresponding sequence and save - or load a preprepared sentence
%--------------------------------------------------------------------------
if DEMO
    
    load(fullfile(PATH,'triangle_square.mat'))
    
else
    
    % SAY: {'triangle','square','triangle','square','triangle','square'};
    %----------------------------------------------------------------------
    spm_voice_read
    Y  = getaudiodata(VOX.audio);
    save(fullfile(PATH,'triangle_square.mat'),'Y');
    [FS,read] = spm_voice_FS(wfile)
    
end

% illustrate candidate intervals (and F0) for the first word
%--------------------------------------------------------------------------
VOX.FS       = spm_voice_FS;
VOX.IT       = 1;
VOX.onsets   = 1;
VOX.interval = 1;
spm_voice_get_word(Y);
VOX.interval = 0;
spm_voice_fundamental(Y,VOX.FS);
VOX.onsets   = 0;


% segment without priors
%--------------------------------------------------------------------------
VOX.depth    = 0;
spm_figure('GetWin','Segmentation: no priors');  clf;
SEG0 = spm_voice_read(Y);
EEG0 = spm_voice_segmentation(Y,SEG0);

% segment with priors
%--------------------------------------------------------------------------
VOX.depth    = 1;
spm_figure('GetWin','Segmentation: with priors'); clf; 
SEG1 = spm_voice_read(Y,P);
EEG1 = spm_voice_segmentation(Y,SEG1);
VOX.depth    = 0;


% illustrate the role of prior precision
%==========================================================================

%% lexical priors; i.e., plausible sequence of words
%--------------------------------------------------------------------------
str{1} = {'is'};
str{2} = {'there'};
str{3} = {'a'};
str{4} = {'triangle','square'};
str{5} = {'below','above','there'};
[i,P]  = spm_voice_i(str);                   % get priors
P      = spm_softmax(log(P));                % and double their precision


% record corresponding sequence and save - or load a preprepared sentence
%--------------------------------------------------------------------------
if DEMO
    
    load(fullfile(PATH,'sentence.mat'))
    
else
    
    % SAY: {'is','there','a','square','above'};
    %----------------------------------------------------------------------
    spm_voice_read
    Y  = getaudiodata(VOX.audio);
    save(fullfile(PATH,'sentence.mat'),'Y');
    
end

%% change VOX.noise to simulate speech in noise (e.g., VOX.noise = 2)
%--------------------------------------------------------------------------
VOX.noise = 1;                                  % noise level (default: 1)
s         = 3;                                  % index of surprising word

%% segment with priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation (P300): with priors'); clf;
SEG1  = spm_voice_read(Y,P);
EEG1  = spm_voice_segmentation(Y,SEG1);

% segment without priors
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation (P300): no priors'); clf; 
SEG0  = spm_voice_read(Y);
EEG0  = spm_voice_segmentation(Y,SEG0);

% add evoked responses with priors, to highlight more exuberant ERPs
%--------------------------------------------------------------------------
subplot(4,1,4), hold on
plot(VOX.PST,EEG1,'-.');

% remove priors from a single word (s)
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation: Bayesian surprise'); clf;
Q     = P; Q(:,s) = ones(nw,1)/nw;
SEG2  = spm_voice_read(Y,Q);
EEG2  = spm_voice_segmentation(Y,SEG2);


% illustrate surprise responses (e.g.,P300)
%==========================================================================

% simulations of P300: responses to a single word (s)
%--------------------------------------------------------------------------
spm_figure('GetWin','P300'); clf;
i   = spm_voice_i(SEG1(s).str);
E0  = EEG0(:,i);
E1  = EEG1(:,i);
E2  = EEG2(:,i);

% peristimulus time for plotting (250 ms before and after offset)
%--------------------------------------------------------------------------
t   = SEG2(s).IT/VOX.FS - 1/4;
x   = [t,t + 1/2,t + 1/2,t];
y   = [-1,-1,1,1];
c   = spm_softmax(rand(3,1))';

% responses with and without priors - all words
%--------------------------------------------------------------------------
subplot(3,2,1), plot(VOX.PST,E0,'-.',VOX.PST,E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Waveforms','FontSize',16)
axis square, spm_axis tight, set(gca,'YLim',[-1 1]), box off

% response differentials - all words
%--------------------------------------------------------------------------
subplot(3,2,3), plot(VOX.PST,E0 - E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Differences','FontSize',16)
axis square, spm_axis tight, set(gca,'YLim',[-1 1]/3), box off

% responses with and without priors - one word
%--------------------------------------------------------------------------
subplot(3,2,2), plot(VOX.PST,E2,'-.',VOX.PST,E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Waveforms','FontSize',16)
axis square, spm_axis tight, set(gca,'YLim',[-1 1]), box off

% show peristimulus time
%--------------------------------------------------------------------------
hold on, h = fill(x,y,c); hold off
set(h,'Facealpha',1/8,'EdgeAlpha',1/8);

% response differentials - one word
%--------------------------------------------------------------------------
subplot(3,2,4), plot(VOX.PST,E2 - E1) 
xlabel('time (sec)'), ylabel('amplitude'), title('Differences','FontSize',16)
axis square, spm_axis tight, set(gca,'YLim',[-1 1]/3), box off

% show peristimulus time
%--------------------------------------------------------------------------
hold on, h = fill(x,y/3,c); hold off
set(h,'Facealpha',1/8,'EdgeAlpha',1/8);


%% illustrate the relationship between belief updating and RMS responses
%==========================================================================
for i = 1:numel(SEG2)
    q     = spm_softmax(SEG2(i).L{1});
    p     = SEG2(i).p;
    KL(i) = q'*(log(q + exp(-8)) - log(p + exp(-8)));
end

% show RMS responses as a function of time
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation: Bayesian surprise');
RMS  = var(EEG2,[],2);
RMS  = 4*RMS/max(RMS);
subplot(4,2,3), hold off, plot(VOX.PST,RMS), hold on
xlabel('time (sec)'), ylabel('Power/KL (nats)')
title('Evoked power','FontSize',16)
spm_axis tight

% and overlay associated belief updating in terms of KL divergence
%--------------------------------------------------------------------------
for i = 1:numel(SEG2)
    t = SEG2(i).IT/VOX.FS + 1/8;
    plot([t,t],[0 KL(i)],'r','LineWidth',4)
end

i     = find(diff(RMS(1:end - 1)) > 0 & diff(RMS(2:end)) < 0);
[R,i] = sort(RMS(i),'descend');
R     = R(1:numel(KL));
KL    = sort(KL,'descend');
B     = [ones(numel(KL),1) KL(:)]\R;

subplot(4,2,4), plot(KL,R,'.r','MarkerSize',16), hold on
kl    = -1:4;
plot(kl,B(1) + kl*B(2),':r','MarkerSize',16)
xlabel('Bayesian suprise (KL - nats)'), ylabel('Evoked power')
title('Belief updating','FontSize',16)
spm_axis square
