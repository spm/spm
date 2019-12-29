function spm_voice(PATH)
% Create lexical and prosody cell arrays from sound file exemplars
% FORMAT spm_voice(PATH)
%
% PATH         -  directory containing sound files of exemplar words
%                 (and a test.wav file in a subdirectory /test)
%
% saves VOX.mat
%
% VOX.LEX(w,k) -  structure array for k variants of w words
% VOX.PRO(p)   -  structure array for p aspects of prosody
% VOX.WHO(w)   -  structure array for w aspects of idenity
%
%  This routine creates structure arrays used to infer the lexical class,
%  prosody and speaker identity of a word. It uses a library of sound
%  files, each containing 32 words spoken with varying prosody. The name of
%  the sound file labels the word in question. These exemplars are then
%  transformed (using a series of discrete cosine transforms) into a set of
%  parameters, which summarise the lexical content and prosody. The inverse
%  transform generates  timeseries that can be played to articulate a word.
%  The transform operates on a word structure xY to create lexical and
%  prosody parameters (Q and P respectively). The accuracy of  lexical
%  inference (i.e., voice the word recognition) is assessed  using the
%  exemplar (training) set and a narrative sound file called '../test.wav'
%  (and associated '../test.txt'). The operation of each subroutine can be
%  examined using graphical outputs by selecting the appropriate options in
%  a voice recognition specific global variable VOX. this structure is
%  saved in the sound file for subsequent use.
%
%  Auxiliary routines will be found at the end of the script. These include
%  various optimisation schemes and illustrations of online voice
%  recognition
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice.m 7761 2019-12-29 17:48:04Z karl $


%% setup options and files
%==========================================================================

% directory of sound files if necessary
%--------------------------------------------------------------------------
clear global VOX
close all
%clear all
clc
if ~nargin
    PATH = 'C:\Users\karl\Dropbox\Papers\Voice recognition\Sound files';
end

global VOX
VOX.analysis = 0;
VOX.graphics = 0;
VOX.interval = 0;
VOX.onsets   = 0;
VOX.mute     = 1;

%% get training corpus
%==========================================================================
[xY,word,NI] = spm_voice_get_xY(PATH); save xY xY word NI


%% set VOX (LEX, PRO and WHO) and parameters; e.g., timbre  = mean(P(:,4));
%==========================================================================
spm_voice_get_LEX(xY,word,NI);


%% articulate every word under all combinations of (5 levels) of prosody
% {VOX.PRO.str}: {'amp','lat','dur','tim','Tu','Tv''Tf',Tw','p0','p1','p2'}
% {VOX.WHO.str}: {'ff0','ff1'}
%--------------------------------------------------------------------------
VOX.mute     = 0;
VOX.graphics = 1;
VOX.RAND     = 1/8;

nw    = numel(VOX.LEX);                       % number of lexical features
np    = numel(VOX.PRO(1).pE);                 % number of prosody features
k     = [1 np];
f     = fix(np/2);
for w = 1:nw
    for i = k
        for j = k
            spm_voice_speak(w,[8;1;f;f;j;f;f;f;i;f;f],[3;3]); pause(1/2)
        end
    end
end

VOX.RAND = 0;


%% invert a test file of 87 words & optimise basis function order (nu,nv)
%==========================================================================
wtest   = fullfile(PATH,'test','test.wav');
ttest   = fullfile(PATH,'test','test.txt');
spm_voice_test(wtest,ttest);


%% save structure arrays in sound file directory
%--------------------------------------------------------------------------
VOX.analysis = 0;
VOX.graphics = 0;
VOX.interval = 0;
VOX.onsets   = 0;
VOX.mute     = 1;
DIR   = fileparts(which('spm_voice.m'));
save(fullfile(DIR,'VOX'),'VOX')


%% read the first few words of a test file
%==========================================================================
VOX.graphics = 0;

% prior words (first words are the spoken words)
%--------------------------------------------------------------------------
clear str
str{1} = {'is'};
str{2} = {'there'};
str{3} = {'a'};
str{4} = {'triangle','square'};
str{5} = {'below','above'};
str{6} = {'no','yes'};
str{7} = {'is','there'};
str{8} = {'is','there'};
[i,P]  = spm_voice_i(str);

% Read test file
%--------------------------------------------------------------------------
spm_voice_read(wtest,P);



%% recover funcdamental and formant formant frequencies
%--------------------------------------------------------------------------
try
    VOX = rmfield(VOX,{'F0','F1'});
end

VOX.formant = 1;
for i = 1:4
   spm_voice_identity(wtest,P);
end


%% illustrate accuracy (i.e., inference) using training corpus
%==========================================================================

% likelihood of training set
%--------------------------------------------------------------------------
[nw,ns] = size(xY);
    
q     = [];
p     = [];
r     = 0;
for i = 1:nw
    for j = 1:16:ns
        
        % evaluate lexical (L) and prosody (M) likelihoods
        %------------------------------------------------------------------
        [L,M] = spm_voice_likelihood(xY(i,j));
        L(:)  = spm_softmax(L(:));
        M     = spm_softmax(M);
        L     = sum(L,2);
        
        % inferred and true lexical outcome
        %------------------------------------------------------------------
        q(:,end + 1) = L;
        p(i,end + 1) = 1;
        r        = r + M;
    end
end

% show results
%--------------------------------------------------------------------------
spm_figure('GetWin','Accuracy (training)'); clf
subplot(2,1,1); imagesc(q)
a     = 1 - sum((p(:) - q(:)) > 1/8)/sum(p(:));
str   = sprintf('Lexical classification accuracy %-2.0f p.c.',100*a);
title(str,'FontSize',16), xlabel('exemplars'), ylabel('words')
set(gca,'Ytick',1:nw),set(gca,'YtickLabel',{VOX.LEX.word})

subplot(3,1,3); imagesc(r'), axis image
j     = numel(VOX.PRO);
k     = numel(VOX.PRO(1).pE);
xlabel('level'), ylabel('attribute'), title('Prosody','FontSize',16)
set(gca,'Xtick',1:k),set(gca,'Ytick',1:j),set(gca,'YtickLabel',{VOX.PRO.str})


%% auxiliary code
%==========================================================================
return

%% auxiliary code - P300 demonstration
%==========================================================================
spm_voice_P300


%% auxiliary code
%==========================================================================

% record a sentence and segment with active listening
%--------------------------------------------------------------------------
% SAY: "Square Square Square Square"
%--------------------------------------------------------------------------
clear VOX
global VOX
load VOX
VOX.audio = audiorecorder(22050,16,1);
[FS,read] = spm_voice_FS(VOX.audio);

clear str
str{1} = {'is'};
str{2} = {'there a'};
str{3} = {'square','triangle','red','green'};
str{4} = {'triangle','square','below','above'};
str{5} = {'below','above'};
[i,P]  = spm_voice_i(str);

% Read test file
%--------------------------------------------------------------------------
VOX.depth = 1;
VOX.LL    = 4;
spm_voice_read(VOX.audio,P);


spm_voice_read(read(VOX.audio),P);


%% optmise regularization with respect to classification accuracy
%==========================================================================
clear all
global VOX
load VOX
PATH  = 'C:\Users\karl\Dropbox\Papers\Voice recognition\Sound files';
wtest = fullfile(PATH,'test','test.wav');
ttest = fullfile(PATH,'test','test.txt');
load(fullfile(PATH,'xY.mat'));

% reporting options
%--------------------------------------------------------------------------
VOX.graphics = 0;
VOX.mute     = 1;
VOX.onsets   = 0;

% search over variables
%--------------------------------------------------------------------------
E     = [1 2 4 8 16]*2;
G     = [1 2 4 8 16]/2;
for i = 1:numel(E)
    
    try VOX = rmfield(VOX,'E'); end
    try VOX = rmfield(VOX,'G'); end
    % VOX.E = E(i);
    VOX.G = G(i);
    P     = spm_voice_get_LEX(xY,word,NI);
    A(i)  = spm_voice_test(wtest,ttest)   
end

subplot(2,3,6), plot(E,A,'.','MarkerSize',16), axis square
title('Accuracy (E)','FontSize',16),drawnow



%% iterative inversion
%==========================================================================
global VOX
load VOX
VOX.analysis = 0;
VOX.graphics = 1;
VOX.interval = 0;
VOX.formant  = 0;
VOX.mute     = 1;
VOX.rf       = 0;
VOX.RF       = 0;

% get parameters for a particular word
%--------------------------------------------------------------------------
xY       = spm_voice_speak(12);
xY.P.lat = -8;
P        = xY.P;

% compose and decompose iteratively
%--------------------------------------------------------------------------
for i = 1:4
    
    % generate timeseries and play
    %----------------------------------------------------------------------
    Y        = spm_voice_iff(xY);
    sound(Y,spm_voice_FS)
    spm_voice_fundamental(Y,spm_voice_FS);
    
    % map back to parameters by inverting
    %----------------------------------------------------------------------
    xY       = spm_voice_ff(Y);
    xY.P.amp = P.amp;
    xY.P.lat = P.lat;
    xY.P.dur = P.dur;

    disp(xY.P)

end



%% Illustrate mapping from formant space to coeficients
%==========================================================================
spm_figure('GetWin','Formant priors'); clf
VOX.analysis = 0;
VOX.graphics = 0;
VOX.interval = 0;
for w = 1:4
    
    % Expectations for this word
    %----------------------------------------------------------------------
    xy      = spm_voice_speak(w);    
    [Q,U,V] = spm_voice_Q(xy.W,xy.P.pch);
    
    subplot(4,2,2*w - 1), imagesc(Q)
    xlabel('time (seconds)'), ylabel('Formants (Hz)')
    title(sprintf('Mean (%s)',VOX.LEX(w).word),'FontSize',16)
    
    % Expectations for this word
    %----------------------------------------------------------------------
    L     = Q;
    for i = 1:numel(Q)
        dQ    = Q;
        dQ(i) = dQ(i) + 1/128;
        xy.W  = U\dQ/V';
        dL    = spm_voice_likelihood(xy,w);
        L(i)  = dL(w);
    end
    
    subplot(4,2,2*w), imagesc(L)
    xlabel('time (seconds)'), ylabel('Formants (Hz)')
    title('Sensitivity','FontSize',16), drawnow
    
end



%% estimate formant frequency via maximum likelihood
%==========================================================================
clear VOX
global VOX
load VOX

if 0
    
    % Series of single words
    %----------------------------------------------------------------------
    load Ann
    str = {'yes','yes','yes','yes'};
else
    
    %  spoken sentence
    %----------------------------------------------------------------------
    load Annagain
    clear str
    str{1} = {'is'};
    str{2} = {'there'};
    str{3} = {'a'};
    str{4} = {'square'};
    str{5} = {'below','above','there'};
end

VOX.formant  = 1;
VOX.FS       = 22050;
sound(Y,VOX.FS)

% estimate fundamental and formant frequencies
%--------------------------------------------------------------------------
[i,P] = spm_voice_i(str);
clear f0 f1
for i = 1:2
   spm_voice_identity(Y,P);
end

% mimic speaker
%--------------------------------------------------------------------------
VOX.mute    = 0;
VOX.rf      = 0;
VOX.RAND    = 1/8;
VOX.DT      = 1/32;
VOX.depth   = 0;

[SEG,p,q,r] = spm_voice_read(Y,P);
spm_voice_segmentation(Y,SEG);

% Now segment the synthetic speech
%--------------------------------------------------------------------------
[xY,y]      = spm_voice_speak(p,q,r);
spm_voice_read(y,P);



%% interactive pitch analysis
%==========================================================================
global VOX
load VOX
VOX.graphics = 1;
VOX.mute     = 0;

F0 = exp(VOX.WHO(1).pE + VOX.R.F0);  % range of fundamental frequency
F1 = exp(VOX.WHO(2).pE + VOX.R.F1);  % range of formant frequency
w  = 14;                             % spoken word
spm_voice_speak(w);

% collect points over F0 and F1
%--------------------------------------------------------------------------
clear ff0 ff1
for i = 1:128
    subplot(4,2,6), imagesc(ones(numel(F0),numel(F1)))
    xlabel('F0'),ylabel('F1')
    [f0,f1] = ginput(1);
    ff0(i)  = fix(f0);
    ff1(i)  = fix(f1);
    try
        spm_voice_speak(w,8,[ff0(i);ff1(i)]);
    catch
        ff0(i)  = [];
        ff1(i)  = [];
        break
    end
end

% linear regression and graphics
%--------------------------------------------------------------------------
Y   = F1(ff1(:));
X   = [ones(numel(ff0),1) F0(ff0(:))];
B   = X\Y;
disp('coefficients: F1 = B(1) + B tF0*B(2)'); disp(B)
disp('standard deviation'); disp(std(Y - X*B))

subplot(4,2,6)
plot(F0(ff0),F1(ff1),'.',F0,F1,':',F0,B(1) + F0*B(2),'-.')
xlabel('F0'),ylabel('F1'),title('selected points')

subplot(4,2,5)
plot(ff0,ff1,'.',1:numel(F0),1:numel(F1),':')
xlabel('F0'),ylabel('F1'),title('relationship between F1 and F2')


%% auxiliary code
%==========================================================================

%% load script and prompt for audio file: 32 words, at one word per second
%--------------------------------------------------------------------------
str   = textread(ttest,'%s');
str   = unique(str);
for s = 1:numel(str)
    clc,  disp([str(s) 'new file']), pause
    for i = 1:32
        clc,disp(str(s));     pause(0.9)
        clc,disp('    -*- '); pause(0.1)
    end
end

%% add word
%==========================================================================

%% prompt for audio file: 32 words, at one word per second
%--------------------------------------------------------------------------
str   = 'please';
audio = audiorecorder(16000,16,1);
[FS,read] = spm_voice_FS(audio);
stop(audio);

% countdown
%--------------------------------------------------------------------------
record(audio,38);
for i = 3:-1:1
    clc,disp(i);          pause(0.5)
    clc,disp('    -*- '); pause(0.5)
end

% record
%--------------------------------------------------------------------------
for i = 1:32
    clc,disp(str);        pause(0.5)
    clc,disp('    -*- '); pause(0.5)
end

% save
%--------------------------------------------------------------------------
PATH  = 'C:\Users\karl\Dropbox\Papers\Voice recognition\Sound files';
Y     = read(audio);
try
    delete(fullfile(PATH,[str '.wav']))
end
audiowrite(fullfile(PATH,[str '.wav']),Y,FS)
