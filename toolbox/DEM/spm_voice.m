function spm_voice(PATH)
% Creates lexical and prosody cell arrays from sound file exemplars
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice.m 7583 2019-05-02 12:10:08Z karl $


%% setup options and files
%==========================================================================

% directory of sound files if necessary
%--------------------------------------------------------------------------
clear global VOX
if ~nargin
    PATH = 'C:\Users\karl\Dropbox\Papers\Voice recognition\Sound files';
end

global VOX
VOX.analysis = 0;
VOX.graphics = 0;
VOX.interval = 0;
VOX.mute     = 1;
VOX.onsets   = 1;


%% get training corpus
%==========================================================================
[xY,word] = spm_voice_get_xY(PATH);

%% get VOX (LEX, PRO and WHO) and parameters; e.g., 1/F0  = mean(P(:,6));
%==========================================================================
P     = spm_voice_get_LEX(xY,word);


%% articulate every word under all combinations of (5 levels) of prosody
VOX.mute  = 0;
% {VOX.PRO.str}: 'amp'  'lat'  'dur'  'tim'  'ff0'  'inf'  'bif'
% {VOX.WHO.str}: 'ff1'
%--------------------------------------------------------------------------
nw    = numel(VOX.LEX);                       % number of words
k     = [3 6];                                % number of prosody features
for w = 1:nw
    for i = k
        for j = k
            spm_voice_speak(w,[8;1;5;5;1;i;j],3);
        end
    end
end


%% invert a test file of 87 words & optimise basis function order (ny,nv)
%==========================================================================
try
    VOX = rmfield(VOX,'nu');
    VOX = rmfield(VOX,'nv');
end
wtest   = fullfile(PATH,'test','test.wav');
ttest   = fullfile(PATH,'test','test.txt');
spm_voice_test(wtest,ttest);

%% save structure arrays in sound file directory
%--------------------------------------------------------------------------
DIR   = fileparts(which('spm_voice.m'));
save(fullfile(DIR,'VOX'),'VOX')


%% read the first few words of a test file
%==========================================================================
VOX.graphics = 0;

% prior words
%--------------------------------------------------------------------------
clear str
str{1} = {'is'};
str{2} = {'there'};
str{3} = {'a'};
str{4} = {'triangle','square'};
str{5} = {'below','above'};
str{6} = {'no','yes'};
str{7} = {'is','there'};
[i,P]  = spm_voice_i(str);

% Read test file
%--------------------------------------------------------------------------
spm_voice_read(wtest,P);


%% illustrate accuracy (i.e., inference) using training corpus
%==========================================================================

% likelihood of training set
%--------------------------------------------------------------------------
[nw,ns] = size(xY);
    
q     = [];
p     = [];
r     = 0;
J     = 0;
for i = 1:nw
    for j = 1:ns
        
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

return


%% auxiliary code
%==========================================================================

% record a sentence and segment with active listening
%--------------------------------------------------------------------------
% SAY: "Square Square Square Square"
%--------------------------------------------------------------------------
VOX.audio = audiorecorder(22050,16,1);
spm_voice_read


%% record for four seconds and illustrate recognition using priors
%==========================================================================
% SAY: "Is there a square below?"
%--------------------------------------------------------------------------
stop(VOX.audio);
record(VOX.audio,4);


%% run through sound file and evaluate likelihoods
%==========================================================================
clear W Q R SEG
VOX.IT = 1;                                    % current index
for s  = 1:size(P,2)
    
    % find next word
    %----------------------------------------------------------------------
    L   = spm_voice_get_word(VOX.audio,P(:,s));
        
    % break if EOF
    %----------------------------------------------------------------------
    if isempty(L), break, end
    
    % identify the most likely word and prosody
    %----------------------------------------------------------------------
    [d,w]  = max(L{1});                        % most likely word
    [d,c]  = max(L{2});                        % most likely prosody
    [d,r]  = max(L{3});                        % most likely identity
    W(1,s) = w(:);                             % lexical class
    Q(:,s) = c(:);                             % prosody classes
    R(:,s) = r(:);                             % speaker classes
    
    % string
    %----------------------------------------------------------------------
    SEG(s).str = VOX.LEX(w,1).word;            % lexical string
    SEG(s).p   = P(:,s);                       % prior
    SEG(s).L   = L;                            % posteriors
    SEG(s).P   = Q(:,s);                       % prosody
    SEG(s).R   = R(:,s);                       % speaker
    SEG(s).I0  = VOX.I0;                       % first
    SEG(s).IT  = VOX.IT;                       % final
    
    disp({SEG.str})
    
end

% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,Q,R);

%% segmentation graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation'); clf;
spm_voice_segmentation(VOX.audio,SEG)


%% optmise spectral defaults with respect to classification accuracy
%==========================================================================
clear all
global VOX
PATH  = 'C:\Users\karl\Dropbox\Papers\Voice recognition\Sound files';
wtest = fullfile(PATH,'test','test.wav');
ttest = fullfile(PATH,'test','test.txt');

% reporting options
%--------------------------------------------------------------------------
VOX.graphics = 0;
VOX.mute     = 1;
VOX.onsets   = 0;

% expansion point (i.e., defaults)
%--------------------------------------------------------------------------
VOX.Nu  = 32;
VOX.Nv  = 8;
VOX.Tu  = 4;
VOX.Tv  = 1;
VOX.E   = 64;
VOX.F0  = 96;
VOX.F1  = 32;
VOX.F2  = 1024;

% search over variables
%--------------------------------------------------------------------------
Pu    = 16:4:64;
Pv    = 82:4:128;
for i = 1:numel(Pu)
    for j = 1:numel(Pv);
        try
            VOX = rmfield(VOX,'nu');
            VOX = rmfield(VOX,'nv');
        end
            
        VOX.F1 = Pu(i);
        VOX.F0 = Pv(j);
        
        [xY,word] = spm_voice_get_xY(PATH);
        P         = spm_voice_get_LEX(xY,word);
        A(i,j)    = spm_voice_test(wtest,ttest)
    end
end

imagesc(Pv,Pu,A), axis square, title('Accuracy','FontSize',16),drawnow


%% iterative inversion
%==========================================================================
VOX.analysis = 1;
VOX.graphics = 1;
VOX.interval = 0;
VOX.mute     = 0;
VOX.F2       = 1;

% get parameters for a particular word
%--------------------------------------------------------------------------
xY          = spm_voice_speak(14);
xY.P.amp    = 0;
xY.P.lat    = -8;
xY.P.ff1    = log(32);
xY.P.inf(1) = 1/100;

% compose and decompose iteratively
%--------------------------------------------------------------------------
for i = 1:8
    
    % generate timeseries and play
    %----------------------------------------------------------------------
    Y  = spm_voice_iff(xY);
    sound(Y,VOX.FS)
    
    % map back to parameters by inverting
    %----------------------------------------------------------------------
    VOX.F1   = exp(xY.P.ff1);
    xY       = spm_voice_ff(Y);
    xY.P.amp = 0;
    xY.P.lat = -8;
  
end



%% estimate formant frequency via maximum likelihood
%==========================================================================
global VOX
load VOX
load Ann
VOX.analysis = 1;
VOX.graphics = 1;
VOX.interval = 0;
VOX.onsets   = 1;
VOX.mute     = 1;
VOX.F2       = 1024;
VOX.FS       = 22050;
VOX.IT       = 1;

sound(Y,VOX.FS)
str       = {'yes','yes','yes','yes'};
w         = spm_voice_i(str);
[i,P]     = spm_voice_i(str);
P         = spm_softmax(log(P));
[L,F0,F1] = spm_voice_identity(Y,P(:,1));


VOX.F1       = F1;
VOX.F0       = F0;
VOX.mute     = 0;
VOX.analysis = 0;
VOX.graphics = 0;

spm_voice_read(Y,P);


% search over variables
%--------------------------------------------------------------------------
P     = spm_softmax(log(P)/2);

clear A
F1    = linspace(32,48,8);
F0    = linspace(200,300,4);
for i = 1:numel(F1)
    for j = 1:numel(F0);
        
        VOX.F1 = F1(i);
        VOX.F0 = F0(j);
        VOX.IT = 1;
        
        if 0
            % first word
            %--------------------------------------------------------------
            L      = spm_voice_get_word(Y,P(:,1));
            A(i,j) = log(L{1}(w(1)) + exp(-16))
        else

            % sentence
            %--------------------------------------------------------------
            SEG   = spm_voice_read(Y,P);
            L     = 0;
            for s = 1:numel(w)
                L = L + log(SEG(s).L{1}(w(s)) + exp(-16));
            end
            A(i,j) = L
        end
    end
end


% results over search
%--------------------------------------------------------------------------
imagesc(F0,F1,A), axis square, title('Accuracy','FontSize',16),drawnow
xlabel('Fundamental frequency'),ylabel(' formant frequency')

% segment using expected weaknesses
%--------------------------------------------------------------------------
A(:)   = spm_softmax(A(:));
F1     = F1*spm_vec(sum(A,2));
F0     = F0*spm_vec(sum(A,1));
VOX.F1 = F1;
VOX.F0 = F0;

VOX.mute     = 0;
VOX.analysis = 0;
VOX.graphics = 0;

spm_voice_read(Y);


%% auxiliary code
%==========================================================================

% load script and prompt for audio file: 32 words, at one word per second
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

