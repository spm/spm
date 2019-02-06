function [LEX,PRO] = spm_voice(path)
% Creates lexical and prosody cell arrays from sound file exemplars
% FORMAT [LEX,PRO] = spm_voice(path)
%
% path     -  directory containing lists of exemplar words
%
% LEX(w,k) -  structure array for k variants of w words
% PRO(p)   -  structure array for p aspects of prosody
%
%  This routine creates a pair of structure array is used to infer the
%  lexical and prosody of a word. It uses a library of sound files, each
%  containing 32 words spoken with varying prosody. The name of the sound
%  file labels the word in question. These exemplars are then transformed
%  (using a series of discrete cosine and Hilbert transforms) into a set of
%  parameters, which summarise the lexical content and prosody. The inverse
%  transform generates  timeseries that can be played to articulate a word.
%  The transform operates on a word structure xY to create lexical and
%  prosody parameters (Q and P respectively).
%  The accuracy of  lexical inference (i.e., voice the word recognition)
%  is assessed  using the exemplar (training) set and a narrative sound
%  file called '../test.wav'  (and associated '../test.txt').
%  The operation of each subroutine can be examined using graphical
%  outputs by selecting the appropriate options in a voice recognition
%  specific global variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice.m 7528 2019-02-06 19:19:49Z karl $



%% setup options and files
%==========================================================================

% directory of sound files if necessary
%--------------------------------------------------------------------------
if ~nargin
    PATH = 'C:\Users\karl\Dropbox\Papers\Voice recognition\Sound files';
end

global voice_options
voice_options.graphics = 1;
voice_options.mute     = 0;
voice_options.onsets   = 0;

%% get the list of words and sampling frequency (FS)
%--------------------------------------------------------------------------
spm_figure('GetWin','voice')
cd(PATH)
wfile      = dir('*.wav');

try
    [Y,FS] = audioread(wfile(1).name,[1,128]);
    read   = @audioread;
catch
    [Y,FS] = wavread(wfile(1).name,[1,128]);
    read   = @wavread;
end


% assemble cell array of word structures for subsequent characterisation
%==========================================================================
nw    = numel(wfile);                               %  number of words
ns    = 32;                                         %  number of samples
for w = 1:nw
    
    % get lexicon name and create structure
    %----------------------------------------------------------------------
    wname         = wfile(w).name;                  %  name of word
    [d,word]      = fileparts(wname);
    LEX(w,1).word = word;
    
    % get fundamental frequency (F0)
    %----------------------------------------------------------------------
    Y     = read(wname);
    Ymax  = max(abs(Y));
    I     = spm_voice_frequency(Y(abs(Y) > Ymax/32),FS);
    F0    = FS/mean(I);
    
    % get the midpoint of words from the (maxima) of successive exemplars
    %----------------------------------------------------------------------
    Y     = spm_conv(abs(Y),FS/4);
    I     = find(diff(diff(Y,1) > 0) < 0);
    [i,j] = sort(Y(I),'descend');
    I     = sort(I(j(1:ns)));
    
    % and plot
    %----------------------------------------------------------------------
    subplot(2,1,1), plot(Y), hold on, plot(I,Y(I),'ro'), hold off
    title(sprintf('Power peaks - %s',word),'FontSize',16)
    xlabel('frequency (hertz)'), ylabel('power'), box off
    
    for s = 1:ns
        
        % retrieve (one second) epoch around midpoint and transform
        %------------------------------------------------------------------
        Y        = read(wname,[-500 500]*FS/1000 + I(s));
        xY(w,s)  = spm_voice_ff(Y/Ymax,FS,F0);
        
        %  fix amplitude for subsequent analysis of prosody
        %------------------------------------------------------------------
        xY(w,s).P.amp = log(1);
        
        %  apply inverse transform and play, if requested
        %------------------------------------------------------------------
        if ~voice_options.mute
            spm_voice_iff(xY(w,s),FS,1/4);
        end
    end
end

%% assemble parameters for subsequent recognition (inversion) analysis
%==========================================================================
clear P Q
for w = 1:size(xY,1)
    for s = 1:size(xY,2)
        Q{w}(:,s) = spm_vec(xY(w,s).Q);
        P{w}(:,s) = spm_vec(xY(w,s).P);
    end
end

% concatenate and illustrate distribution of prosody parameters
%==========================================================================
P     = full(spm_cat(P)');
str   = {'amplitude','duration','f0','f1','timbre','inf 1','inf 2','inf 3'};
for i = 1:size(P,2)
    subplot(3,3,i)
    hist(exp(P(:,i)),32), axis square
    title(sprintf('%s mean: %.2f',str{i},mean(exp(P(:,i)))))
end



%% Eigenmodes of lexical (U) and prosody (V) parameters
%==========================================================================

% lexical
%--------------------------------------------------------------------------
[U,S] = spm_svd(cov(spm_cat(Q)'));
S     = log(diag(S));
s     = find(S > (max(S) - 4),1,'last');
subplot(2,2,1), bar(S)
title('Eigenvalues - lexical','FontSize',16), ylabel('log value')
hold on, plot([s,s],[0,max(S)],'r:'), hold off
xlabel(sprintf('eigenbasis (%i)',s)), axis square

% prosody (based on correlation)
%--------------------------------------------------------------------------
[V,S] = spm_svd(corr(P));
S     = log(diag(S));
subplot(2,2,2), bar(S)
title('Eigenvalues - prosidy','FontSize',16)
xlabel('eigenbasis'), ylabel('log value'), axis square

subplot(2,1,2), bar(V')
title('Eigenmodes - prosidy','FontSize',16)
xlabel('prosody mode'), ylabel('amplitude'), axis square
str   = {'amplitude','duration','frequency','formant','timbre','inf 1','inf 2','inf 3'};
legend(str)


%% prosody: use K means clustering of eigenmodes for prosody categories
%==========================================================================

% save mean and eigenmodes of prosody profile in prosody array
%--------------------------------------------------------------------------
p0       = mean(P);
PRO(1).P = spm_unvec(p0(:),xY(1).P);
PRO(1).C = cov(P);
PRO(1).U = full(V);

% categorise  principal eigenmodes into three levels
%--------------------------------------------------------------------------
K     = 2;
for p = 1:K
    
    % save lexical exemplars
    %----------------------------------------------------------------------
    W          = bsxfun(@minus,P,p0)*PRO(1).U(:,p);
    [pL,pE,pC] = spm_kmeans(W,3,'fixed-points');
    [d,i]      = sort(pL,'descend');
    [d,j]      = max(abs(V(:,p)));
    pP         = 1./squeeze(pC);
    pL         = log(pL + exp(-8));
    
    
    % save lexical exemplars
    %----------------------------------------------------------------------
    PRO(p).str = str(j);
    PRO(p).pE  = pE(i);
    PRO(p).pP  = pC(i);
    PRO(p).pP  = pP(i);
    PRO(p).pL  = pL(i);
    
end


%% lexical: use K means clustering of eigenmodes for prosody categories
%==========================================================================

%  evaluate word specific means and covariances in eigenmode space
%--------------------------------------------------------------------------
LEX      = LEX(:,1);                              % lexical structure array
N        = 1;                                     % number of variants
K        = 64;                                    % number of eigenmodes
p0       = mean(spm_cat(Q),2);                    % average word
LEX(1).Q = spm_unvec(p0,xY(1).Q);
LEX(1).C = cov(spm_cat(Q)');
LEX(1).U = U(:,1:K);

for w = 1:nw
    W          = bsxfun(@minus,Q{w},p0)'*LEX(1).U;
    [pL,pE,pC] = spm_kmeans(W,N,'fixed-points');
    [i,k]      = sort(pL,'descend');
    
    for i = k
        LEX(w,i).pE = pE(i,:)';
        LEX(w,i).pC = pC(:,:,i);
        LEX(w,i).pP = spm_inv(pC(:,:,i));
    end
end

% covariance normalisation; based on within and between word variance
%--------------------------------------------------------------------------
SSQ   = trace(LEX(1).U'*cov(spm_cat(Q)')*LEX(1).U);
SSP   = trace(cov(spm_cat({LEX.pE})));
SSR   = (SSQ - SSP);
for w = 1:nw
    for k = 1:N
        
        pC          = LEX(w,k).pC;
        pC          = pC*K/trace(pC) + eye(K,K)/4;
        pC          = SSR*pC/trace(pC);
        pP          = spm_inv(pC);
        
        LEX(w,k).pC = pC;
        LEX(w,k).pP = pP;
    end
end

%%  illustrate voice recognition (i.e., inference) based upon eigenmodes
%==========================================================================


% likelihood of training set
%--------------------------------------------------------------------------
q     = [];
p     = [];
r     = [];
J     = 0;
for i = 1:size(xY,1)
    for j = 1:size(xY,2)
        
        % evaluate lexical (L) and prosody (M) likelihoods
        %------------------------------------------------------------------
        [L,M] = spm_voice_likelihood(xY(i,j),LEX,PRO);
        L(:)  = spm_softmax(L(:));
        M     = spm_softmax(M);
        L     = sum(L,2);
        
        % inferred and true lexical outcome
        %------------------------------------------------------------------
        q(:,end + 1) = L;
        p(i,end + 1) = 1;
        r        = [r,M];
        
        % joint distibution of lexical and prosody outcomes
        %------------------------------------------------------------------
        J  = J + spm_cross({L},num2cell(M,1));
        
    end
end

% illustrate classification accuracy
%--------------------------------------------------------------------------
PRO(1).J = J/sum(J(:));

% show results
%--------------------------------------------------------------------------
subplot(2,1,1); imagesc(q)
a     = 1 - sum((p(:) - q(:)) > 1/8)/sum(p(:));
str   = sprintf('Lexical classification accuracy %-2.0f p.c.',100*a);
title(str,'FontSize',16), xlabel('exemplars'), ylabel('words')

subplot(4,1,3); imagesc(r)
pr    = sum(r,2);
pr    = 100*pr/sum(pr);
str   = sprintf('Prosody proportions: p.c.: %.0f %.0f %.0f ',pr);
xlabel('exemplars'), ylabel('prosody'), title(str,'FontSize',16)

%% articulate every word under all combinations of (3 levels) of prosody
%--------------------------------------------------------------------------
for w = 1:nw
    for i = 1:3
        for j = 1:3
            spm_voice_speak(w,[i;j],LEX,PRO,1/4);
        end
    end
end


%%  apply to test narrative of 87 words
%--------------------------------------------------------------------------
spm_voice_test('../test.wav','../test.txt',LEX,PRO);


%% save lexical and prosody arrays in sound file directory
%--------------------------------------------------------------------------
save LEX LEX PRO

return



%% auxiliary code
%==========================================================================

% load script and prompt for audio file: 32 words, at one word per second
%--------------------------------------------------------------------------
tfile = '../test.txt';
str   = textread(tfile,'%s');
str   = unique(str);
for s = 1:numel(str)
    clc,  disp([str(s) 'new file']), pause
    for i = 1:32
        clc,disp(str(s));     pause(0.9)
        clc,disp('    -*- '); pause(0.1)
    end
end


%% check formant quality following gigenreduction
%==========================================================================
for w = 1:nw
    for s = 1:8
        
        % play each word after eigenreduction
        %----------------------------------------------------------------
        xy      = xY(w,s);
        xQ      = spm_vec(xy.Q);
        xQ      = xQ - mean(xQ);
        xQ      = spm_vec(LEX(1).Q) + LEX(1).U *LEX(1).U'*xQ
        xy.Q    = spm_unvec(xQ,xy.Q);
        spm_voice_iff(xy,FS,1/2);
    end
end

