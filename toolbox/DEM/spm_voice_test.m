function spm_voice_test(wfile,sfile,LEX,PRO)
% Reads and translates a sound file to assess recognition accuracy
% FORMAT spm_voice_test(wfile,sfile,LEX,PRO)
%
% wfile  - .wav file
% sfile  - .txt file
% LEX    - lexical structure array
% PRO    - prodidy structure array
%
%  This routine tests, recognition on a small test corpus specified in
%  terms of a sound file and text file of successive words. It assesses
%  the accuracy of inference in relation to the known words and then plays
%  them back with and without prosody (or lexical content)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_test.m 7528 2019-02-06 19:19:49Z karl $


% create lexical structures for subsequent word recognition
%==========================================================================
str   = textread(sfile,'%s');                 % word string to test
ns    = numel(str);                           % number of words to test

% get sampling frequency (FS)
%--------------------------------------------------------------------------
try
    [Y,FS] = audioread(wfile,[1,128]);
    read   = @audioread;
catch
    [Y,FS] = wavread(wfile,[1,128]);
    read   = @wavread;
end

%% get data features from a wav file
%==========================================================================

% get fundamental frequency
%--------------------------------------------------------------------------
Y     = read(wfile);
I     = spm_voice_frequency(Y(abs(Y) > 1/32),FS);
F0    = FS/mean(I);

% identify the midpoints of ns words in terms of spectral maxima
%--------------------------------------------------------------------------
Y     = spm_conv(abs(Y),FS/4);
I     = find(diff(diff(Y,1) > 0) < 0);
[d,j] = sort(Y(I),'descend');
I     = sort(I(j(1:ns)));


%% run through sound file and evaluate likelihoods
%==========================================================================
R     = [];                                    % lexical likelihood
dI    = (-1:1)*FS/32;                          % sampling latencies
for s = 1:ns

    % retrieve epoch and decompose at fundamental frequency
    %----------------------------------------------------------------------
    for i = 1:numel(dI)
        Y     = read(wfile,round([-1/2 1/2]*FS + I(s) + dI(i)));
        xy(i) = spm_voice_ff(Y,FS,F0);
    end
    
    % identify the most likely word and place in structure
    %----------------------------------------------------------------------
    [L,M]  = spm_voice_likelihood(xy,LEX,PRO); % log likelihoods
    L(:)   = spm_softmax(L(:));                % likelihoods
    M      = spm_softmax(M);                   % likelihoods
    m      = squeeze(sum(sum(L,1),2));         % marginalise over lexicon
    L      = squeeze(sum(sum(L,2),3));         % marginalise over sampling
    M      = spm_dot(M,{m});                   % marginalise prosody
    
    % identify the most likely word and prosody
    %----------------------------------------------------------------------
    [d,w]    = max(L);                         % most likely word
    [d,p]    = max(M);                         % most likely  prosodies
    R(:,s)   = L;                              % lexical likelihoods
    W(1,s)   = w(:);                           % lexical class
    P(:,s)   = p(:);                           % prosody classes
    str{s,2} = LEX(w,1).word;                  % lexical string

end


%% illustrate classification accuracy
%==========================================================================

% display true and inferred strings
%--------------------------------------------------------------------------
disp(str)

% assess classification accuracy
%--------------------------------------------------------------------------
clear q p
word  = {LEX(:,1).word};
nw    = length(word);
c     = zeros(nw,nw);
q     = zeros(nw,ns);
p     = zeros(nw,ns);
for s = 1:ns
    w      = strmatch(str{s,2},word,'exact');
    q(w,s) = 1;
    w      = strmatch(str{s,1},word,'exact');
    p(w,s) = 1;
    c(:,w) = c(:,w) + log(R(:,s) + exp(-6));
end

% graphics
%--------------------------------------------------------------------------
a  = (sum(p(:) ~= q(:))/2)/ns;
subplot(4,1,1), imagesc(1 - p), title('True','FontSize',16),
set(gca,'YTick',1:nw,'YTickLabel',word)
subplot(4,1,2), imagesc(1 - R), title('Inferred','FontSize',16)
set(gca,'YTick',1:nw,'YTickLabel',word)
xlabel('word')
subplot(2,2,3), imagesc(spm_softmax(c))
set(gca,'XTick',1:nw,'XTickLabel',word)
set(gca,'YTick',1:nw,'YTickLabel',word)
titlestr = sprintf('Classification accuracy %-2.0f p.c.',100*(1 - a));
title(titlestr,'FontSize',16), axis square
subplot(2,2,4), imagesc(P)
title('Prosody','FontSize',16), axis square
set(gca,'YTick',1:size(P,2))
xlabel('word'), ylabel('mode')


%% articulate: prosody without lexical content
%--------------------------------------------------------------------------
DT     = diff([1;I]/FS) - 2/3;
spm_voice_speak([],P,LEX,PRO,DT);

%% articulate: with no prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,[],LEX,PRO,DT);


%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,P,LEX,PRO,DT);


