function [L] = spm_voice_test(wfile,sfile)
% Read and translate a sound file to assess recognition accuracy
% FORMAT [L] = spm_voice_test(wfile,sfile)
%
% wfile   - .wav file
% sfile   - .txt file
%
% rqeuires
% VOX.LEX - lexical structure array
% VOX.PRO - prodidy structure array
% VOX.WHO - speaker structure array
%
% L       - accuracy (log likelihood)
%
%  This routine tests, recognition on a small test corpus specified in
%  terms of a sound file and text file of successive words. It assesses
%  the accuracy of inference in relation to the known words and then plays
%  them back with and without prosody (or lexical content)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_test.m 7750 2019-12-05 17:54:29Z spm $


% create lexical structures for subsequent word recognition
%==========================================================================
global VOX
str   = textread(sfile,'%s');                 % word string to test
word  = {VOX.LEX(:,1).word};                  % words in lexicon
ns    = numel(str);                           % number of words to test

% get sampling frequency (FS)
%--------------------------------------------------------------------------
try
    [Y,FS] = audioread(wfile,[1,1]);
    read   = @audioread;
catch
    [Y,FS] = wavread(wfile,[1,1]);
    read   = @wavread;
end

%% get data features from a wav file
%==========================================================================

% get F0 and the midpoint of words (maxima of acoutics power)
%--------------------------------------------------------------------------
G      = spm_voice_check(read(wfile),FS,1/4);
I      = find((diff(G(1:end - 1)) > 0) & (diff(G(2:end)) < 0));
[d,j]  = sort(G(I),'descend');
I      = sort(I(j(1:ns)));


%% run through sound file and evaluate likelihoods
%==========================================================================
xY    = {};
for s = 1:ns
    
    % retrieve epoch and decompose at fundamental frequency
    %----------------------------------------------------------------------
    Y      = read(wfile,round([-1/2 1/2]*FS + I(s)));
    j      = spm_voice_onsets(Y,FS);
    xY{s}  = spm_voice_ff(Y(j{end}),FS);

end


%% illustrate classification accuracy
%==========================================================================
R     = [];                                    % lexical likelihood
for s = 1:ns
    
    % identify the most likely word and place in structure
    %----------------------------------------------------------------------
    [L,M]    = spm_voice_likelihood(xY{s});    % log likelihoods

    % posteriors and pointer
    %--------------------------------------------------------------------------
    L        = spm_softmax(L);                 % posteriors
    M        = spm_softmax(M);                 % likelihood
    
    % identify the most likely word and prosody
    %----------------------------------------------------------------------
    [d,w]    = max(L);                         % most likely word
    [d,p]    = max(M);                         % most likely prosodies
    R(:,s)   = L;                              % lexical likelihoods
    W(1,s)   = w(:);                           % lexical class
    P(:,s)   = p(:);                           % prosody classes
    str{s,2} = VOX.LEX(w).word;                % lexical string
    
end


%% classification accuracy
%==========================================================================
spm_figure('GetWin','Accuracy (test)'); clf

% display true and inferred strings
%--------------------------------------------------------------------------
disp(str)

% assess classification accuracy
%--------------------------------------------------------------------------
clear q p
word  = {VOX.LEX.word};
nw    = length(word);
c     = zeros(nw,nw);
q     = zeros(nw,ns);
p     = zeros(nw,ns);
L     = 0;
for s = 1:ns
    w      = strmatch(str{s,2},word,'exact');
    q(w,s) = 1;
    w      = strmatch(str{s,1},word,'exact');
    p(w,s) = 1;
    c(:,w) = c(:,w) + log(R(:,s) + exp(-16));
    L      = L + log(R(w,s) + exp(-16));
end

% graphics
%--------------------------------------------------------------------------
a  = (sum(p(:) ~= q(:))/2)/ns;
a  = ceil(100*(1 - a));
subplot(4,1,1), imagesc(1 - p), title('True','FontSize',16),
set(gca,'YTick',1:nw,'YTickLabel',word)
subplot(4,1,2), imagesc(1 - R), title('Inferred','FontSize',16)
set(gca,'YTick',1:nw,'YTickLabel',word)
xlabel('word')
subplot(2,3,4), imagesc(spm_softmax(c))
set(gca,'XTick',1:nw)
set(gca,'YTick',1:nw,'YTickLabel',word)
titlestr = sprintf('%-2.0f p.c. accuracy: %-2.1f',a,L);
title(titlestr,'FontSize',16), axis square
subplot(2,3,5), imagesc(P)
title('Prosody','FontSize',16)
set(gca,'YTick',1:size(P,1)),set(gca,'YTickLabel',{VOX.PRO.str})
xlabel('word'), ylabel('mode'), axis square

if exist('nu','var')
    subplot(2,3,6), imagesc(nv,nu,LL), title('Accuracy','FontSize',16)
    xlabel('order (intervals)'), ylabel('order (formants)'), axis square
    colorbar
end
drawnow
