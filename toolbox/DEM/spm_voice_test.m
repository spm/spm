function [L] = spm_voice_test(wfile,sfile)
% Reads and translates a sound file to assess recognition accuracy
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_test.m 7558 2019-03-28 12:39:16Z karl $


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
G      = spm_conv(abs(read(wfile)),FS/4);
I      = find((diff(G(1:end - 1)) > 0) & (diff(G(2:end)) < 0));
[i,j]  = sort(G(I),'descend');
I      = sort(I(j(1:ns)));

% have orders been optimised with respect to the likelihood?
%--------------------------------------------------------------------------
OPT    = ~isfield(VOX,'nu');

%% run through sound file and evaluate likelihoods
%==========================================================================
xY    = {};
for s = 1:ns
    
    % retrieve epoch and decompose at fundamental frequency
    %----------------------------------------------------------------------
    clear xy
    Y     = read(wfile,round([-1/2 1/2]*FS + I(s)));
    j     = spm_voice_onsets(Y,FS);
    nj    = numel(j);
    for i = 1:nj
        xy(i,1) = spm_voice_ff(Y(j{i}),FS);
    end

    % store in xY
    %----------------------------------------------------------------------
    xY{s} = xy;
    
end

%% grid search to maximise classication accuracy
%==========================================================================
if OPT
    nu    = 4:size(xY{1}(1).Q,1);                              % order (Hz)
    nv    = 4:size(xY{1}(1).Q,2);                              % order (ms)
    LL    = zeros(numel(nu),numel(nv));
    for i = 1:numel(nu)
        for j = 1:numel(nv);
            VOX.nu = nu(i);
            VOX.nv = nv(j);
            
            % evaluate the log likelihood of correct word
            %--------------------------------------------------------------
            for s = 1:ns
                L       = spm_voice_likelihood(xY{s});     % log likelihoods
                L(:)    = spm_softmax(L(:));               % likelihoods
                L       = squeeze(sum(sum(L,2),3));        % marginal
                w       = strmatch(str{s,1},word,'exact'); % correct word
                LL(i,j) = LL(i,j) + log(L(w) + exp(-8));   % log likelihood
            end
            
        end
    end
    
    % optimal order of DCT basis functions
    %----------------------------------------------------------------------
    [i,j]  = find(LL == max(LL(:)),1);
    VOX.nu = nu(i);
    VOX.nv = nv(j);
else
    VOX.nu = 8;
    VOX.nv = 8;
end

%% illustrate classification accuracy
%==========================================================================
R     = [];                                    % lexical likelihood
for s = 1:ns
    
    % identify the most likely word and place in structure
    %----------------------------------------------------------------------
    [L,M]  = spm_voice_likelihood(xY{s});      % log likelihoods
    L(:)   = spm_softmax(L(:));                % likelihoods
    M      = spm_softmax(M);                   % likelihoods
    if ndims(M) > 2
        m  = squeeze(sum(sum(L,1),2));         % marginalise over lexicon
        M  = spm_dot(M,{m});                   % marginalise prosody
    end
    L      = squeeze(sum(sum(L,2),3));         % marginalise over sampling
    
    % identify the most likely word and prosody
    %----------------------------------------------------------------------
    [d,w]    = max(L);                         % most likely word
    [d,p]    = max(M);                         % most likely  prosodies
    R(:,s)   = L;                              % lexical likelihoods
    W(1,s)   = w(:);                           % lexical class
    P(:,s)   = p(:);                           % prosody classes
    str{s,2} = VOX.LEX(w,1).word;              % lexical string
    
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
word  = {VOX.LEX(:,1).word};
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
    c(:,w) = c(:,w) + log(R(:,s) + exp(-6));
    L      = L + log(R(w,s) + exp(-6));
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

if ~OPT, return, end

subplot(2,3,6), imagesc(nv,nu,LL), title('Accuracy','FontSize',16)
xlabel('order (intervals)'), ylabel('order (formants)'), axis square
drawnow

return


%% articulate: prosody without lexical content
%--------------------------------------------------------------------------
spm_voice_speak([],P,1,LEX,PRO,WHO);


%% articulate: with no prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,[],1,LEX,PRO,WHO);


%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
spm_voice_speak(W,P,[1;2],LEX,PRO,WHO);


