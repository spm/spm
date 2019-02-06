function spm_voice_read(LEX,PRO)
% Reads and translates a wav file
% FORMAT [S] = DEM_birdsong(file)
%
% wfile  - .wav file
% LEX    - lexical dictionary array
% ns     -  number of word samples
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_read.m 7528 2019-02-06 19:19:49Z karl $


% create lexical structures for subsequent word recognition
%==========================================================================
wfile = '../test.wav';

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

%% get data features from a succession (1 Hz) exemplars
%--------------------------------------------------------------------------
r     = [];
I     = 1;
dI    = (-1:1)*FS/16;
for s = 1:ns
    % find next initial time point
    %----------------------------------------------------------------------
    c = 1;
    while c
        try
            Y = read(wfile,round([0 1/2]*FS + I(s)));
        catch
            break
        end
        [y,i] = max(abs(Y));
        if y > 1/16 && i > FS/4
            I(s) = I(s);
            c    = 0;
        else
            I(s) = I(s) + FS/8;
        end
    end
    
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
    
    % identify the most likely word and place in structure
    %----------------------------------------------------------------------
    %     DT       = exp(xy(i).P(:,2)) + exp(xy(i).P(:,3)) + dI(i)/FS;
    %     I(s + 1) = I(s) + FS*DT;
    
    
    % if ~strcmp(STR{s},str{s}), break, end
    %     spm_voice_iff(xy(i),FS,1/32);
    %     disp(str{s}),pause(1)
    
end


%% articulate: with lexical content and prosody
%--------------------------------------------------------------------------
DT     = diff([1;I]/FS) - 2/3;
spm_voice_speak(W,P,LEX,PRO,DT);

