function [O] = spm_voice_get_word(wfile,P)
% Retrieves likelihoods from an audio device or file
% FORMAT [O] = spm_voice_get_word(wfile,P)
%
% wfile  - .wav file or audiorecorder object
% P      - lexical  prior [optional]
%
% O{1}   - lexical likelihood (or posterior if priors are specified)
% O{2}   - prosody likelihood
% O{3}   - speaker likelihood
%
% requires the following in the global variable VOX:
% LEX    - lexical structure array
% PRO    - prodidy structure array
% WHO    - speaker structure array
% FS     - sampling    frequency (Hz)
% F0     - fundamental frequency (Hz)
% I0     - index or pointer to current epoch (i.e., CurrentSample)
%
% This routine evaluates the likelihood of a word, prosody and identity by
% inverting successive epochs of data from an audiofile or device
% starting at the index.  Based on the word with the greatest
% posterior, it updates the index, ready for the next word. Priors over the
% words can be specified to implement an Occam's window; thereby
% restricting the number of lexical entries evaluated - and augmenting the
% likelihoods to return the posterior probability over words.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_word.m 7561 2019-03-30 10:39:07Z karl $

%% get  parameters from VOX
%==========================================================================
global VOX
VOX.onsets = 1;

% get source (recorder) and FS
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    FS    = get(wfile,'SampleRate');
    read  = @getaudiodata;
    
    % ensure 2 second of data has been accumulated
    %----------------------------------------------------------------------
    dt    = (get(wfile,'TotalSamples') - VOX.IT)/FS;
    pause(2 - dt);
    
else
    
    % sound file (read)
    %----------------------------------------------------------------------
    try
        xI     = audioinfo(wfile);
        FS     = xI.SampleRate;
        read   = @audioread;
    catch
        [Y,FS] = wavread(wfile,[1 1]);
        read   = @wavread;
    end
end


%% find word and evaluate likelihoods
%==========================================================================

% find next word
%--------------------------------------------------------------------------
for i = 1:4
    
    % find next spectral peak (I)
    %----------------------------------------------------------------------
    Y = read(wfile);
    n = numel(Y);
    j = fix((0:FS) + VOX.IT);
    G = spm_conv(abs(Y(j(j < n))),FS*VOX.C);
    I = find((diff(G(1:end - 1)) > 0) & (diff(G(2:end)) < 0));
    I = I(G(I) > (VOX.U + min(G)));
    
    
    % advance pointer if silence
    %----------------------------------------------------------------------
    if isempty(I)
        
        % advance pointer 500 ms
        %------------------------------------------------------------------
        VOX.IT = VOX.IT + FS/2;
        
        % ensure 2 second of data has been accumulated
        %------------------------------------------------------------------
        if isa(wfile,'audiorecorder')
            dt = (get(wfile,'TotalSamples') - VOX.IT)/FS;
            pause(2 - dt);
        end
        
    else
        
        % move pointer to 500ms before peak
        %------------------------------------------------------------------
        I  = VOX.IT + I - FS/2;
        I  = I(I < (I(1) + FS/2));
        I  = max(I,1);
        break
    end
end

% break if EOF
%--------------------------------------------------------------------------
if isempty(I)
    O  = {};
    return
end

%% get 1 second segment and remove previous word
%--------------------------------------------------------------------------
for w = 1:numel(I)
    
    % get intervals (j) for this peak
    %----------------------------------------------------------------------
    for i = 1:8
        try
            j    = fix((0:FS) + I(w));
            y    = Y(j(j < n));
            j    = logical(j < VOX.IT);
            y(j) = 1/512;
            j    = spm_voice_onsets(y,FS);
            break
        catch
            I(w) + I(w) + FS/64;
        end
    end
    
    % retrieve epochs and decompose at fundamental frequency
    %----------------------------------------------------------------------
    nj    = numel(j);
    J     = zeros(2,0);
    xy    = struct('Y',[],'Q',[],'P',struct([]),'i',[]);
    for i = 1:nj
        xy(i,1) = spm_voice_ff(y(j{i}),FS);
        J(:,i)  = I(w) + [j{i}(1);j{i}(end)];
    end
    
    
    % identify the most likely word
    %======================================================================
    [L,M,N] = spm_voice_likelihood(xy);
    L(:)    = spm_softmax(L(:));                % likelihoods
    M       = spm_softmax(M);                   % likelihoods
    N       = spm_softmax(N);                   % likelihoods
    m       = squeeze(sum(sum(L,1),2));         % marginalise over lexicon
    L       = squeeze(sum(sum(L,2),3));         % marginalise over sampling
    if size(M,3) > 1
        M   = spm_dot(M,{m});                   % marginalise prosody
        N   = spm_dot(N,{m});                   % marginalise identity
    end
    
    % posteriors and pointer averaged under lexical marginal
    %----------------------------------------------------------------------
    O{1,w} = L;
    O{2,w} = M;
    O{3,w} = N;
    q(:,w) = J*m;
    p(:,w) = L;
    
end

%% Ambiguity
%--------------------------------------------------------------------------
[d,m]  = max(sum(p.*log(p)));


%% return posteriors and advance pointer
%--------------------------------------------------------------------------
O      = O(:,m);
VOX.I0 = fix(q(1,m));
VOX.IT = fix(q(2,m));

return


%% identify the most likely word
%==========================================================================
% [L,M,N] = spm_voice_likelihood(xy);
% L       = spm_softmax(L);                   % likelihoods
% M       = spm_softmax(M);                   % likelihoods
% N       = spm_softmax(N);                   % likelihoods
% F       = sum(L.*log(L));
% [F,m]   = max(F(1,1,:));
%
% L       = L(:,1,m);                         % marginalise over sampling
% M       = M(:,:,m);                         % marginalise prosody
% N       = N(:,:,m);                         % marginalise identity
% I       = J(m,:);
%
% % advance pointer based on marginal over samples
% %--------------------------------------------------------------------------
% O{1}    = L;
% O{2}    = M;
% O{3}    = N;
% VOX.I0  = fix(I(1));
% VOX.IT  = fix(I(2));


