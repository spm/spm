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
% $Id: spm_voice_get_word.m 7546 2019-03-18 11:02:22Z karl $

%% get  parameters from VOX
%==========================================================================
global VOX
VOX.onsets = 1;

I0  = VOX.I0;
IT  = VOX.IT;


% get source (recorder) and FS
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    FS    = get(wfile,'SampleRate');
    read  = @getaudiodata;
    
    % ensure 2 second of data has been accumulated
    %----------------------------------------------------------------------
    dt    = (get(wfile,'TotalSamples') - I0)/FS;
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
for i = 1:2
    
    % check for content
    %----------------------------------------------------------------------
    Y = read(wfile);
    n = numel(Y);
    j = round(0:2*FS) + I0;
    G = spm_voice_filter(Y(j(j < n)),FS);
    G = spm_conv(G,FS/6);
    I = find((diff(G(1:end - 1)) > 0) & (diff(G(2:end)) < 0));
    I = I(I > IT(2) & G(I) > 1/128);
    
    % advance pointer if silence
    %----------------------------------------------------------------------
    if isempty(I)
        
        % advance pointer one second
        %------------------------------------------------------------------
        I0 = I0 + FS;
        
        % ensure 2 second of data has been accumulated
        %------------------------------------------------------------------
        if isa(wfile,'audiorecorder')
            dt = (get(wfile,'TotalSamples') - I0)/FS;
            pause(2 - dt);
        end
        
    else
        
        % move to centre of next word
        %------------------------------------------------------------------
        I0 = I0 + I(1);
        break
    end
end

% break if EOF
%--------------------------------------------------------------------------
if isempty(I)
    O       = {};
    VOX.I0  = I0;
    return
end

% retrieve epoch and decompose at fundamental frequency
%--------------------------------------------------------------------------
for i = 1:size(VOX.I,1)
    j       = VOX.I(i,1)*FS:VOX.I(i,2)*FS;
    j       = round(j) + I0;
    xy(i,1) = spm_voice_ff(Y(j(j < n)),FS);
end


% identify the most likely word
%==========================================================================
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

% expected sampling interval
%--------------------------------------------------------------------------
IT      = m'*VOX.I*FS;


% advance pointer based on marginal over samples
%--------------------------------------------------------------------------
[m,i]   = max(m);
O{1}    = L;
O{2}    = M;
O{3}    = N;
VOX.I0  = I0;
VOX.IT  = IT;

return

