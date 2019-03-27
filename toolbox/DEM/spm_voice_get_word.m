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
% $Id: spm_voice_get_word.m 7557 2019-03-27 17:11:16Z karl $

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
    pause(1 - dt);

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
    G = spm_conv(abs(Y(j(j < n))),FS/4);
    I = find((diff(G(1:end - 1)) > 0) & (diff(G(2:end)) < 0));
    I = I(G(I) > 1/128);
    
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
            pause(1 - dt);
        end
        
    else
        
        % move pointer to 300ms before peak
        %------------------------------------------------------------------
        I  = VOX.IT + I(1) - FS/2;
        break
    end
end

% break if EOF
%--------------------------------------------------------------------------
if isempty(I)
    O  = {};
    return
end

% get 1 second segment and remove previous word
%--------------------------------------------------------------------------
j    = fix((0:FS) + I);
y    = Y(j(j < n));
j    = logical(j < VOX.IT);
y(j) = 1/1024;

% retrieve epoch and decompose at fundamental frequency
%--------------------------------------------------------------------------
clear xy
j     = spm_voice_onsets(y,FS);
nj    = numel(j);
for i = 1:nj
    xy(i,1) = spm_voice_ff(y(j{i}),FS);
    J(i,:)  = I + [j{i}(1),j{i}(end)];
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

% advance pointer based on marginal over samples
%--------------------------------------------------------------------------
I       = m'*J;
O{1}    = L;
O{2}    = M;
O{3}    = N;
VOX.I0  = fix(I(1));
VOX.IT  = fix(I(2));

return

