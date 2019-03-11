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
% requires the following in the global variable voice_options:
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
% $Id: spm_voice_get_word.m 7540 2019-03-11 10:44:51Z karl $

%% get  parameters from voice_options
%==========================================================================
global voice_options
voice_options.onsets = 1;

LEX = voice_options.LEX;
PRO = voice_options.PRO;
WHO = voice_options.WHO;
I0  = voice_options.I0;


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
Y     = read(wfile);

% find next word
%--------------------------------------------------------------------------
DI    = FS/4;
IS    = numel(Y);
I     = [];
for i = 1:4
    
    % check for content
    %----------------------------------------------------------------------
    j = (0:FS) + I0;
    try
        I = spm_voice_check(Y(j(j < IS)),FS);
    end
    
    % advance pointer if silence
    %----------------------------------------------------------------------
    if isempty(I) && (I0 + DI + FS) < IS
        I0 = I0 + DI;
    else
        break
    end
end

% break if EOF
%--------------------------------------------------------------------------
if isempty(I)
    O = {};
    return
end

% get sampling epoch
%--------------------------------------------------------------------------
I0     = I0 + I(1);
j      = (0:FS) + I0;
I      = spm_voice_onset(Y(j(j < IS)),FS);
I0     = I0 + I(1);
j      = (0:FS) + I0;
[I,dI] = spm_voice_onset(Y(j(j < IS)),FS);
I0     = I0 + I(1);
dI     = [dI;I(end)] - I(1);


% retrieve epoch and decompose at fundamental frequency
%--------------------------------------------------------------------------
clear xy
for i = 1:numel(dI)
    y     = Y(round((0:dI(i)) + I0));
    xy(i) = spm_voice_ff(y,FS);
end

% identify the most likely word
%==========================================================================
[L,M,N] = spm_voice_likelihood(xy,LEX,PRO,WHO);
L(:)    = spm_softmax(L(:));                % likelihoods
M       = spm_softmax(M);                   % likelihoods
N       = spm_softmax(N);                   % likelihoods
m       = squeeze(sum(sum(L,1),2));         % marginalise over lexicon
L       = squeeze(sum(sum(L,2),3));         % marginalise over sampling
if numel(dI) > 1
    M   = spm_dot(M,{m});                   % marginalise prosody
    N   = spm_dot(N,{m});                   % marginalise identity
end

% advance pointer based on marginal over samples
%--------------------------------------------------------------------------
[m,i] = max(m);
O{1}  = L;
O{2}  = M;
O{3}  = N;

voice_options.I0 = I0 + xy(i).i(end);

return

