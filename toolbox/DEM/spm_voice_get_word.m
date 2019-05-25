function [O] = spm_voice_get_word(wfile,P)
% Evaluates the likelihood of the next word in a file or object
% FORMAT [O] = spm_voice_get_word(wfile,P)
%
% wfile  - .wav file, audiorecorder object or (double) time series
% P      - lexical prior [optional]
%
% O{1}   - lexical likelihood (or posterior if priors are specified)
% O{2}   - prosody likelihood
% O{3}   - speaker likelihood
%
% requires the following in the global variable VOX:
%
% LEX    - lexical structure array
% PRO    - prodidy structure array
% WHO    - speaker structure array
% FS     - sampling    frequency (Hz)
% F0     - fundamental frequency (Hz)
% IT     - index or pointer to offset of last word (i.e., CurrentSample)
%
% and updates:
% I0     - index or pointer to onset  of current word
% IT     - index or pointer to offset of current word
%
% This routine evaluates the likelihood of a word, prosody and identity by
% inverting successive epochs of data from an audiofile or device starting
% at VOX.IT.  Based on the word with the least variational free energy , it
% updates the index, ready for the next word. Priors over the words can be
% specified to implement an Occam's window (of three nats); thereby
% restricting the number of lexical entries evaluated - and augmenting the
% likelihoods to give the posterior probability over words.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_word.m 7598 2019-05-25 13:09:47Z karl $


%% log prior over lexical content 
%==========================================================================
global VOX
if nargin < 2
    nw = numel(VOX.LEX);
    LP = ones(nw,1)/nw;
else
    LP = log(P + exp(-8));
end

% within Ockham's window W
%--------------------------------------------------------------------------
W      = find(LP > (max(LP) - 3));

%% get onset of next word 
%==========================================================================
[Y,I,FS] = spm_voice_get_next(wfile);

% break if EOF
%--------------------------------------------------------------------------
if isempty(I), O  = {}; return, end


%% get 1 second segment and remove previous word
%==========================================================================

% get intervals (j) for this peak
%--------------------------------------------------------------------------
n    = numel(Y);
j    = fix((0:FS + FS/4) + I);
y    = Y(j(j < n & j > 1));
j    = logical(j < VOX.IT);
y(j) = 0;
j    = spm_voice_onsets(y,FS);

% retrieve epochs and decompose at fundamental frequency
%--------------------------------------------------------------------------
clear xy J
for i = 1:numel(j)
    xy(i,1)       = spm_voice_ff(y(j{i}),FS);
    xy(i,1).P.lat = log((I + j{i}(1) - VOX.IT)/FS);
    J(i,:)        = I + [j{i}(1),j{i}(end)];
end

% select the most likely action (interval)
%==========================================================================
[L,M,N] = spm_voice_likelihood(xy,W);
L       = bsxfun(@plus,L,LP);
Q       = spm_softmax(L);
F       = sum(Q.*(L - log(Q + exp(-16))));
D       = (J(:,2) - J(:,1))/FS;
F       = F(:) - log(D);
[n,m]   = max(F);

% posteriors
%--------------------------------------------------------------------------
O{1}    = spm_softmax(L(:,m));              % posteriors
O{2}    = spm_softmax(M(:,:,m));            % likelihood
O{3}    = spm_softmax(N(:,:,m));            % likelihood

% update pointers based upon selected word
%--------------------------------------------------------------------------
VOX.I0  = fix(J(m,1));
VOX.IT  = fix(J(m,2));
VOX.J   = j;
VOX.I   = I;

return


