function [O,F] = spm_voice_get_word(wfile,P)
% Evaluates the likelihood of the next word in a file or object
% FORMAT [O,F] = spm_voice_get_word(wfile,P)
%
% wfile  - .wav file, audiorecorder object or (double) time series
% P      - lexical prior probability [optional]
%
% O{1}   - lexical likelihood (or posterior if priors are specified)
% O{2}   - prosody likelihood
% O{3}   - speaker likelihood
%
% F      - maximum free energy (i.e., log evidence)
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
% $Id: spm_voice_get_word.m 7617 2019-06-13 12:01:17Z karl $


%% log prior over lexical content 
%==========================================================================
global VOX
if nargin < 2;
    nw = numel(VOX.LEX);
    P  = ones(nw,1)/nw;
end
LP = log(P + eps);

% within Ockham's window W
%--------------------------------------------------------------------------
W  = find(LP(:,1) > (max(LP(:,1)) - 3));


%% get onset of next word 
%==========================================================================
[Y,I,FS] = spm_voice_get_next(wfile);

% break if EOF
%--------------------------------------------------------------------------
if isempty(I), O  = {}; F = 0; return, end


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
nj   = numel(j);

% retrieve epochs and decompose at fundamental frequency
%--------------------------------------------------------------------------
clear xy J
for i = 1:nj
    xy(i,1)       = spm_voice_ff(y(j{i}),FS);
    xy(i,1).P.lat = log((I + j{i}(1) - VOX.IT)/FS);
    J(i,:)        = I + [j{i}(1),j{i}(end)];
end

% select the most likely action (interval)
%==========================================================================
[L,M,N] = spm_voice_likelihood(xy,W);
LI      = log(VOX.FI*P(:,1) + eps)';
L       = L + repmat(LP(:,1), 1,size(L,2));
L       = L + repmat(LI(1:nj),size(L,1),1);
Q       = spm_softmax(L);
F       = sum(Q.*(L - log(Q + eps)));

% if priors over the next word, accumulate free energy for each interval
%--------------------------------------------------------------------------
G       = spm_zeros(F);
if size(P,2) > 1 && nj > 1
    for i = 1:nj
        VOX.IT = fix(J(i,2));
        P      = spm_softmax(LP(:,2:end));
        [O,Fi] = spm_voice_get_word(wfile,P);
        G(i)   = Fi;
    end
    if all(G)
        F = F + G;
    elseif size(P,2) < VOX.depth
        F = 0;
    end
end
[F,m] = max(F);

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


