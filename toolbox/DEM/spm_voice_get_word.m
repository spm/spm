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
% $Id: spm_voice_get_word.m 7576 2019-04-23 09:22:44Z karl $

%% get peak identification parameters from VOX
%==========================================================================

% get source (recorder) and FS
%--------------------------------------------------------------------------
[FS,read] = spm_voice_FS(wfile);

% defaults
%--------------------------------------------------------------------------
global VOX
try, VOX.C;  catch, VOX.C  = 1/16;  end              % smoothing for peaks
try, VOX.U;  catch, VOX.U  = 1/128; end              % threshold for peaks
try, VOX.IT; catch, VOX.IT = 1;     end              % final index

% ensure 2 second of data has been accumulated
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    dt     = (get(wfile,'TotalSamples') - VOX.IT)/FS;
    pause(2 - dt);
end

%% log prior over lexical content 
%==========================================================================
if nargin < 2
    nw = numel(VOX.LEX);
    LP = ones(nw,1)/nw;
else
    LP = log(P + exp(-8));
end

% within Ockham's window W
%--------------------------------------------------------------------------
W      = find(LP > (max(LP) - 3));

%% find word and evaluate likelihoods
%==========================================================================

% find next word (waiting for a couple of seconds if necessary)
%--------------------------------------------------------------------------
for i = 1:4
    
    % find next spectral peak (I)
    %----------------------------------------------------------------------
    Y = read(wfile);
    n = numel(Y);
    j = fix((0:FS) + VOX.IT);
    G = spm_voice_check(Y(j(j < n)),FS,VOX.C);
    I = find((diff(G(1:end - 1)) > 0) & (diff(G(2:end)) < 0));
    I = I(G(I) > VOX.U & I > FS/8);
    
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
        I  = VOX.IT + I(1) - FS/2;
        break
    end
end

% break if EOF
%--------------------------------------------------------------------------
if isempty(I), O  = {}; return, end


%% get 1 second segment and remove previous word
%==========================================================================

% get intervals (j) for this peak
%----------------------------------------------------------------------
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
F       = sum(Q.*(L - log(Q + exp(-8))));
[d,m]   = max(F(1,1,:));

% posteriors
%--------------------------------------------------------------------------
O{1}    = spm_softmax(L(:,:,m));            % posteriors
O{2}    = spm_softmax(M(:,:,m));            % likelihood
O{3}    = spm_softmax(N(:,:,m));            % likelihood

% update pointers based upon selected word
%--------------------------------------------------------------------------
VOX.I0  = fix(J(m,1));
VOX.IT  = fix(J(m,2));

return


