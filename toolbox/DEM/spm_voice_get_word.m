function [O,I,J,F] = spm_voice_get_word(wfile,P)
% Evaluate the likelihood of the next word in a file or object
% FORMAT [O,I,J,F] = spm_voice_get_word(wfile,P)
%
% wfile  - .wav file, audiorecorder object or (double) time series
% P      - lexical prior probability [optional]
%
% O{1}   - lexical likelihood (or posterior if priors are specified)
% O{2}   - prosody likelihood
% O{3}   - speaker likelihood
%
% I      - interval index (1/2 sec. before spectral peak)
% J      - interval onset and offset
% F      - maxmium F (log evidence)
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
% IT     - index or pointer to offset of current word
%
% This routine evaluates the likelihood of a word, prosody and identity by
% inverting successive epochs of data from an audiofile or device starting
% at VOX.IT.  Based on the word with the least variational free energy, it
% updates the index, ready for the next word. Priors over the words can be
% specified to implement an Occam's window (of 3 nats); thereby
% restricting the number of lexical entries evaluated - and augmenting the
% likelihoods to give the posterior probability over words.
% If called with more than one prior over lexical content, this routine
% will perform a tree search and return the likelihoods (and intervals)
% with the path of greatest log evidence (i.e., free energy).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_get_word.m 7750 2019-12-05 17:54:29Z spm $


%% log prior over lexical content
%==========================================================================
global VOX
if nargin < 2 || size(P,2) < 1
    nw = numel(VOX.LEX);
    P  = ones(nw,1)/nw;
end
LP     = log(P + eps);                        % log prior for lexicon
try LL = VOX.LL; catch, LL = -512; end        % log prior for empty segment
try LW = VOX.LW; catch, LW = 0;    end        % flag for last word

%% get onset of next word
%==========================================================================
[Y,I,FS] = spm_voice_get_next(wfile);


% specify output arguments and break if EOF
%--------------------------------------------------------------------------
if isempty(I)
    
    O = {};
    J = [];
    F = LL;
    return
    
elseif I < 1
    
    % deal with negative I (inital index)
    %----------------------------------------------------------------------
    i    = 1:FS/32;
    Y(i) = 0;
    I    = 1;
    
end



%% get 1.5 second segment and remove previous word
%==========================================================================

% get intervals (j) for this peak
%--------------------------------------------------------------------------
n     = numel(Y);
j     = fix((0:FS + FS/2) + I);
y     = Y(j(j < n));
j     = logical(j < VOX.IT);
y(j)  = 0;

% get intervals for this segment
%--------------------------------------------------------------------------
j     = spm_voice_onsets(y,FS,1/16);
if LW 
    j = j(end);                                     % last word
end
nj    = numel(j);

% (deep) search
%==========================================================================

% loop over intervals
%--------------------------------------------------------------------------
np    = size(P,2);                                  % number of priors
F     = zeros(nj,1);                                % free energy
J     = zeros(nj,2);                                % onsets and offsets
for i = 1:nj
    
    % retrieve epoch and decompose
    %----------------------------------------------------------------------
    xy       = spm_voice_ff(y(j{i}),FS);            % transform
    lat      = I + j{i}(1) - VOX.IT;                % latency
    xy.P.lat = log(max(lat/FS,1/32));               % log latency
    J(i,:)   = I + [j{i}(1),j{i}(end)];             % onsets and offsets
    
    % Ockham's window (W)
    %------------------------------------------------------------------
    W        = find(LP(:,1) > (max(LP(:,1)) - 3));
    
    % select the most likely action (interval)
    %==================================================================
    [L,M,N]  = spm_voice_likelihood(xy,W);          % likelihoods 
    L        = L + LP(:,1);                         % prior (lexical)
    L        = L + log(VOX.FI(i,:)*P(:,1) + eps);   % prior (peaks)
    Q        = spm_softmax(L);                      % posterior
    F(i)     = sum(Q.*(L - log(Q + eps)));          % free energy
    
    % posteriors
    %----------------------------------------------------------------------
    O{i}     = {L,M,N};                             % log-likelihoods
    
    %  deep search if there are priors over future words
    %----------------------------------------------------------------------
    if np > 1
        
        % move the pointer to end of current interval
        %------------------------------------------------------------------
        VOX.IT = fix(J(i,2));
        Pi     = spm_softmax(LP(:,2:end));

        [Oi,Ii,Ji,Fi] = spm_voice_get_word(wfile,Pi);
        
        % accumulate free energy
        %------------------------------------------------------------------
        F(i) = Fi + F(i);                           
        
    end
    
end

% retain interval with the greatest free energy
%==========================================================================
[F,i] = max(F);

% likelihood and pointer based upon selected interval
%----------------------------------------------------------------------
O      = O{i};
I      = [I; J(:,2)];
J      = fix(J(i,:));
VOX.IT = fix(J(2));
