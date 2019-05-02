function [L,F0,F1] = spm_voice_identity(wfile,P)
% Evaluates the likelihood of the next word in a file or object
% FORMAT [L,F0,F1] = spm_voice_identity(wfile,P)
%
% wfile  - .wav file, audiorecorder object or (double) time series
% P      - lexical prior [optional]
%
% L      - joint likelihood over F0 and F1 ( identity attributes)
% F0     - fundamental frequency
% F1     - format frequency
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
% this routine estimates the joint probability over fundamental and formant
% frequencies based upon a spoken word source, using the average spectral
% content of a spoken word.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_identity.m 7583 2019-05-02 12:10:08Z karl $

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
if isempty(I), L  = {}; return, end


%% get 1 second segment and remove previous word
%==========================================================================

% get intervals (j) for this peak
%----------------------------------------------------------------------
j    = fix((0:FS + FS/4) + I);
y    = Y(j(j < n & j > 1));
j    = logical(j < VOX.IT);
y(j) = 0;
J    = spm_voice_onsets(y,FS);

% get F0
%--------------------------------------------------------------------------
F0   = spm_voice_fundamental(Y,FS);


% likelihood search over fundamental and formant frequencies
%--------------------------------------------------------------------------
F1    = linspace(24,48,16);
F     = zeros(numel(F1),1);
for i = 1:numel(F1)
    
    VOX.F0 = F0;
    VOX.F1 = F1(i);

    % sentence
    %------------------------------------------------------------------
    clear xy
    for k = 1:numel(J)
        xy(k,1) = spm_voice_ff(y(J{k}),FS);
    end
    
    % free energy
    %------------------------------------------------------------------
    L    = spm_voice_likelihood(xy,W);
    L    = bsxfun(@plus,L,LP);
    Q    = spm_softmax(L);
    Fi   = sum(Q.*(L - log(Q + exp(-8))));
    F(i) = spm_softmax(Fi(:))'*Fi(:);

    
end

% plot(F1,F)
% posteriors
%--------------------------------------------------------------------------
L      = F;
L(:)   = spm_softmax(F(:));
F1     = F1*spm_vec(sum(L,2));

return




