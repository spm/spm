function [xY] = spm_voice_speak(w,p,q)
% inverse decomposition at fundamental frequency
% FORMAT spm_voice_speak(w,p,q)
%
% w      - lexcial index (1 x number of words)
% p      - prosody index (k x number of words)
% q      - prosody index (2 x number of words)

% requires the following in the global variable voice_options:
% LEX    - lexical structure array
% PRO    - prodidy structure array
% WHO    - speaker structure array
%
% This routine recomposes and plays a timeseries, specified as a sequence
% of words that can be articulated with a particular prosody.  This routine
% plays the same role as spm_voice_iff but uses the dictionaries supplied
% by the lexical and prosody structures to enable a categorical
% specification of a spoken phrase. In other words, it allows one to map
% from discrete state space of lexical content and prosody to continuous
% time outcomes.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_speak.m 7535 2019-03-03 20:27:09Z karl $

% check for empty indices (that will invoke average lexical or prosody)
%--------------------------------------------------------------------------
global voice_options
LEX = voice_options.LEX;
PRO = voice_options.PRO;
WHO = voice_options.WHO;

if isempty(w), w = zeros(0,1); end
if isempty(p), p = zeros(0,1); end
if isempty(q), q = zeros(0,1); end

n  = max([size(w,2),size(p,2),size(q,2)]);

if size(w,2) < n, w = repmat(w(:,1),1,n); end
if size(p,2) < n, p = repmat(p(:,1),1,n); end
if size(q,2) < n, q = repmat(q(:,1),1,n); end

% assemble word structure arrays
%==========================================================================
for s = 1:n
    
    % lexical parameters
    %----------------------------------------------------------------------
    Q     = 0;
    for i = 1:size(w,1)
        Q = Q + LEX(w(1,s),1).qE;
    end
    xY(s).Q = spm_unvec(spm_vec(LEX(1).Q) + Q,LEX(1).Q);
    
    % prosody parameters
    %----------------------------------------------------------------------
    P     = spm_vec(spm_zeros(PRO(1).P));
    for i = 1:size(p,1)
        j    = PRO(1).i(i);
        P(j) = P(j) + PRO(i).pE(p(i,s));
    end
    
    % and identity
    %----------------------------------------------------------------------
    for i = 1:size(q,1)
        j    = WHO(1).i(i);
        P(j) = P(j) + WHO(i).pE(q(i,s));
    end
    xY(s).P = spm_unvec(spm_vec(PRO(1).P) + P,PRO(1).P);
    
end

% assemble sequence
%--------------------------------------------------------------------------
if n > 1
    
    % turn off graphics, get FS and convert into audio signal 
    %----------------------------------------------------------------------
    g = voice_options.graphics;
    voice_options.graphics = 0;
    try, FS = voice_options.FS; catch, FS  = 22050; end
    
    for s = 1:n
        y{s} = spm_voice_iff(xY(s));
    end
    Y   = zeros(spm_length(y),1);
    i0  = 0;
    for s = 1:n
        ni    = numel(y{s});
        ii    = i0 + (1:ni)';
        Y(ii) = Y(ii) + y{s};
        i0    = ii(end) - round(ni/4);
    end
    Y   = Y(1:ii(end));
    
    
    % send to speaker
    %----------------------------------------------------------------------
    sound(full(Y),FS);
    voice_options.graphics = g;
    
else
    
    % and send to spm_voice_iff
    %----------------------------------------------------------------------
    spm_voice_iff(xY,1/2);
    
end

return


