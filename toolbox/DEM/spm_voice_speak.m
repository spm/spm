function spm_voice_speak(w,p,LEX,PRO,DT)
% inverse decomposition at fundamental frequency
% FORMAT spm_voice_speak(w,p,LEX,PRO,DT)
%
% w      - lexcial index (1 x number of words)
% p      - prosody index (k x number of words)
% LEX    - lexical structure array
% PRO    - prosody structure array
% DT     - optional latencies (seconds)
%
% This routine recomposes and plays a timeseries, specified as a sequence of
% words that can be articulated with a particular prosody.  This routine
% plays the same role as spm_voice_iff but uses the  dictionaries supplied
% by the lexical and prosody structures to enable a categorical
% specification of a spoken phrase. In other words, it allows one to map
% from discrete state space of lexical content and prosody to continuous
% time outcomes.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_speak.m 7528 2019-02-06 19:19:49Z karl $


% check for empty indices (that will invoke average lexical or prosody)
%--------------------------------------------------------------------------
if isempty(w), w = zeros(0,1); end
if isempty(p), p = zeros(0,1); end
n  = max(size(w,2),size(p,2));
if size(w,2) < n
    w = repmat(w(:,1),1,n);
elseif size(p,2) < n
    p = repmat(p(:,1),1,n);
end


% assemble word structure arrays
%==========================================================================
for s = 1:n
    
    % lexical parameters
    %----------------------------------------------------------------------
    Q     = 0;
    for i = 1:size(w,1)
        Q = Q + LEX(1).U*LEX(w(1,s),1).pE;
    end
    xY(s).Q = spm_unvec(spm_vec(LEX(1).Q) + Q,LEX(1).Q);
    
    % prosody parameters
    %----------------------------------------------------------------------
    P     = 0;
    for i = 1:size(p,1)
        P = P + PRO(1).U(:,i)*PRO(i).pE(p(i,s));
    end
    xY(s).P = spm_unvec(spm_vec(PRO(1).P) + P,PRO(1).P);
    
end

% sent to spm_voice_iff
%--------------------------------------------------------------------------
if nargin > 4
    spm_voice_iff(xY,16000,DT);
else
    spm_voice_iff(xY,16000);
end

