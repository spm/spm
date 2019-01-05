function [L,i,str] = spm_voice_likelihood(xY,LEX)
% returns the lexical likelihood
% FORMAT [L,i] = spm_voice_likelihood(xY,LEX);
%
% xY   - lexical data structure
% LEX  - lexical dictionary array
%
% L    - likelihood over k words
% i    - index of most likely word
% str  - most likely word string
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_likelihood.m 7512 2019-01-05 21:15:16Z karl $

% likelihood of XY
%==========================================================================
for w = 1:size(LEX,1)
    for k = 1:size(LEX,2)
        P.y    = xY.y;
        P.ci   = xY.ci;
        P.ni   = xY.ni;
        P.ti   = xY.ti;
        e      = LEX(w,1).U'*(spm_vec(LEX(w,k).pE) - spm_vec(P));
        L(w,k) = -e'*LEX(w,k).pP*e/2;
    end
end

% likelihood over lexical entries
%--------------------------------------------------------------------------
L(:)  = spm_softmax(L(:));
L     = sum(L,2);
[d,i] = max(L);
str   = LEX(i,1).str;






