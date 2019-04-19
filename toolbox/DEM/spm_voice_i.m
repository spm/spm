function [i,P] = spm_voice_i(str)
% Gets indices or word strings from Lexington 
% FORMAT spm_voice_i(str)
% FORMAT spm_voice_i(w)
%
% str  - string or cell array
% i    - index
% P    - corresponding array of prior probabilities

% requires the following in the global variable VOX:
% LEX    - lexical structure array
%
%  This routine returns the index or indices of a word if supplied with a
%  string or cell array. Alternatively, it returns the string
%  corresponding to an index all vector of indices.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_i.m 7574 2019-04-19 20:38:15Z karl $

% get timeseries from audio recorder(or from a path
%--------------------------------------------------------------------------


% words in lexicon
%==========================================================================
global VOX
word  = {VOX.LEX(:,1).word};                  % words in lexicon

% return cell array of indexed words
%--------------------------------------------------------------------------
if isnumeric(str)
    i = word(str);
    return
end

% return indices
%--------------------------------------------------------------------------
if iscell(str)
    if nargout < 2
        for w = 1:numel(str)
            i(w) = spm_voice_i(str{w});
        end
        return
    else
        
        % return indices and prior probabilities
        %------------------------------------------------------------------
        nw    = numel(str);
        P     = zeros(numel(VOX.LEX),nw) + 1/128;
        for w = 1:nw
            i      = spm_voice_i(str{w});
            P(i,w) = 1;
        end
        
        % sum to one constraint
        %------------------------------------------------------------------
        P     = bsxfun(@rdivide,P,sum(P));
        return
    end
end

% get indices
%--------------------------------------------------------------------------
i    = find(strcmp(str,word));






