function [i,P] = spm_voice_i(str)
% Get indices, word strings or priors from lexicon
% FORMAT [str] = spm_voice_i(i)
% FORMAT [i  ] = spm_voice_i(str)
% FORMAT [i,P] = spm_voice_i(str)
%
% str  - string or cell array
% i    - index in lexicon (VOX.LEX)
% P    - corresponding array of prior probabilities

% requires the following in the global variable VOX:
% LEX  - lexical structure array
%
%  This routine returns the index or indices of a word if supplied with a
%  string or cell array. Alternatively, it returns the string corresponding
%  to an index or vector of indices. if called with the second argument, it
%  returns a prior probability matrix where the specified words (for
%  strings) have a prior odds ratio of 64.
%
%  NB: If a string is not in the lexicon, 0 is returned.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_i.m 7750 2019-12-05 17:54:29Z spm $


% get timeseries from audio recorder(or from a path
%--------------------------------------------------------------------------

% words in lexicon
%==========================================================================
global VOX
word  = {VOX.LEX.word};                  % words in lexicon

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
        P     = zeros(numel(VOX.LEX),nw) + 1/64;
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
if isempty(i), i = 0; end
