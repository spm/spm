function [i] = spm_voice_i(str)
% Gets indices or word strings from Lexington 
% FORMAT spm_voice_i(str)
% FORMAT spm_voice_i(w)
%
% str  - string or cell array
% i    - index

% requires the following in the global variable VOX:
% LEX    - lexical structure array
%
%  This routine returns the index or indices of a word if supplied with a
%  string or cell array. Alternatively, it returns the string
%  corresponding to an index all vector of indices.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_i.m 7562 2019-04-01 09:49:28Z karl $

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

% return cell array of indexed words
%--------------------------------------------------------------------------
if iscell(str)
    for w = 1:numel(str)
        i(w) = spm_voice_i(str{w});
    end
    return
end

% get indices
%--------------------------------------------------------------------------
i    = find(strcmp(str,word));






