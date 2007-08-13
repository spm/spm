function [tok] = tokenize(str, sep, rep)

% TOKENIZE cuts a string into pieces, returning a cell array
%
% Use as
%   t = tokenize(str, sep)
%   t = tokenize(str, sep, rep)
% where str is a string and sep is the separator at which you want
% to cut it into pieces.
%
% Using the optional boolean flag rep you can specify whether repeated
% seperator characters should be squeezed together (e.g. multiple
% spaces between two words). The default is rep=1, i.e. repeated
% seperators are treated as one.

% Copyright (C) 2003-2006, Robert Oostenveld
%
% $Log: tokenize.m,v $
% Revision 1.5  2006/05/10 07:15:22  roboos
% allow for squeezing multiple separators into one
%
% Revision 1.4  2006/01/30 14:56:41  roboos
% fixed another stupid bug in previous cvs commit
%
% Revision 1.3  2006/01/30 14:55:44  roboos
% fixed stupid bug in previous cvs commit
%
% Revision 1.2  2006/01/30 13:38:18  roboos
% replaced dependency on strtok by more simple code
% changed from dos to unix
%
% Revision 1.1  2005/05/23 13:47:51  roboos
% old implementation, new addition to CVS for fieldtrip release
%

tok = {};
f = find(str==sep);
f = [0, f, length(str)+1];
for i=1:(length(f)-1)
  tok{i} = str((f(i)+1):(f(i+1)-1));
end

if nargin<3 || rep
  % remove empty cells, which occur if the separator is repeated (e.g. multiple spaces)
  tok(find(cellfun('isempty', tok)))=[];
end

