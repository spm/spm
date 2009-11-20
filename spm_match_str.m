function [sel1, sel2] = spm_match_str(a, b)
% MATCH_STR looks for matching labels in two listst of strings
% and returns the indices into both the 1st and 2nd list of the matches.
% They will be ordered according to the first input argument.
%  
% [sel1, sel2] = match_str(strlist1, strlist2)
%
% The strings can be stored as a char matrix or as an vertical array of
% cells, the matching is done for each row.
%_______________________________________________________________________
% Copyright (C) 2000, Robert Oostenveld

% Robert Oostenveld
% $Id: spm_match_str.m 3589 2009-11-20 17:17:41Z guillaume $

% ensure that both are cell-arrays
if isempty(a)
  a = {};
elseif ~iscell(a)
  a = cellstr(a);
end
if isempty(b)
  b = {};
elseif ~iscell(b)
  b = cellstr(b);
end

% ensure that both are column vectors
a = a(:);
b = b(:);

% regardless of what optimizations are implemented, the code should remain
% functionally compatible to the original, which is 
% for i=1:length(a)
%   for j=1:length(b)
%     if strcmp(a(i),b(j))
%       sel1 = [sel1; i];
%       sel2 = [sel2; j];
%     end
%   end
% end

% replace all unique strings by a unique number and use the fact that
% numeric comparisons are much faster than string comparisons
[dum1, dum2, c] = unique([a; b]);
a = c(1:length(a));
b = c((length(a)+1):end);

sel1 = [];
sel2 = [];
for i=1:length(a)
  % s = find(strcmp(a(i), b));  % for string comparison
  s = find(a(i)==b);            % for numeric comparison
  sel1 = [sel1; repmat(i, size(s))];
  sel2 = [sel2; s];
end
