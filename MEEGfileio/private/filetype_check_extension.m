function [val] = filetype_check_extension(filename, ext)

% FILETYPE_CHECK_EXTENSION helper function to determine the file type
% by performing a case insensitive string comparison of the extension.

% Copyright (C) 2003-2006 Robert Oostenveld
%
% $Log: filetype_check_extension.m,v $
% Revision 1.2  2007/03/21 17:22:08  roboos
% use recursion instead of copying the same code twice
%
% Revision 1.1  2006/08/28 10:10:11  roboos
% moved subfunctions for filetype detection into seperate files
%

if iscell(filename)
  % compare the extension of multiple files
  val = zeros(size(filename));
  for i=1:length(filename)
    val(i) = filetype_check_extension(filename{i}, ext);
  end
else
  % compare the extension of a single file
  if length(filename)<length(ext)
    val = 0;
  else
    val = strcmpi(filename((end-length(ext)+1):end), ext);
  end
end
return

