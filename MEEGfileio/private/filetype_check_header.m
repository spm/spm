function [val] = filetype_check_header(filename, head)

% FILETYPE_CHECK_HEADER helper function to determine the file type
% by reading the first number of bytes of a file and comparing them 
% to a known string (c.f. magic number).

% Copyright (C) 2003-2006 Robert Oostenveld
%
% $Log: filetype_check_header.m,v $
% Revision 1.2  2007/03/21 17:22:46  roboos
% added support for checking the header of all files contained within a directory
%
% Revision 1.1  2006/08/28 10:10:11  roboos
% moved subfunctions for filetype detection into seperate files
%

if iscell(filename)
  % compare the header of multiple files
  val = zeros(size(filename));
  for i=1:length(filename)
    val(i) = filetype_check_header(filename{i}, head);
  end
elseif isdir(filename)
  val = 0;
else
  fid = fopen(filename, 'rb');
  if fid<0
    warning(sprintf('could not open %s', filename));
    val = 0;
  else
    [str, siz] = fread(fid, length(head), 'char');
    fclose(fid);
    if siz~=length(head)
      error('could not read the header from the file');
    end
    val = all(char(str(:))==head(:));
  end
end
return

