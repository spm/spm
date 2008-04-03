function [val] = filetype_check_header(filename, head, offset)

% FILETYPE_CHECK_HEADER helper function to determine the file type
% by reading the first number of bytes of a file and comparing them
% to a known string (c.f. magic number).

% Copyright (C) 2003-2006 Robert Oostenveld
%
% $Log: filetype_check_header.m,v $
% Revision 1.6  2008/01/15 07:57:35  roboos
% fixed problem with fclose, which could happen twice
%
% Revision 1.5  2008/01/14 11:13:16  roboos
% fixed bug: files would not be closed, resulting in problems after opening too many files
%
% Revision 1.4  2007/12/17 16:15:22  roboos
% added support for checking the file header against multiple options
%
% Revision 1.3  2007/12/12 14:39:47  roboos
% added optional offset
%
% Revision 1.2  2007/03/21 17:22:46  roboos
% added support for checking the header of all files contained within a directory
%
% Revision 1.1  2006/08/28 10:10:11  roboos
% moved subfunctions for filetype detection into seperate files
%

if nargin<3
  offset = 0;
end

if iscell(filename)
  % compare the header of multiple files
  val = zeros(size(filename));
  for i=1:length(filename)
    val(i) = filetype_check_header(filename{i}, head, offset);
  end
elseif isdir(filename)
  % a directory cannot have a header
  val = 0;
else
  % read the first few bytes from the file and compare them to the desired header
  fid = fopen(filename, 'rb');
  if fid<0
    warning(sprintf('could not open %s', filename));
    val = 0;
  else
    fseek(fid, offset, 'cof');
    if iscell(head)
      for i=1:length(head)
        len(i) = length(head{i});
      end
      [str, siz] = fread(fid, max(len), 'char=>char');
      fclose(fid);
      for i=1:length(head)
        val = strncmp(str, head{i}, len(i));
        if val
          break
        end
      end
    else
      [str, siz] = fread(fid, length(head), 'char=>char');
      fclose(fid);
      if siz~=length(head)
        error('could not read the header from the file');
      end
      val = all(str(:)==head(:));
    end
  end
end
return

