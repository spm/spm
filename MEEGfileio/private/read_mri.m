function [mri] = read_mri(fn, s, offset);

% READ_MRI reads all MRI slices of size [Nx,Ny] from a 'uint8' binary file
%
% [mri] = read_mri(filename, s, offset)
%    size 	: [Nx, Ny, Nz]
%    offset	: the number of slices to skip
%    global fb	: give feedback
%
% See also READ_MRI2, READ_SHORT

% Copyright (C) 1999, Robert Oostenveld
%
% $Log: read_mri.m,v $
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

global fb;

fid = fopen(fn, 'rb');
if fid~=-1

  if nargin<2 | isempty(s)
    % assume MRI size of 256x256
    s = [256, 256, inf];
  elseif size(s,2)==1
    % assume square matrix size
      s = [s(1), s(1), inf];
  elseif size(s,2)==2
    % read all available slices
      s = [s(1), s(2), inf];
  end
  
  if nargin<3 | isempty(offset)
    % read the MRI starting at the first slice
    offset = 0;
  end

  if not(isempty(fb)) & fb==1
    disp([s, offset]);
  end

  i = 0;

  % seek to the first requested slice
  fseek(fid, s(1)*s(2)*offset, 'bof');

  while feof(fid)==0 & i<s(3)
    % read a single slice of the binary data
    if not(isempty(fb)) & fb==1
      disp(i);
    end
    [buf, count] = fread(fid, s([1, 2]), 'uint8');
    if size(buf) == s([1, 2])
      % add this slice to the MRI stack
      i = i + 1;
      mri(:,:,i) = uint8(buf);
    end
  end

  fclose(fid);

  if not(isempty(fb)) & fb==1
    imagesc(mri(:,:,i))
    %imagesc(permute(mri(:,:,i), [2,1]))
    %set(gca, 'YDir', 'normal');
  end

else
  error('unable to open file');
end

