function [t] = neuralynx_timestamp(filename, num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for reading a single timestamp of a single channel Neuralynx file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename, 'rb', 'ieee-le');
headersize = 16384;
switch filetype(filename)
  case 'neuralynx_ncs'
    recordsize = 1044;  % in bytes
  case 'neuralynx_nse'
    recordsize = 112;   % in bytes
  case 'neuralynx_nts'
    recordsize = 8;     % in bytes
end
if ~isinf(num)
  % read the timestamp of the indicated record
  fseek(fid, headersize + (num-1)*recordsize, 'bof');
  t = fread(fid, 1, 'uint64=>uint64');
else
  % read the timestamp of the last record
  fseek(fid, -recordsize, 'eof');
  t = fread(fid, 1, 'uint64=>uint64');
end
fclose(fid);