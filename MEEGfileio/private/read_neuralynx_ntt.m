function [ntt] = read_neuralynx_ntt(filename, begrecord, endrecord)

% READ_NEURALYNX_NTT reads a single tetrode file
%
% Use as
%   [ntt] = read_neuralynx_ntt(filename)
%   [ntt] = read_neuralynx_ntt(filename, begrecord, endrecord)

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: read_neuralynx_ntt.m,v $
% Revision 1.4  2006/12/13 15:46:31  roboos
% read and keep timestamps as uint64
%
% Revision 1.3  2006/12/12 11:31:31  roboos
% cleaned up the code, made code more consistent with other neuralynx functions, moved subfunctions to seperate files, use numeric arrays instead of cell-arrays for storing the data
%
% Revision 1.2  2006/03/29 15:00:46  roboos
% changed some whitespace
%
% Revision 1.1  2005/05/20 06:39:35  roboos
% new implementation of single file readers
%

if nargin<2
  begrecord = 1;
end
if nargin<3
  endrecord = inf;
end

% The file starts with a 16*1024 bytes header in ascii, followed by a
% number of records (c.f. trials).

% The format of a tetrode record is
% int64 TimeStamp
% int32 ScNumber
% int32 CellNumber
% int32 Param[0] ƒ ƒ int32 Param[7]
% int16 ChanW[0]
% int16 ChanX[0]
% int16 ChanY[0]
% int16 ChanZ[0] ƒ ƒ
% int16 ChanW[31]
% int16 ChanX[31]
% int16 ChanY[31]
% int16 ChanZ[31]

hdr = neuralynx_getheader(filename);

fid = fopen(filename, 'rb', 'ieee-le');
% determine the length of the file
fseek(fid, 0, 'eof');
headersize = 16384;
recordsize = 304;
NRecords   = floor((ftell(fid) - headersize)/recordsize);

if begrecord<1
  error('cannot read before the first record');
elseif endrecord>NRecords
  endrecord = NRecords;
end

if endrecord>=begrecord
  % rewind to the first record to be read
  status = fseek(fid, headersize + (begrecord-1)*recordsize, 'bof');
  if status~=0
    error('cannot jump to the requested record');
  end

  numrecord    = (endrecord-begrecord+1);
  TimeStamp    = zeros(1,numrecord,'uint64');
  ScNumber     = zeros(1,numrecord);
  CellNumber   = zeros(1,numrecord);
  Param        = zeros(8,numrecord);
  Samp         = zeros(4,32,numrecord);

  for k=1:numrecord
    TimeStamp(k)  = fread(fid, 1, 'uint64=>uint64');
    ScNumber(k)   = fread(fid, 1, 'int32');
    CellNumber(k) = fread(fid, 1, 'int32');
    Param(k,:)    = fread(fid, 8, 'int32');
    Samp(:,:,k)   = fread(fid, [4 32], 'int16'); % chan W, X, Y, Z
  end

  ntt.TimeStamp  = TimeStamp;
  ntt.ScNumber   = ScNumber;
  ntt.CellNumber = CellNumber;
  ntt.Param      = Param;
  % apply the scaling factor from ADBitVolts and convert to uV
  ntt.dat        = Samp * 1e6 * hdr.ADBitVolts;
end
fclose(fid);

% store the header info in the output structure
ncs.NRecords = NRecords;
ntt.hdr      = hdr;
