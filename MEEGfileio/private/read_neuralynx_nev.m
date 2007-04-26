function [nev] = read_neuralynx_nev(filename);

% READ_NEURALYNX_NEV reads the event information from the *.nev file in a
% Neuralynx dataset directory
%
% Use as
%   nev = read_neuralynx_hdr(datadir)
%   nev = read_neuralynx_hdr(eventfile)
%
% The output structure contains all events and timestamps.

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: read_neuralynx_nev.m,v $
% Revision 1.1  2006/12/13 15:43:49  roboos
% renamed read_neuralynx_event into xxx_nev, consistent with the file extension
%
% Revision 1.3  2005/09/05 13:08:52  roboos
% implemented a faster way of reading all events, needed in the case of many triggers (e.g. each refresh)
%
% Revision 1.2  2005/06/24 06:57:32  roboos
% added PktStart to the record header reading, this shifts all fields by two bytes (thanks to Thilo)
%
% Revision 1.1  2005/05/19 07:09:58  roboos
% new implementation
%

if filetype(filename, 'neuralynx_ds')
    % replace the directory name by the filename
    filename = fullfile(filename, 'Events.Nev');
end

% The file starts with a 16*1024 bytes header in ascii, followed by a
% number of records (c.f. trials).
%
% The format of an event record is
%   int16 PktId
%   int16 PktDataSize
%   int64 TimeStamp
%   int16 EventId
%   int16 TTLValue
%   int16 CRC
%   int32 Dummy
%   int32 Extra[0]
%   ...
%   int32 Extra[7]
%   char EventString[0]
%   ...
%   char EventString[127]
% PktId is usually 0x1002.
% PktDataSize is random data.
% Dummy is random data.
% CRC may contain random data.
% Extra is user-defined data.
% TTLValue is the value sent to the computer on a parallel input port.

% read all event records
fid = fopen(filename, 'rb', 'ieee-le');
fseek(fid, 16*1024, 'bof');
nev = [];

% this is the slow way of reading it
% while ~feof(fid)
%   nev(end+1).PktStart     = fread(fid, 1, 'int16');
%   nev(end  ).PktId        = fread(fid, 1, 'int16');
%   nev(end  ).PktDataSize  = fread(fid, 1, 'int16');
%   nev(end  ).TimeStamp    = fread(fid, 1, 'int64');
%   nev(end  ).EventId      = fread(fid, 1, 'int16');
%   nev(end  ).TTLValue     = fread(fid, 1, 'int16');
%   nev(end  ).CRC          = fread(fid, 1, 'int16');
%   nev(end  ).Dummy        = fread(fid, 1, 'int32');
%   nev(end  ).Extra        = fread(fid, 8, 'int32');
%   nev(end  ).EventString  = fread(fid, 128, 'char');
% end

% this is a faster way of reading it
% and it is still using the automatic type conversion from Matlab
fp = 16*1024;
fseek(fid, fp+ 0, 'bof'); PktStart       = fread(fid, inf, 'uint16', 184-2);
num = length(PktStart);
fseek(fid, fp+ 2, 'bof'); PktId          = fread(fid, num, 'uint16', 184-2);
fseek(fid, fp+ 4, 'bof'); PktDataSize    = fread(fid, num, 'uint16', 184-2);
fseek(fid, fp+ 6, 'bof'); TimeStamp      = fread(fid, num, 'uint64', 184-8);
fseek(fid, fp+14, 'bof'); EventId        = fread(fid, num, 'uint16', 184-2);
fseek(fid, fp+16, 'bof'); TTLValue       = fread(fid, num, 'uint16', 184-2);
fseek(fid, fp+18, 'bof'); CRC            = fread(fid, num, 'uint16', 184-2);
fseek(fid, fp+22, 'bof'); Dummy          = fread(fid, num, 'int32' , 184-4);
% read each of the individual extra int32 values and concatenate them
fseek(fid, fp+24+0*4, 'bof'); Extra1     = fread(fid, num, 'int32', 184-4);
fseek(fid, fp+24+1*4, 'bof'); Extra2     = fread(fid, num, 'int32', 184-4);
fseek(fid, fp+24+2*4, 'bof'); Extra3     = fread(fid, num, 'int32', 184-4);
fseek(fid, fp+24+3*4, 'bof'); Extra4     = fread(fid, num, 'int32', 184-4);
fseek(fid, fp+24+4*4, 'bof'); Extra5     = fread(fid, num, 'int32', 184-4);
fseek(fid, fp+24+5*4, 'bof'); Extra6     = fread(fid, num, 'int32', 184-4);
fseek(fid, fp+24+6*4, 'bof'); Extra7     = fread(fid, num, 'int32', 184-4);
fseek(fid, fp+24+7*4, 'bof'); Extra8     = fread(fid, num, 'int32', 184-4);
Extra = [Extra1 Extra2 Extra3 Extra4 Extra5 Extra6 Extra7 Extra8];
% read the complete data excluding header as char and cut out the piece with the EventString content
fseek(fid, fp, 'bof'); EventString        = fread(fid, [184 num], 'char');
EventString = char(EventString(57:184,:)');
fclose(fid);

% restructure the data into a struct array
% first by making cell-arrays out of it ...
PktStart      = mat2cell(PktStart   , ones(1,num), 1);
PktId         = mat2cell(PktId      , ones(1,num), 1);
PktDataSize   = mat2cell(PktDataSize, ones(1,num), 1);
TimeStamp     = mat2cell(TimeStamp  , ones(1,num), 1);
EventId       = mat2cell(EventId    , ones(1,num), 1);
TTLValue      = mat2cell(TTLValue   , ones(1,num), 1);
CRC           = mat2cell(CRC        , ones(1,num), 1);
Dummy         = mat2cell(Dummy      , ones(1,num), 1);
Extra         = mat2cell(Extra      , ones(1,num), 8);
EventString   = mat2cell(EventString , ones(1,num), 128);
% ... and then convert the cell-arrays into a single structure
nev = struct(...
'PktStart'    , PktStart     , ...
'PktId'       , PktId        , ...
'PktDataSize' , PktDataSize  , ...
'TimeStamp'   , TimeStamp    , ...
'EventId'     , EventId      , ...
'TTLValue'    , TTLValue     , ...
'CRC'         , CRC          , ...
'Dummy'       , Dummy        , ...
'Extra'       , Extra        , ...
'EventString' , EventString);

% remove null values and convert to strings
for i=1:length(nev)
  nev(i).EventString = nev(i).EventString(find(nev(i).EventString));
  nev(i).EventString = char(nev(i).EventString(:)');
end

