function [trg] = read_eep_trg(filename);

% READ_EEP_TRG reads triggers from an EEProbe *.trg file
%
% This function returns an 1xN array with the triggers
%
% [trg] = read_eep_trg(filename)
%
% An EEP trigger file is formatted like
%   0.00195312 256
%   0.000     10350  __
%   17.033   2242926   1
%   20.535   2701934   5
%   21.096   2775406  13
%   21.098   2775662   8
%   ...
% where the first column is the trigger latency in seconds, the second
% column is the byte offset in the file and the third column is the triggercode. 
% The triggercode can be numeric or a string. The first line of the file contains the
% sample duration and ???.
%
% See also READ_EEP_CNT, READ_EEP_REJ, READ_EEP_AVR

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_eep_trg.m,v $
% Revision 1.4  2004/03/29 15:18:52  roberto
% minor changes to the naming of input arguments
%
% Revision 1.3  2003/03/13 14:25:45  roberto
% converted from DOS to UNIX
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

trg = [];

fid = fopen(filename, 'rb');
if fid<0
   return
end

header = fgetl(fid);
while ~feof(fid)
  tmp = fscanf(fid, '%f %d %s', 3);
  if ~isempty(tmp)
    new.time   = tmp(1)*1000;			% in ms
    new.offset = tmp(2)+1;			% offset 1
    new.code   = char(tmp(3:end));		% string
    new.type   = str2double(new.code);		% numeric
    trg = [trg; new];
  end
end

fclose(fid);  

