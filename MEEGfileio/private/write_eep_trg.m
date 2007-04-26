function write_eep_trg(fn, trg, freq);

% WRITE_EEP_TRG writes triggers to an EEProbe *.trg file
%
% [trg] = write_eep_trg(fn, trg, freq)
%
% with
%   fn     filename including the .trg extension
%   trg    structure array with the triggers in ms
%   freq   sampling frequency of the continuous data
%
% where trg contains the fields
%   trg.time	trigger latency in ms
%   trg.offset	byte offset in the file (not neccesary)
%   trg.code	trigger code as a string
%   trg.type    trigger code as a number
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
% See also READ_EEP_TRG, READ_EEP_CNT, READ_EEP_REJ, READ_EEP_AVR

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: write_eep_trg.m,v $
% Revision 1.2  2003/07/23 15:02:27  roberto
% added check on valid input for read_ctf_meg4, other changes unknown
%
% Revision 1.1  2003/03/18 12:15:49  roberto
% new implementation
%

fid = fopen(fn, 'wb');
if fid<0
   return
end

tmp = fprintf(fid, '%12.10f %d\n', 1/freq, freq); % I do not know what the 2nd number means

if ~isfield(trg, 'offset')
  % add a fake trigger offset
  for i=1:length(trg)
    trg(i).offset = 0;
  end
end

if ~isfield(trg, 'code') & isfield(trg, 'type')
  % convert numeric trigger codes into strings
  for i=1:length(trg)
    if isnan(trg(i).code)
      trg(i).code = '__';
    else
      trg(i).code = sprintf('%d', trg(i).type);
    end
  end
end

for i=1:length(trg)
  fprintf(fid, '%f\t%d\t%s\n', trg(i).time/1000, trg(i).offset, trg(i).code);
end

fclose(fid);  


