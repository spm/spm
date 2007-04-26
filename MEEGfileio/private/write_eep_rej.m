function write_eep_rej(fn, rej);

% WRITE_EEP_REJ writes rejection marks to an EEProbe *.rej file
%
% write_eep_rej(filename, reject)
% 
% The rejection marks consist of an Nx2 matrix with the begin and end
% latency of N rejection marks. The latency is in miliseconds.
%
% An EEP rejection file is formatted like
%   0.0000-0.3640
%   2.4373-3.5471
%   ... 
% where rejection begin and end are given in seconds. This function 
% converts the latency from miliseconds to seconds.
%
% See also READ_EEP_TRG, READ_EEP_CNT, READ_EEP_AVR

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: write_eep_rej.m,v $
% Revision 1.1  2003/03/20 10:42:04  roberto
% new implementation
%

fid = fopen(fn, 'wb');
if fid<0
   return
end

% convert to seconds
rej = rej/1000;

for i=1:size(rej,1)
  fprintf(fid, '%f-%f\n', rej(i,1), rej(i,2));
end

fclose(fid);
