function [rej] = read_eep_rej(filename);

% READ_EEP_REJ reads rejection marks from an EEProbe *.rej file
%
% This function returns a Nx2 matrix with the begin and end latency
% of N rejection marks. The latency is in miliseconds.
%
% rej = read_eep_rej(filename)
%
% An EEP rejection file is formatted like
%   0.0000-0.3640
%   2.4373-3.5471
%   ... 
% where rejection begin and end are given in seconds. This function 
% converts the latency in miliseconds.
%
% See also READ_EEP_TRG, READ_EEP_CNT, READ_EEP_AVR

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_eep_rej.m,v $
% Revision 1.6  2004/05/27 09:30:27  roberto
% fixed bug in naming of variable
%
% Revision 1.5  2004/03/29 15:18:52  roberto
% minor changes to the naming of input arguments
%
% Revision 1.4  2003/03/24 12:34:33  roberto
% minor changes
%
% Revision 1.3  2003/03/13 14:25:45  roberto
% converted from DOS to UNIX
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

rej = [];

fid = fopen(filename, 'rb');
if fid<0
   return 
end

while ~feof(fid)
  tmp = fscanf(fid, '%f-%f', 2);
  if ~isempty(tmp)
    rej = [rej; tmp'];
  end
end

% convert to ms
rej = 1000*rej;

fclose(fid);  

