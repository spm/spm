function [dat, lab] = read_bham(filename)

% READ_BHAM reads the EEG data files as recorded by Praamstra in Birmingham
% the datafiles are in a particular ascii format
%
% [dat, lab] = read_bham(filename)

% Copyright (C) 2000, Robert Oostenveld
% 
% $Log: read_bham.m,v $
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

fid = fopen(filename, 'rt');

lablen = 6;
line   = fgetl(fid);
numelc = 0;
while ~isempty(line)
  numelc = numelc + 1;
  [t, r] = strtok(line);
  lab(numelc,:) = [blanks(lablen-length(t)), t];
  line = r;
end

buf = fscanf(fid, '%f');
dat = zeros(numelc, length(buf)/numelc);
dat(:) = buf;

