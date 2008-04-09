function data = read_asa_msr(fn);

% READ_ASA_MSR reads EEG or MEG data from an ASA data file
% converting the units to uV or fT

% Copyright (C) 2002, Robert Oostenveld
% 
% $Log: read_asa_msr.m,v $
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

Npnt = read_asa(fn, 'NumberPositions=', '%d');
Ntime = read_asa(fn, 'NumberTimesteps=', '%d');
UnitT = read_asa(fn, 'UnitTime', '%s');
UnitM = read_asa(fn, 'UnitMeas', '%s');
Timesteps = read_asa(fn, 'Timesteps', '%s');
lab  = read_asa(fn, 'Labels', '%s', Npnt);


val = read_asa(fn, 'Values', '%f');
if any(size(val)~=[Npnt,Ntime])
  msm_file = read_asa(fn, 'Values', '%s');
  [path, name, ext] = fileparts(fn);
  fid = fopen(fullfile(path, msm_file), 'rb');
  val = fread(fid, [Ntime, Npnt], 'float32')';
  fclose(fid);
end

tmp = sscanf(Timesteps, '%f(%f)%f');
time = linspace(tmp(1), tmp(3), Ntime);

if strcmp(lower(UnitT),'ms')
  time = 1*time;
elseif strcmp(lower(UnitT),'s')
  time = 1000*time;
elseif ~isempty(UnitT)
  error(sprintf('Unknown unit of time (%s)', UnitT));
end

if strcmp(lower(UnitM),'uv')
  val = 1*val;
elseif strcmp(lower(UnitM),'µv')
  val = 1*val;
elseif strcmp(lower(UnitM),'mv')
  val = 1000*val;
elseif strcmp(lower(UnitM),'v')
  val = 1000000*val;
elseif strcmp(lower(UnitM),'ft')
  val = 1*val;
elseif strcmp(lower(UnitM),'pt')
  val = 1000*val;
elseif ~isempty(UnitM)
  error(sprintf('Unknown unit of measurement (%s)', UnitM));
end

if length(size(lab))==2
  lab = tokenize(lab{1});
end

data.time  = time;
data.data  = val;
data.label = lab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [element] = tokenize(line);
count = 0;
while ~isempty(line)
  [tok, line] = strtok(line);
  if ~isempty(tok)
    count = count+1;
    element{count} = tok;
  else
    break;
  end
end

