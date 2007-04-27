function [data] = appenddata(cfg, varargin);

% APPENDDATA combines multiple datasets that have been preprocessed separately 
% into a single large dataset.
%
% Use as
%   data = appenddata(cfg, data1, data2, data3, ...)
% where the configuration can be empty.
%
% If the input datasets all have the same channels, the trials will be
% concatenated. This is usefull for example if you have different
% experimental conditions, which, besides analyzing them separately, for
% some reason you also want to analyze together.
%
% If the input datasets have different channels, but the same number of
% trials, the channels will be concatenated within each trial. This is
% usefull for example if the data that you want to analyze contains both
% MEG and EMG channels which require different preprocessing options.
%
% See also PREPROCESSING

% undocumented options:
%   none

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: appenddata.m,v $
% Revision 1.12  2007/01/09 09:50:39  roboos
% some changes in whitespace
%
% Revision 1.11  2007/01/04 17:06:33  roboos
% implemented support for appending spike data (i.e. timestamps) to a continuously sampled data set (i.e. one obtained from PREPROCESSING)
%
% Revision 1.10  2006/07/24 11:29:29  roboos
% use private/findcfg function for locating the trl and event in the nested (previous) cfgs
%
% Revision 1.9  2006/05/02 19:17:04  roboos
% search for th trl in cfg and cfg.previous etc., and also append it
%
% Revision 1.8  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.7  2006/01/30 12:15:01  roboos
% ensure consistent function declaration: "function [x] = funname()"
%
% Revision 1.6  2005/08/05 09:16:22  roboos
% removed the obsolete data.offset
%
% Revision 1.5  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.4  2005/03/14 11:17:29  jansch
% fixed small bug in handling of dimensions in the case of concatenating trials
%
% Revision 1.3  2005/02/16 17:07:13  roboos
% fixed error that occurred when input data label was not a column-cell-array
%
% Revision 1.2  2005/01/27 12:06:33  roboos
% fixed bug in checking of identical time axes
%
% Revision 1.1  2005/01/17 14:39:43  roboos
% new implementation
%

% set the defaults
cfg = [];

Ndata = nargin-1;
if Ndata<2
  error('you must give at least two datasets to append');
end

% determine the type of input data, supported are continuous (i.e. segmented and preprocessed continuous data) and spike data
iscontinuous = zeros(1,Ndata);
isspike      = zeros(1,Ndata);
for i=1:Ndata
  iscontinuous(i) = isfield(varargin{i}, 'trial') && isfield(varargin{i}, 'time') && iscell(varargin{i}.trial) && iscell(varargin{i}.time);
  isspike(i)      = isfield(varargin{i}, 'waveform') && isfield(varargin{i}, 'timestamp') && iscell(varargin{i}.waveform) && iscell(varargin{i}.timestamp);
end

% determine the dimensions of the data
Nchan  = zeros(1,Ndata);
Ntrial = zeros(1,Ndata);
label  = {};
for i=1:Ndata
  Nchan(i) = length(varargin{i}.label);
  if iscontinuous(i)
    Ntrial(i) = length(varargin{i}.trial);
    fprintf('input dataset %d, %d channels, %d trials\n', i, Nchan(i), Ntrial(i));
  elseif isspike(i)
    fprintf('input dataset %d, %d channels with spikes\n', i, Nchan(i));
  end
  label = [label(:); varargin{i}.label(:)];
end

% try to locate the trial definition (trl) in the nested configuration
for i=1:Ndata
  if isfield(varargin{i}, 'cfg')
    trl{i} = findcfg(varargin{i}.cfg, 'trl');
  else
    trl{i} = [];
  end
  if isempty(trl{i}) && iscontinuous(i)
    % a trial definition is expected in each continuous data set
    warning(sprintf('could not locate the trial definition ''trl'' in data structure %d', i));
  end
end

if ~all(iscontinuous) && ~all(isspike)
  % the input consists of both continuous data and spike data
  % convert the spike data into a continous representation
  firstcontinuous = find(iscontinuous);
  firstcontinuous = firstcontinuous(1);
  % copy the timeaxis from the first continuous dataset, and determine the correspondence between samples and timestamps
  FirstTimeStamp     = varargin{firstcontinuous}.hdr.FirstTimeStamp;      % this corresponds with sample 1 in the original continuous data
  LastTimetamp       = varargin{firstcontinuous}.hdr.LastTimeStamp;       % this corresponds with sample hdr.nSample*hdr.nTrial
  TimeStampPerSample = varargin{firstcontinuous}.hdr.TimeStampPerSample;
  spikesel = find(isspike);
  for i=1:length(spikesel)
    data = varargin{spikesel(i)};
    for j=1:Ntrial(firstcontinuous)
      % the data will contain a one for each spike and be zero elsewhere
      data.time{j}  = varargin{firstcontinuous}.time{j};
      data.trial{j} = zeros(Nchan(spikesel(i)), length(data.time{j}));
    end
    for j=1:Nchan(spikesel(i))
      % for each spike in this channel, determine its corresponding sample in the original continuous data
      samplenum = round(double(data.timestamp{j} - FirstTimeStamp)/TimeStampPerSample) + 1;
      for k=1:size(trl{firstcontinuous},1)
        % determine which spikes are in this trial
        trlbeg = trl{firstcontinuous}(k,1);
        trlend = trl{firstcontinuous}(k,2);
        trlsmp = (samplenum(samplenum>=trlbeg & samplenum<=trlend) - trlbeg + 1);
        data.trial{k}(j,trlsmp) = data.trial{k}(j,trlsmp) + 1;  % increment with one 
      end
    end
    data = rmfield(data, 'waveform');
    data = rmfield(data, 'timestamp');
    data = rmfield(data, 'unit');
    varargin{spikesel(i)} = data;
    % update the local information belonging to this data set
    Ntrial(spikesel(i)) = Ntrial(firstcontinuous);
    trl{spikesel(i)}    = trl{firstcontinuous};
  end
end

if all(Nchan==length(unique(label)))
  % concatenate the trials
  fprintf('concatenating the trials over all datasets\n');
  data = varargin{1};
  data.trial  = {};
  data.time   = {};
  for i=1:Ndata
    data.trial  = cat(2, data.trial,  varargin{i}.trial(:)');
    data.time   = cat(2, data.time,   varargin{i}.time(:)');
  end
  % also concatenate the trial specification
  cfg.trl = cat(1, trl{:});

elseif length(unique(label))==sum(Nchan)
  % concatenate the channels in each trial
  fprintf('concatenating the channels within each trial\n');
  data = varargin{1};
  if ~all(diff(Ntrial)==0)
    error('not all datasets have the same number of trials')
  else
    Ntrial = Ntrial(1);
  end
  for i=2:Ndata
    for j=1:Ntrial
      if ~all(data.time{j}==varargin{i}.time{j})
        error('there is a difference in the time axes of the input data');
      end
      data.trial{j} = [data.trial{j}; varargin{i}.trial{j}];
    end
    data.label = [data.label(:); varargin{i}.label(:)];
  end

else
  % labels are inconsistent, cannot determine how to concatenate the data
  error('cannot determine how the data should be concatenated');
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: appenddata.m,v 1.12 2007/01/09 09:50:39 roboos Exp $';
% remember the configuration details of the input data
cfg.previous = [];
for i=1:Ndata
  try, cfg.previous{i} = varargin{i}.cfg; end
end
% remember the exact configuration details in the output 
data.cfg = cfg;

fprintf('output dataset, %d channels, %d trials\n', length(data.label), length(data.trial));

