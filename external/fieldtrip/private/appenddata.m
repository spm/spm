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
% Revision 1.15  2007/11/13 15:25:07  roboos
% removed handling of spike data, since that can be done in appendspike
%
% Revision 1.14  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.13  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
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

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i}, 'datatype', 'raw', 'feedback', 'no');
end

% set the defaults
cfg = [];

Ndata = nargin-1;
if Ndata<2
  error('you must give at least two datasets to append');
end

% determine the dimensions of the data
Nchan  = zeros(1,Ndata);
Ntrial = zeros(1,Ndata);
label  = {};
for i=1:Ndata
  Nchan(i) = length(varargin{i}.label);
  Ntrial(i) = length(varargin{i}.trial);
  fprintf('input dataset %d, %d channels, %d trials\n', i, Nchan(i), Ntrial(i));
  label = [label(:); varargin{i}.label(:)];
end

% try to locate the trial definition (trl) in the nested configuration
for i=1:Ndata
  if isfield(varargin{i}, 'cfg')
    trl{i} = findcfg(varargin{i}.cfg, 'trl');
  else
    trl{i} = [];
  end
  if isempty(trl{i})
    % a trial definition is expected in each continuous data set
    warning(sprintf('could not locate the trial definition ''trl'' in data structure %d', i));
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
cfg.version.id = '$Id: appenddata.m,v 1.15 2007/11/13 15:25:07 roboos Exp $';
% remember the configuration details of the input data
cfg.previous = [];
for i=1:Ndata
  try, cfg.previous{i} = varargin{i}.cfg; end
end
% remember the exact configuration details in the output 
data.cfg = cfg;

fprintf('output dataset, %d channels, %d trials\n', length(data.label), length(data.trial));

