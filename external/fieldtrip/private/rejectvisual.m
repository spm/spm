function [data] = rejectvisual(cfg, data);

% REJECTVISUAL shows the preprocessed data in all channels and/or trials to
% allow the user to make a visual selection of the data that should be
% rejected. The data can be displayed in a "summary" mode, in which case
% the variance (or another metric) in each channel and each trial is
% computed. Alternatively, all channels can be shown at once allowing
% paging through the trials, or all trials can be shown, allowing paging
% through the channels.
%
% Use as
%   [data] = rejectvisual(cfg, data)
%
% The configuration can contain
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see CHANNELSELECTION for details
%   cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.latency     = [begin end] in seconds, or 'minperlength', 'maxperlength',
%                     'prestim', 'poststim' (default = 'maxperlength')
%   cfg.method      = string, describes how the data should be shown, this can be
%                     'summary'  show a single number for each channel and trial (default)
%                     'channel'  show the data per channel, all trials at once
%                     'trial'    show the data per trial, all channels at once
%   cfg.keepchannel = string, determines how to deal with channels that are
%                     not selected, can be
%                     'no'     completely remove unselected channels from the data (default)
%                     'yes'    keep unselected channels in the output data
%                     'nan'    fill the channels that are unselected with NaNs
%   cfg.metric      = string, describes the metric that should be computed in summary mode
%                     for each channel in each trial, can be
%                     'var'    variance within each channel (default)
%                     'min'    minimum value in each channel
%                     'max'    maximum value each channel
%                     'absmax' maximum absolute value in each channel
%                     'range'  range from min to max in each channel
%   cfg.alim        = value that determines the amplitude scaling for the
%                     channel and trial display, if empty then the amplitude
%                     scaling is automatic (default = [])
%
% Optionally, the raw data is preprocessed (filtering etc.) prior to
% displaying it or prior to computing the summary metric. The
% preprocessing and the selection of the latency window is NOT applied
% to the output data.
%
% The following settings are usefull for identifying EOG artifacts:
%   cfg.bpfilter    = 'yes'
%   cfg.bpfilttype  = 'but'
%   cfg.bpfreq      = [1 15]
%   cfg.bpfiltord   = 4
%   cfg.rectify     = 'yes'
%
% The following settings are usefull for identifying muscle artifacts:
%   cfg.bpfilter    = 'yes'
%   cfg.bpfreq      = [110 140]
%   cfg.bpfiltord   = 10
%   cfg.bpfilttype  = 'but'
%   cfg.rectify     = 'yes'
%   cfg.boxcar      = 0.2
%
% See also REJECTARTIFACT, REJECTCOMPONENT

% Undocumented local options:
% cfg.feedback
%
% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.blc
% cfg.blcwindow
% cfg.boxcar
% cfg.bpfilter
% cfg.bpfiltord
% cfg.bpfilttype
% cfg.bpfreq
% cfg.derivative
% cfg.detrend
% cfg.dftfilter
% cfg.dftfreq
% cfg.hilbert
% cfg.hpfilter
% cfg.hpfiltord
% cfg.hpfilttype
% cfg.hpfreq
% cfg.implicitref
% cfg.lnfilter
% cfg.lnfiltord
% cfg.lnfreq
% cfg.lpfilter
% cfg.lpfiltord
% cfg.lpfilttype
% cfg.lpfreq
% cfg.medianfilter
% cfg.medianfiltord
% cfg.rectify
% cfg.refchannel
% cfg.reref

% Copyright (C) 2005-2006, Markus Bauer, Robert Oostenveld
%
% $Log: rejectvisual.m,v $
% Revision 1.19  2007/12/18 17:52:20  sashae
% added option for trial selection, replaced some old code by call to findcfg
%
% Revision 1.18  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.17  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.16  2007/03/20 16:57:16  roboos
% added a fixme statement
%
% Revision 1.15  2007/01/18 10:08:18  roboos
% remove data.offset if present (old datasets)
%
% Revision 1.14  2007/01/11 13:53:13  roboos
% imlemented cfg.alim, which allows manual specificationh of the amplitude limits in the channel and trial display
%
% Revision 1.13  2007/01/10 11:46:01  roboos
% implemented selection of time window using cfg.latency
%
% Revision 1.12  2006/11/30 13:58:11  roboos
% implemented two new interactive methods: channel and trial browsing
%
% Revision 1.11  2006/11/07 08:29:56  roboos
% also allow trls with more than 3 columns
%
% Revision 1.10  2006/11/01 08:22:37  roboos
% made the default for cfg.channel consistent with documentation -> all instead of MEG
%
% Revision 1.9  2006/10/04 07:05:57  roboos
% added option keepchannel, can be no|yes|nan
%
% Revision 1.8  2006/08/29 14:26:27  roboos
% fixed bug: remove rejected trials from cfg.trl, and add the original trl as trlold to make it compatible with REJECTARTIFACT (thanks to Floris)
%
% Revision 1.7  2006/08/15 15:49:12  roboos
% implemented the deselection of channels in datasets, i.e. bad channels are removed (thanks to Markus Bauer)
%
%
% Revision 1.7  2006/06/14 12:44:00  marbau
% implemented the deselection of channels in datasets
%
% Revision 1.6  2006/06/14 12:44:00  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.5  2006/06/14 11:53:08  roboos
% switched to using cfg.preproc substructure
%
% Revision 1.4  2006/06/12 12:06:31  roboos
% added interactive option for plotting selected trials in a seperate figure (i.e. for detailled visual inspection)
%
% Revision 1.3  2006/06/12 11:20:30  roboos
% improved search for trl
%
% Revision 1.2  2006/06/12 07:51:10  roboos
% use local copy of offset, do not add to output data
%
% Revision 1.1  2006/05/17 14:45:42  roboos
% renamed rejecttrial into rejectvisual, added modified data (i.e. trials rejected) as output, implemented visual selection of bad trials using mouse, implemented manual selection of bad trials using keyboard
%
% Revision 1.2  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.1  2005/11/11 14:39:28  roboos
% new implementation based on old code from Markus
%

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

if ~isfield(cfg, 'channel'),     cfg.channel = 'all';          end
if ~isfield(cfg, 'trials'),       cfg.trials = 'all';          end
if ~isfield(cfg, 'latency'),     cfg.latency = 'maxperlength'; end
if ~isfield(cfg, 'keepchannel'), cfg.keepchannel = 'no';       end
if ~isfield(cfg, 'feedback'),    cfg.feedback = 'textbar';     end
if ~isfield(cfg, 'method'),      cfg.method = 'summary';       end
if ~isfield(cfg, 'alim'),        cfg.alim = [];                end

% for backward compatibility
if ~isfield(cfg, 'metric') && any(strcmp(cfg.method, {'var', 'min', 'max', 'absmax', 'range'}))
  cfg.metric = cfg.method;
  cfg.method = 'summary';
end

if ~isfield(cfg, 'metric')
  cfg.metric = 'var';
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
end

% determine the duration of each trial
for i=1:length(data.time)
  begsamplatency(i) = min(data.time{i});
  endsamplatency(i) = max(data.time{i});
end

% determine the latency window which is possible in all trials
minperlength = [max(begsamplatency) min(endsamplatency)];
maxperlength = [min(begsamplatency) max(endsamplatency)];

% latency window for averaging and variance computation is given in seconds
if (strcmp(cfg.latency, 'minperlength'))
  cfg.latency = [];
  cfg.latency(1) = minperlength(1);
  cfg.latency(2) = minperlength(2);
elseif (strcmp(cfg.latency, 'maxperlength'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = maxperlength(2);
elseif (strcmp(cfg.latency, 'prestim'))
  cfg.latency = [];
  cfg.latency(1) = maxperlength(1);
  cfg.latency(2) = 0;
elseif (strcmp(cfg.latency, 'poststim'))
  cfg.latency = [];
  cfg.latency(1) = 0;
  cfg.latency(2) = maxperlength(2);
end

% ensure that the preproc specific options are located in the cfg.preproc substructure
cfg = createsubcfg(cfg, 'preproc');

if strcmp(cfg.method, 'channel')
  fprintf('showing the data per channel, all trials at once\n');
  [chansel, trlsel, cfg] = rejectvisual_channel(cfg, data);
elseif strcmp(cfg.method, 'trial')
  fprintf('showing the data per trial, all channels at once\n');
  [chansel, trlsel, cfg] = rejectvisual_trial(cfg, data);
elseif strcmp(cfg.method, 'summary')
  fprintf('showing a summary for all channels and trials\n');
  [chansel, trlsel, cfg] = rejectvisual_summary(cfg, data);
end

fprintf('%d trials marked as GOOD, %d trials marked as BAD\n', sum(trlsel), sum(~trlsel));
fprintf('%d channels marked as GOOD, %d channels marked as BAD\n', sum(chansel), sum(~chansel));

% trl is not specified in the function call, but the data is given ->
% try to locate the trial definition (trl) in the nested configuration
if isfield(data, 'cfg')
  trl  = findcfg(data.cfg, 'trl');
else
  trl  = [];
end
if isempty(trl)
  % a trial definition is expected in each continuous data set
  warning('could not locate the trial definition ''trl'' in the data structure');
end
trlold=trl;
 
% construct an artifact matrix from the trl matrix
if ~isempty(trl)
  % remember the sample numbers (begin and end) of each trial and each artifact
  % updating the trl and creating a trlold makes it compatible with REJECTARTIFACT
  if ~strcmp(cfg.trials, 'all')
    trl=trl(cfg.trials,:);
  end
  cfg.artifact = trl(~trlsel,1:2);
  cfg.trl      = trl( trlsel,:);
  cfg.trlold   = trlold;
else
  % since sample numbers are unknown, it is not possible to remember them here
  cfg.artifact = [];
  cfg.trl      = [];
  cfg.trlold   = [];
end

% remove the selected trials from the data
data.time  = data.time(trlsel);
data.trial = data.trial(trlsel);

% remove the offset vector if present (only applies to datasets that have been preprocessed a long time ago)
if isfield(data, 'offset')
  data = rmfield(data, 'offset');
end


if ~all(chansel)
  switch cfg.keepchannel
    case 'no'
      % remove channels that are not selected
      for i=1:length(data.trial)
        data.trial{i} = data.trial{i}(chansel,:);
      end
      data.label = data.label(chansel);
    case 'nan'
      % fill the data from the bad channels with nans
      for i=1:length(data.trial)
        data.trial{i}(~chansel,:) = nan;
      end
    case 'yes'
      % keep all channels, also when they are not selected
  end
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
cfg.version.id = '$Id: rejectvisual.m,v 1.19 2007/12/18 17:52:20 sashae Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
data.cfg = cfg;

