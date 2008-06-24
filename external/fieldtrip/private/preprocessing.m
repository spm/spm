function [data] = preprocessing(cfg, data);

% PREPROCESSING reads MEG and/or EEG data according to user-specified trials
% and applies several user-specified preprocessing steps to the signals.
%
% Use as
%   [data] = preprocessing(cfg)
% or
%   [data] = preprocessing(cfg, data)
%
% The first input argument "cfg" is the configuration structure, which
% contains all details for the dataset filenames, trials and the
% preprocessing options. You can only do preprocessing after defining the
% segments of data to be read from the file (i.e. the trials), which is for
% example done based on the occurence of a trigger in the data.
%
% If you are calling PREPROCESSING with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see DEFINETRIAL
%   cfg.padding      = length to which the trials are padded for filtering
%   cfg.continuous   = 'yes' or 'no' whether the file contains continuous data
%                      (default is determined automatic)
%
% Instead of specifying the dataset, you can also explicitely specify the
% name of the file containing the header information and the name of the
% file containing the data, using
%   cfg.datafile     = string with the filename
%   cfg.headerfile   = string with the filename
%
% If you are calling PREPROCESSING with also the second input argument
% "data", then that should contain data that was already read from file in
% a previous call to PREPROCESSING. In that case only the configuration
% options below apply.
%
% The channels that will be read and/or preprocessed are specified with
%   cfg.channel      = Nx1 cell-array with selection of channels (default = 'all'),
%                      see CHANNELSELECTION for details
%
% The preprocessing options for the selected channels are specified with
%   cfg.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.bsfilter      = 'no' or 'yes'  bandstop filter
%   cfg.lnfilter      = 'no' or 'yes'  line noise removal using notch filter
%   cfg.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.lpfreq        = lowpass  frequency in Hz
%   cfg.hpfreq        = highpass frequency in Hz
%   cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.bsfreq        = bandstop frequency range, specified as [low high] in Hz
%   cfg.lnfreq        = line noise frequency in Hz, default 50Hz
%   cfg.dftfreq       = line noise frequencies for DFT filter, default [50 100 150] Hz
%   cfg.lpfiltord     = lowpass  filter order
%   cfg.hpfiltord     = highpass filter order
%   cfg.bpfiltord     = bandpass filter order
%   cfg.bsfiltord     = bandstop filter order
%   cfg.lnfiltord     = line noise notch filter order
%   cfg.lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.bsfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.lpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.hpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.bpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.bsfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   cfg.medianfiltord = length of median filter
%   cfg.blc           = 'no' or 'yes'
%   cfg.blcwindow     = [begin end] in seconds, the default is the complete trial
%   cfg.detrend       = 'no' or 'yes', this is done on the complete trial
%   cfg.polyremoval   = 'no' or 'yes', this is done on the complete trial
%   cfg.polyorder     = polynome order (default = 2)
%   cfg.hilbert       = 'no', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'no')
%   cfg.rectify       = 'no' or 'yes'
%   cfg.precision     = 'single' or 'double' (default = 'double')
%
% Preprocessing options that you should only use for EEG data are
%   cfg.reref         = 'no' or 'yes' (default = 'no')
%   cfg.refchannel    = cell-array with new EEG reference channel(s)
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.montage       = 'no' or a montage structure (default = 'no')
%
% Preprocessing options that you should only use when you are calling PREPROCESSING with 
% also the second input argument "data" are
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')

% Undocumented local options:
% cfg.artfctdef
% cfg.removemcg
%
% This function depends on PREPROC which has the following options:
% cfg.absdiff
% cfg.boxcar
% cfg.polyremoval, documented
% cfg.polyorder, documented
% cfg.blc, documented
% cfg.blcwindow, documented
% cfg.bpfilter, documented
% cfg.bpfiltord, documented
% cfg.bpfilttype, documented
% cfg.bpfreq, documented
% cfg.bsfilter, documented
% cfg.bsfiltord, documented
% cfg.bsfilttype, documented
% cfg.bsfreq, documented
% cfg.derivative
% cfg.detrend, documented
% cfg.dftfilter, documented
% cfg.dftfreq, documented
% cfg.hilbert, documented
% cfg.hpfilter, documented
% cfg.hpfiltord, documented
% cfg.hpfilttype, documented
% cfg.hpfreq, documented
% cfg.implicitref, documented
% cfg.lnfilter, documented
% cfg.lnfiltord, documented
% cfg.lnfilttype, documented
% cfg.lnfreq, documented
% cfg.lpfilter, documented
% cfg.lpfiltord, documented
% cfg.lpfilttype, documented
% cfg.lpfreq, documented
% cfg.medianfilter, documented
% cfg.medianfiltord, documented
% cfg.rectify, documented
% cfg.refchannel, documented
% cfg.reref, documented

% Copyright (C) 2003-2007, Robert Oostenveld, SMI, FCDC
%
% $Log: preprocessing.m,v $
% Revision 1.91  2008/06/24 12:47:12  roboos
% added cfg.montage as alternative for rereferencing
%
% Revision 1.90  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.89  2008/05/06 15:43:46  sashae
% change in trial selection, cfg.trials can be logical
%
% Revision 1.88  2008/03/04 16:37:04  roboos
% automatically read continuous or trial based data, this reads all data into memory
%
% Revision 1.87  2008/02/06 16:23:05  sashae
% added option for trial selection (only pertains to call of preprocessing with already preprocessed data)
%
% Revision 1.86  2007/11/21 23:57:02  roboos
% added default cfg.bsfilter=no
%
% Revision 1.85  2007/11/14 13:15:44  roboos
% fixed bug: also use filter padding for bandstopfilter, thanks to Saskia
%
% Revision 1.84  2007/11/08 10:09:50  roboos
% allow other versions of the hilbert transformed signal to be computed (e.g. complex, real, imag)
%
% Revision 1.83  2007/09/24 10:11:13  ingnie
% relevant fields of input data are copied to output data when preprocessing is done on data which is already read, try catch replaced be if isfield
%
% Revision 1.82  2007/09/20 12:58:23  roboos
% added possibility for preprocessing data that was already preprocessed (i.e. filter sequentially)
%
% Revision 1.81  2007/08/01 15:12:51  ingnie
% added documentation on polyremoval option
%
% Revision 1.80  2007/07/03 15:57:13  roboos
% added hack to support 32 bit Neuroscan cnt format, thanks to Nathan Weisz
%
% Revision 1.79  2007/05/30 13:24:43  roboos
% renamed the hidden option cfg.datatype to cfg.continuous=yes/no and added it to the documentation
%
% Revision 1.78  2007/05/02 15:23:07  roboos
% added option for bandstopfiltering
%
% Revision 1.77  2007/01/09 09:49:53  roboos
% allow numeric channel selections, give warning when cfg.rejectxxx is present (since deprecated)
%
% Revision 1.76  2007/01/04 17:07:37  roboos
% only call definetrial and rejectartifact if really needed, give a warning in that case (is deprecated)
%
% Revision 1.75  2006/11/29 09:06:36  roboos
% renamed all cfg options with "sgn" into "channel", added backward compatibility when required
% updated documentation, mainly in the artifact detection routines
%
% Revision 1.74  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.73  2006/08/31 07:56:05  roboos
% added onepass-reverse filter to documentation
%
% Revision 1.72  2006/06/14 12:51:01  roboos
% updated documentation
%
% Revision 1.71  2006/06/14 12:43:54  roboos
% removed the documentation for cfg.lnfilttype, since that option is not supported by preproc
%
% Revision 1.70  2006/06/14 11:52:06  roboos
% changed some comments, minor change to the handling of cfg.padding
%
% Revision 1.69  2006/06/13 14:48:09  ingnie
% updated documentation
%
% Revision 1.68  2006/05/30 20:57:15  roboos
% updated documentation
%
% Revision 1.67  2006/04/27 11:18:34  chrhes
% reinstated the default setting for the 'medianfilter' option since that was still broken, and changed a few words in the documentation.
%
% Revision 1.66  2006/04/25 20:20:46  roboos
% moved some of the sanity checks from preprocessing to private/preproc
% reinserted the default of some of the cfg settings, since that was broken
%
% Revision 1.65  2006/04/25 17:03:17  ingnie
% updated documentation and removed defaults that were also present in preproc.m from code
%
% Revision 1.64  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.63  2005/11/23 10:44:11  roboos
% added bp/lp/hp/lnfilttype to the documentation
%
% Revision 1.62  2005/09/14 07:47:11  roboos
% added support for single precision data as output
%
% Revision 1.61  2005/09/02 13:51:39  roboos
% added defaulf cfg.dftfreq=[50 100 150], added new option to documentation
%
% Revision 1.60  2005/08/05 09:16:22  roboos
% removed the obsolete data.offset
%
% Revision 1.59  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.58  2005/05/04 07:28:41  roboos
% changed fprintf into progress, added cfg.feedback default
%
% Revision 1.57  2005/04/13 08:48:38  roboos
% changed tab to space in a comment
%
% Revision 1.56  2005/01/25 13:59:57  roboos
% added check for cfg.datatype=continous and extended the call to read_fcdc_data with the boundary check for non-continous data
%
% Revision 1.55  2005/01/21 15:23:31  roboos
% added default for medianfilter and medianfilterord
%
% Revision 1.54  2005/01/21 09:53:11  roboos
% implemented median filter in preproc, updated help
%
% Revision 1.53  2005/01/10 16:52:11  roboos
% modifield the piece of code that does the actual preprocessing to use the
% private/preproc function. This ensures consistency with all other functions
% that perform filtering, etc.
%
% Revision 1.52  2004/12/06 12:43:03  roboos
% added warning for inconsistent EEG-rereferencing configuration
%
% Revision 1.51  2004/11/17 08:58:02  roboos
% re-ordered help and included description of cfg.detrend
% changed order of preprocessing: hilbert and rectify as last options
% removed the modification of the padding for integer number of line-noise cycles, this only worked if padding was unequal to 0
% the dftfilter function now takes care of ensuring that the sine wave is estimated on a integer number of line-noise cycles
%
% Revision 1.50  2004/11/15 09:17:14  roboos
% moved dataset to filename conversion into separate function (private/dataset2files)
%
% Revision 1.49  2004/10/01 09:53:45  roboos
% added one line of comments
%
% Revision 1.48  2004/09/22 10:20:27  roboos
% converted to use external subfunctions time2offset and offset2time
% and add offset field to data structure if it is missing
%
% Revision 1.47  2004/09/14 11:10:48  jansch
% fixed bug in sanity check concerning bpfiltering in conjunction with
% hilbert-transformation
%
% Revision 1.46  2004/08/20 06:57:19  roboos
% removed channel specific processing options for EMG and EEG, these are now applied to all channels
% updated and restructured help, added some extra configuration checks
%
% Revision 1.45  2004/08/05 08:58:51  roboos
% added reference to DEFINETRIAL in help
%
% Revision 1.44  2004/06/28 15:03:28  olejen
% removed empty line
%
% Revision 1.43  2004/06/28 15:02:01  roboos
% added empty line to demonstrate Ole
%
% Revision 1.42  2004/06/24 12:09:50  roberto
% removed old version information
%
% Revision 1.41  2004/06/24 10:13:43  roberto
% cosmetical changes in comments only
%
% Revision 1.40  2004/06/03 15:50:20  roberto
% removed experimental event handling, should be part of definetrial and further users own responsibility
%
% Revision 1.39  2004/05/27 15:48:56  roberto
% added cfg.removeeog and added handle for functino that should be implemented by Ali

% set the defaults
if ~isfield(cfg, 'channel'),      cfg.channel = {'all'};        end
if ~isfield(cfg, 'removemcg'),    cfg.removemcg = 'no';         end
if ~isfield(cfg, 'removeeog'),    cfg.removeeog = 'no';         end
if ~isfield(cfg, 'feedback'),     cfg.feedback = 'text';        end
if ~isfield(cfg, 'precision'),    cfg.precision = 'double';     end
if ~isfield(cfg, 'padding'),      cfg.padding = 0;              end % padding is only done when filtering

% these options relate to the actual preprocessing, it is neccessary to specify here because of padding
if ~isfield(cfg, 'lnfilter'),     cfg.lnfilter = 'no';          end
if ~isfield(cfg, 'dftfilter'),    cfg.dftfilter = 'no';         end
if ~isfield(cfg, 'lpfilter'),     cfg.lpfilter = 'no';          end
if ~isfield(cfg, 'hpfilter'),     cfg.hpfilter = 'no';          end
if ~isfield(cfg, 'bpfilter'),     cfg.bpfilter = 'no';          end
if ~isfield(cfg, 'bsfilter'),     cfg.bsfilter = 'no';          end
if ~isfield(cfg, 'medianfilter'), cfg.medianfilter = 'no';      end
% these options relate to the actual preprocessing, it is neccessary to specify here because of channel selection
if ~isfield(cfg, 'reref'),        cfg.reref = 'no';             end
if ~isfield(cfg, 'refchannel'),   cfg.refchannel = {};          end
if ~isfield(cfg, 'implicitref'),  cfg.implicitref = [];         end

% support for the following options was removed on 20 August 2004 in Revision 1.46
if isfield(cfg, 'emgchannel'), error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emghpfreq'),  error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emgrectify'), error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'emghilbert'), error('EMG specific preprocessing is not supported any more'); end
if isfield(cfg, 'eegchannel'), error('EEG specific preprocessing is not supported any more'); end
if isfield(cfg, 'resamplefs'), error('resampling is not supported any more, see RESAMPLEDATA'); end

if nargin>1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do preprocessing of data that has already been read
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % the input data must be raw
  data = checkdata(data, 'dataformat', 'raw', 'hasoffset', 'yes');

  if cfg.padding>0
    error('cfg.padding should be zero, since filter padding is only possible while reading the data from file');
  elseif isfield(cfg, 'trl')
    error('cfg.trl should not be present, since no data will be read');
  elseif isfield(cfg, 'dataset')
    error('cfg.dataset should not be present, since no data will be read');
  elseif isfield(cfg, 'datafile')
    error('cfg.datafile should not be present, since no data will be read');
  elseif isfield(cfg, 'headerfile')
    error('cfg.headerfile should not be present, since no data will be read');
  end
  
  % set the defaults
  if ~isfield(cfg, 'trials'), cfg.trials = 'all'; end

  % select trials of interest
  if ~strcmp(cfg.trials, 'all')
    if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
    fprintf('selecting %d trials\n', length(cfg.trials));
    data.trial  = data.trial(cfg.trials);
    data.time   = data.time(cfg.trials);
    data.offset = data.offset(cfg.trials);
  end

  % translate the channel groups (like 'all' and 'MEG') into real labels
  cfg.channel = channelselection(cfg.channel, data.label);
  rawindx = match_str(data.label, cfg.channel);

  % this will contain the newly processed data
  dataout = [];

  progress('init', cfg.feedback, 'preprocessing');
  ntrl = length(data.trial);
  for i=1:ntrl
    progress(i/ntrl, 'preprocessing trial %d from %d\n', i, ntrl);
    % do the preprocessing on the selected channels
    [dataout.trial{i}, dataout.label, dataout.time{i}, cfg] = preproc(data.trial{i}(rawindx,:), data.label(rawindx), data.fsample, cfg, data.offset(i));
  end % for all trials
  progress('close');

  % remember the configuration details of the input data
  if isfield(data, 'cfg'); cfg.previous = data.cfg; end
  
  % update the trial definition (trl) in case of trial selection
  if ~strcmp(cfg.trials, 'all')
    % try to locate the trl in the nested configuration
    if isfield(data, 'cfg')
      trl = findcfg(data.cfg, 'trl');
    else
      trl = [];
    end
    if isempty(trl)
      % a trial definition is expected in each continuous data set
      warning('could not locate the trial definition ''trl'' in the data structure');
    else
      cfg.trlold=trl;
      cfg.trl=trl(cfg.trials,:);
    end
  end
  
  % take along relevant fields of input data to output data
  if isfield(data, 'hdr');      dataout.hdr     = data.hdr;         end
  if isfield(data, 'fsample');  dataout.fsample = data.fsample;     end
  if isfield(data, 'grad');     dataout.grad    = data.grad;        end
  if isfield(data, 'elec');     dataout.elec    = data.elec;        end

  % replace the input data with the output data
  data = dataout;

else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the data from file and do the preprocessing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % if neccessary convert dataset into headerfile and datafile
  cfg = dataset2files(cfg);

  % check whether the data and header file are given
  if ~isfield(cfg, 'datafile') || ~isfield(cfg, 'headerfile')
    error('the datafile and/or headerfile is not correctly specified');
  end

  % read the header
  hdr = read_header(cfg.headerfile);

  % this option relates to reading over trial boundaries in a pseudo-continuous dataset
  if ~isfield(cfg, 'continuous')
    if isfield(cfg, 'datatype') && strcmp(cfg.datatype, 'continuous')
      cfg.continuous = 'yes';
    elseif hdr.nTrials==1
      cfg.continuous = 'yes';
    else
      cfg.continuous = 'no';
    end
  end

  % this should be a cell array
  if ~iscell(cfg.channel) && ischar(cfg.channel)
    cfg.channel = {cfg.channel};
  end

  % this should be a cell array
  if ~iscell(cfg.refchannel) && ischar(cfg.refchannel)
    cfg.refchannel = {cfg.refchannel};
  end

  % do a sanity check for the re-referencing
  if strcmp(cfg.reref, 'no') && ~isempty(cfg.refchannel)
    warning('no re-referencing is performed');
    cfg.refchannel = {};
  end

  % translate the channel groups (like 'all' and 'MEG') into real labels
  cfg.channel = channelselection(cfg.channel, hdr.label);

  if ~isempty(cfg.implicitref)
    % add the label of the implicit reference channel to these cell-arrays
    cfg.channel    = {cfg.channel{:} cfg.implicitref};
  end
  cfg.refchannel = channelselection(cfg.refchannel, cfg.channel);

  % determine the length in samples to which the data should be padded before filtering is applied
  % the filter padding is done by reading a longer segment of data from the original data file
  if cfg.padding>0
    if strcmp(cfg.lnfilter, 'yes') || ...
        strcmp(cfg.dftfilter, 'yes') || ...
        strcmp(cfg.lpfilter, 'yes') || ...
        strcmp(cfg.hpfilter, 'yes') || ...
        strcmp(cfg.bpfilter, 'yes') || ...
        strcmp(cfg.bsfilter, 'yes') || ...
        strcmp(cfg.medianfilter, 'yes')
      padding = round(cfg.padding * hdr.Fs);
    else
      % no filtering will be done, hence no padding is neccessary
      padding = 0;
    end
    % update the configuration (in seconds) for external reference
    cfg.padding = padding / hdr.Fs;
  else
    % no padding was requested
    padding = 0;
  end

  if ~isfield(cfg, 'trl')
    % treat the data as continuous if possible, otherwise define all trials as indicated in the header
    if strcmp(cfg.continuous, 'yes')
      trl = zeros(1, 3);
      trl(1,1) = 1;
      trl(1,2) = hdr.nSamples*hdr.nTrials;
      trl(1,3) = 0;
    else
      trl = zeros(hdr.nTrials, 3);
      for i=1:hdr.nTrials
        trl(i,1) = (i-1)*hdr.nSamples + 1;
        trl(i,2) = (i  )*hdr.nSamples    ;
        trl(i,3) = -hdr.nSamplesPre;
      end
    end
    cfg.trl = trl;
  end

  if any(strmatch('reject',       fieldnames(cfg))) || ...
      any(strmatch('rejecteog',    fieldnames(cfg))) || ...
      any(strmatch('rejectmuscle', fieldnames(cfg))) || ...
      any(strmatch('rejectjump',   fieldnames(cfg)))
    % this is only for backward compatibility
    error('you should call REJECTARTIFACT prior to PREPROCESSING, please update your scripts');
  end

  ntrl = size(cfg.trl,1);
  if ntrl<1
    error('no trials were selected for preprocessing, see DEFINETRIAL for help');
  end

  % compute the template for MCG and the QRS latency indices, and add it to the configuration
  if strcmp(cfg.removemcg, 'yes')
    cfg = template_mcg(cfg);
    mcgchannel = channelselection(cfg.artfctdef.mcg.channel, hdr.label);
    mcgindx    = match_str(cfg.channel, mcgchannel);
    for i=1:length(mcgchannel)
      fprintf('removing mcg on channel %s\n', mcgchannel{i});
    end
  end

  % determine the channel numbers of interest for preprocessing
  [chnindx, rawindx] = match_str(cfg.channel, hdr.label);

  progress('init', cfg.feedback, 'reading and preprocessing');
  for i=1:ntrl
    progress(i/ntrl, 'reading and preprocessing trial %d from %d\n', i, ntrl);
    % non-zero padding is used for filtering and line noise removal
    nsamples = cfg.trl(i,2)-cfg.trl(i,1)+1;
    if nsamples>padding
      % the trial is already longer than the total lenght requested
      begsample  = cfg.trl(i,1);
      endsample  = cfg.trl(i,2);
      begpadding = 0;
      endpadding = 0;
    else
      % begpadding+nsamples+endpadding = total length of raw data that will be read
      begpadding = ceil((padding-nsamples)/2);
      endpadding = floor((padding-nsamples)/2);
      begsample  = cfg.trl(i,1) - begpadding;
      endsample  = cfg.trl(i,2) + endpadding;
      if begsample<1
        warning('cannot apply enough padding at begin of file');
        begpadding = begpadding - (1 - begsample);
        begsample  = 1;
      end
      if endsample>(hdr.nSamples*hdr.nTrials)
        warning('cannot apply enough padding at end of file');
        endpadding = endpadding - (endsample - hdr.nSamples*hdr.nTrials);
        endsample  = hdr.nSamples*hdr.nTrials;
      end
    end

    % ONLY RELEVANT FOR NEUROSCAN CNT
    if ~isfield(cfg, 'nsdf')
      hdr.nsdf=16;
    else
      hdr.nsdf=cfg.nsdf;
    end

    % read the raw data with padding on both sides of the trial
    dat = read_data(cfg.datafile, hdr, begsample, endsample, rawindx, strcmp(cfg.continuous, 'yes'));

    % do the preprocessing on the padded trial data and remove the padding after filtering
    [cutdat{i}, label, time{i}, cfg] = preproc(dat, hdr.label(rawindx), hdr.Fs, cfg, cfg.trl(i,3), begpadding, endpadding);

    % ONLY RELEVANT FOR NEUROSCAN CNT
    hdr=rmfield(hdr,'nsdf');

  end % for all trials
  progress('close');

  % collect the results
  data.hdr                = hdr;                  % header details of the datafile
  data.label              = cfg.channel;          % labels of channels that have been read
  data.trial              = cutdat;               % cell-array with TIMExCHAN
  data.time               = time;                 % vector with the timeaxis for each individual trial
  data.fsample            = hdr.Fs;
  if isfield(hdr, 'grad')
    data.grad             = hdr.grad;             % gradiometer system in head coordinates
  end

end % if nargin>1

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: preprocessing.m,v 1.91 2008/06/24 12:47:12 roboos Exp $';

% remember the exact configuration details in the output
data.cfg = cfg;

