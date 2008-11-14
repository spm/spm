function [data] = checkdata(data, varargin)

% CHECKDATA checks the input data of the main FieldTrip functions, e.g. whether
% the type of data strucure corresponds with the required data. If neccessary
% and possible, this function will adjust the data structure to the input
% requirements (e.g. change dimord, average over trials, convert inside from
% index into logical).
%
% If the input data does NOT correspond to the requirements, this function
% is supposed to give a elaborate warning message and if applicable point
% the user to external documentation (link to website).
%
% Use as
%   [data] = checkdata(data, ...)
%
% Optional input arguments should be specified as key-value pairs and can include
%   feedback           = yes, no
%   datatype           = raw, freq, timelock, comp, spike, source, volume, dip
%   dimord             = any combination of time, freq, chan, refchan, rpt, subj, chancmb, rpttap
%   senstype           = ctf151, ctf275, ctf151_planar, ctf275_planar, neuromag122, neuromag306, bti148, bti248, bti248_planar, magnetometer, electrode
%   ismeg              = yes, no
%   inside             = logical, index
%   hastrials          = yes, no
%   hasoffset          = yes, no (only applies to raw data)
%   hascumtapcnt       = yes, no (only applies to freq data)
%   hasdof             = yes, no
%   cmbrepresentation  = sparse, full (applies to covariance and cross-spectral density)
%
% For some options you can specify multiple values, e.g.
%   [data] = checkdata(data, 'megtype', {'ctf151', 'ctf275'}), e.g. in megrealign
%   [data] = checkdata(data, 'datatype', {'timelock', 'freq'}), e.g. in sourceanalysis

% Copyright (C) 2007-2008, Robert Oostenveld
%
% $Log: checkdata.m,v $
% Revision 1.2  2008/11/14 10:45:47  jansch
% extended the fixcsd subfunction
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.28  2008/09/15 13:24:17  roboos
% add dimord to source data, remove x/y/zgrid for volume and/or source, also if no conversion between the two is done
%
% Revision 1.27  2008/08/20 19:00:23  jansch
% also enable correct trial-handling in case of isvolume. before this fix,
% sourcestatistics crashed (reproduced with a grandaverage as input). more
% specifically, it crashed in setsubfield, because the trials were not handled
% correctly.
%
% Revision 1.26  2008/07/30 07:41:10  roboos
% when reshaping source parameters, also loop over source.trial(...) structure array and not only source.avg
%
% Revision 1.25  2008/07/25 12:52:44  roboos
% fixed bug, hasdof would remove cumtapcnt (thanks to Jurrian)
%
% Revision 1.24  2008/07/25 07:15:49  roboos
% cleaned up documentation and some minor aspects of the code (sprintf/error)
% added hasdof, with default value of 'no'
%
% Revision 1.23  2008/07/22 11:39:55  ingnie
% moved the "inside" representation code to prior to reshaping, otherwise the logical volume inside would not be reshaped properly
%
% Revision 1.22  2008/07/22 11:11:01  roboos
% changed the source2volume function so that it works with any type of source positions, as long as they span a 3D volume
%
% Revision 1.21  2008/07/21 11:01:47  roboos
% added source2volume conversion
% update isXXXdatatype after all possible conversions
% ensure that the input is either source or volume, not both
% ensure consistent data dimensions in case of volume data (reshape into 3D representation)
% ensure consistent data dimensions in case of source data (reshape into linear vector representation)
%
% Revision 1.20  2008/05/08 10:45:34  jansch
% changed function-call to sensortype into senstype. changed variable-name
% senstype into stype
%
% Revision 1.19  2008/04/01 10:47:05  jansch
% remove nans in reconstructed timecourses in freq2raw
%
% Revision 1.18  2008/01/29 16:05:18  sashae
% added default setting hasoffset='no'
%
% Revision 1.17  2008/01/18 15:41:11  sashae
% implemented hascumtapcnt, adjusted hasoffset
%
% Revision 1.16  2008/01/10 21:30:30  roboos
% changed some formatting, no functional change
%
% Revision 1.15  2007/12/12 10:33:59  roboos
% switched to local subfunction for the conversion from timelock to raw (instead of seperate data2raw function)
% added support for converting from raw to timelocked, with some constraints
%
% Revision 1.14  2007/11/23 14:40:16  ingnie
% Fixed bug in converting data that occured when datatype had multiple options.
% Added loop over all cells of datatype.
%
% Revision 1.13  2007/11/05 16:01:51  roboos
% fixed bug in error message for sensortype
%
% Revision 1.12  2007/07/25 07:30:06  ingnie
% added freq2raw functionality
%
% Revision 1.11  2007/07/03 16:11:35  roboos
% add fample if required for offset
%
% Revision 1.10  2007/05/30 13:27:00  roboos
% fixed bug for senstype, which did not support multiple options
%
% Revision 1.9  2007/05/30 12:00:24  ingnie
% fixed bug in senstype
%
% Revision 1.8  2007/05/29 16:55:56  ingnie
% added options 'senstype', 'ismeg', 'inside', 'hastrials' and 'hasoffset'
%
% Revision 1.7  2007/05/29 12:54:40  roboos
% updated documentation
%
% Revision 1.6  2007/05/16 11:30:06  roboos
% added check for dimord
%
% Revision 1.5  2007/05/03 07:15:05  roboos
% fixed bug: missing end in loop (thanks to Ali)
%
% Revision 1.4  2007/05/02 16:02:28  roboos
% added dim as a requirement to volume data
% changed the program flow (i.e. the locig) for converting data between different types
% implemented conversion between timelock and raw data
% implemented conversion between source and raw data (will only work for projected lcmv virtual-channel time courses)
%
% Revision 1.3  2007/04/18 10:21:16  roboos
% convert volume to source data if possible
%
% Revision 1.2  2007/04/05 15:35:08  roboos
% less strict on volume data
%
% Revision 1.1  2007/04/03 15:30:05  roboos
% renamed latest copy of checkinput to checkdata, to accomodate the upcoming checkcfg function
%
% Revision 1.7  2007/04/03 06:55:48  jansch
% fixed typo fourier -> fourierspctrm
%
% Revision 1.6  2007/04/02 12:05:02  jansch
% fixed typo
%
% Revision 1.5  2007/03/30 16:50:02  ingnie
% fixed bug; hasfreq also when data does not have time field
%
% Revision 1.4  2007/03/28 16:00:35  roboos
% rewrote most of the function together with Ingrid, implemented feedback and datatype checkdata
%
% Revision 1.3  2007/03/27 09:39:56  ingnie
% only fixdimord implemented
%
% Revision 1.2  2007/02/27 13:35:01  roboos
% added some comments for ingnie to work on
%
% Revision 1.1  2007/02/27 09:24:11  roboos
% initial version that only serves as placeholder for the documentation and some example code
%

% in case of an error this function could use dbstack for more detailled
% user feedback
%
% this function should replace/encapsulate
%   fixdimord
%   fixinside
%   fixprecision
%   fixvolume
%   data2raw
%   raw2data
%   grid2transform
%   transform2grid
%   fourier2crsspctrm
%   freq2cumtapcnt
%   sensortype
%   time2offset
%   offset2time
%
% other potential uses for this function:
%   time -> offset in freqanalysis
%   average over trials
%   csd as matrix

% get the optional input arguments
feedback      = keyval('feedback',      varargin); if isempty(feedback), feedback = 'no'; end
datatype      = keyval('datatype',      varargin);
dimord        = keyval('dimord',        varargin);
stype         = keyval('senstype',      varargin); % senstype is a function name which should not be masked
ismeg         = keyval('ismeg',         varargin);
inside        = keyval('inside',        varargin); % can be logical or index
hastrials     = keyval('hastrials',     varargin);
hasoffset     = keyval('hasoffset',     varargin); if isempty(hasoffset), hasoffset = 'no'; end
hascumtapcnt  = keyval('hascumtapcnt',  varargin);
hasdof        = keyval('hasdof',        varargin); if isempty(hasdof), hasdof = 'no'; end
cmbrepresentation = keyval('cmbrepresentation',  varargin);

% determine the type of input data
% this can be raw, freq, timelock, comp, spike, source, volume, dip
israw      = isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell');
isfreq     = isfield(data, 'label') && isfield(data, 'freq') && (isfield(data, 'powspctrm') || isfield(data, 'crsspctrm') || isfield(data, 'fourierspctrm'));
istimelock = isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && isfield(data, 'avg') && isnumeric(data.avg);
iscomp     = isfield(data, 'topo') || isfield(data, 'topolabel');
isspike    = isfield(data, 'label') && isfield(data, 'waveform') && isa(data.waveform, 'cell') && isfield(data, 'timestamp') && isa(data.timestamp, 'cell');
isvolume   = isfield(data, 'transform') && isfield(data, 'dim');
issource   = isfield(data, 'pos');
isdip      = isfield(data, 'dip');

if ~isequal(feedback, 'no')
  if israw
    nchan = length(data.label);
    ntrial = length(data.trial);
    fprintf('the input is raw data with %d channels and %d trials\n', nchan, ntrial);
  elseif isfreq
    nchan = length(data.label);
    nfreq = length(data.freq);
    if isfield(data, 'time'), ntime = num2str(length(data.time));, else ntime = 'no'; end
    fprintf('the input is freq data with %d channels, %d frequencybins and %s timebins\n', nchan, nfreq, ntime);
  elseif istimelock
    nchan = length(data.label);
    ntime = length(data.time);
    fprintf('the input is timelock data with %d channels and %d timebins\n', nchan, ntime);
  elseif iscomp
    ncomp = length(data.label);
    nchan = length(data.topolabel);
    fprintf('the input is component data with %d components and %d original channels\n', ncomp, nchan);
  elseif isspike
    nchan = length(data.label);
    fprintf('the input is spike data\n');
  elseif isvolume
    fprintf('the input is volume data with dimensions [%d %d %d]\n', data.dim(1), data.dim(2), data.dim(3));
  elseif issource
    nsource = size(data.pos, 1);
    fprintf('the input is source data with %d positions\n', nsource);
  elseif isdip
    fprintf('the input is dipole data\n');
  end
end % give feedback

if isfreq || istimelock || iscomp
  % ensure consistency between the dimord string and the axes that describe the data dimensions
  data = fixdimord(data);
end

if issource && isvolume
  % it should be either one or the other
  % the choise here is to represent it as volume description since that is simpler to handle
  % remove the unwanted fields
  if isfield(data, 'pos'),    data = rmfield(data, 'pos');    end
  if isfield(data, 'xgrid'),  data = rmfield(data, 'xgrid');  end
  if isfield(data, 'ygrid'),  data = rmfield(data, 'ygrid');  end
  if isfield(data, 'zgrid'),  data = rmfield(data, 'zgrid');  end
  issource = 0;
end

if ~isempty(datatype)
  if ~isa(datatype, 'cell')
    datatype = {datatype};
  end

  okflag = 0;
  for i=1:length(datatype)
    % check that the data matches with one or more of the required datatypes
    switch datatype{i}
      case 'raw'
        okflag = okflag + israw;
      case 'freq'
        okflag = okflag + isfreq;
      case 'timelock'
        okflag = okflag + istimelock;
      case 'comp'
        okflag = okflag + iscomp;
      case 'spike'
        okflag = okflag + isspike;
      case 'volume'
        okflag = okflag + isvolume;
      case 'source'
        okflag = okflag + issource;
      case 'dip'
        okflag = okflag + isdip;
    end % switch datatype
  end % for datatype

  if ~okflag
    % try to convert the data
    for iCell = 1:length(datatype)
      if isequal(datatype(iCell), {'source'}) && isvolume
        data = volume2source(data);
        isvolume = 0;
        issource = 1;
        okflag = 1;
      elseif isequal(datatype(iCell), {'volume'}) && issource
        data = source2volume(data);
        isvolume = 1;
        issource = 0;
        okflag = 1;
      elseif isequal(datatype(iCell), {'raw'}) && issource
        data = data2raw(data);
        issource = 0;
        israw = 1;
        okflag = 1;
      elseif isequal(datatype(iCell), {'raw'}) && istimelock
        data = timelock2raw(data);
        istimelock = 0;
        israw = 1;
        okflag = 1;
      elseif isequal(datatype(iCell), {'timelock'}) && israw
        data = raw2timelock(data);
        israw = 0;
        istimelock = 1;
        okflag = 1;
      elseif isequal(datatype(iCell), {'raw'}) && isfreq
        data = freq2raw(data);
        isfreq = 0;
        israw = 1;
        okflag = 1;
      end
    end % for iCell
  end % if okflag

  if ~okflag
    % construct an error message
    if length(datatype)>1
      str = sprintf('%s, ', datatype{1:(end-2)});
      str = sprintf('%s%s or %s', str, datatype{end-1}, datatype{end});
    else
      str = datatype{1};
    end
    str = sprintf('This function requires %s data as input.', str);
    error(str);
  end % if okflag
end

if ~isempty(dimord)
  if ~isa(dimord, 'cell')
    dimord = {dimord};
  end

  if isfield(data, 'dimord')
    okflag = any(strcmp(data.dimord, dimord));
  else
    okflag = 0;
  end

  if ~okflag
    % construct an error message
    if length(dimord)>1
      str = sprintf('%s, ', dimord{1:(end-2)});
      str = sprintf('%s%s or %s', str, dimord{end-1}, dimord{end});
    else
      str = dimord{1};
    end
    str = sprintf('This function requires data with a dimord of %s.', str);
    error(str);
  end % if okflag
end

if ~isempty(stype)
  if ~isa(stype, 'cell')
    stype = {stype};
  end

  if isfield(data, 'grad') || isfield(data, 'elec')
    if any(strcmp(senstype(data), stype));
      okflag = 1;
    else
      okflag = 0;
    end
  else
    okflag = 0;
  end

  if ~okflag
    % construct an error message
    if length(stype)>1
      str = sprintf('%s, ', stype{1:(end-2)});
      str = sprintf('%s%s or %s', str, stype{end-1}, stype{end});
    else
      str = stype{1};
    end
    str = sprintf('This function requires %s data as input, but you are giving %s data.', str, senstype(data));
    error(str);
  end % if okflag
end

if ~isempty(ismeg)
  if isequal(ismeg,'yes')
    okflag = isfield(data, 'grad');
  elseif isequal(ismeg,'no')
    okflag = ~isfield(data, 'grad');
  end

  if ~okflag && isequal(ismeg,'yes')
    error('This function requires MEG data with a ''grad'' field');
  elseif ~okflag && isequal(ismeg,'no')
    error('This function should not be given MEG data with a ''grad'' field');
  end % if okflag
end

if ~isempty(inside)
  % TODO absorb the fixinside function into this code
  data   = fixinside(data, inside);
  okflag = isfield(data, 'inside');

  if ~okflag
    % construct an error message
    error('This function requires data with an ''inside'' field.');
  end % if okflag
end

%if isvolume
%  % ensure consistent dimensions of the volumetric data
%  % reshape each of the volumes that is found into a 3D array
%  param = parameterselection('all', data);
%  dim   = data.dim;
%  for i=1:length(param)
%    tmp  = getsubfield(data, param{i});
%    tmp  = reshape(tmp, dim);
%    data = setsubfield(data, param{i}, tmp);
%  end
%end

if issource || isvolume,
  param = parameterselection('all', data);
  % ensure consistent dimensions of the source reconstructed data
  % reshape each of the volumes that is found into a linear vector
  if issource, dim   = [size(data.pos, 1) 1]; end
  % ensure consistent dimensions of the volumetric data
  % reshape each of the volumes that is found into a 3D array
  if isvolume, dim   = data.dim;              end
  for i=1:length(param)
    if any(param{i}=='.')
      % the parameter is nested in a substructure, which can have multiple elements (e.g. source.trial(1).pow, source.trial(2).pow, ...)
      % loop over the substructure array and reshape for every element
      tok  = tokenize(param{i}, '.');
      sub1 = tok{1};  % i.e. this would be 'trial'
      sub2 = tok{2};  % i.e. this would be 'pow'
      tmp1 = getfield(data, sub1);
      for j=1:numel(tmp1)
        tmp2 = getfield(tmp1(j), sub2);
        tmp2 = reshape(tmp2, dim);
        tmp1(j) = setfield(tmp1(j), sub2, tmp2);
      end
      data = setfield(data, sub1, tmp1);
    else
      tmp  = getfield(data, param{i});
      tmp  = reshape(tmp, dim);
      data = setfield(data, param{i}, tmp);
    end
  end
  if isfield(data, 'xgrid'),  data = rmfield(data, 'xgrid');  end
  if isfield(data, 'ygrid'),  data = rmfield(data, 'ygrid');  end
  if isfield(data, 'zgrid'),  data = rmfield(data, 'zgrid');  end
  if issource, data.dimord = 'pos'; end
end

if isequal(hastrials,'yes')
  okflag = isfield(data, 'trial');

  if ~okflag
    error('This function requires data with a ''trial'' field');
  end % if okflag
end

if isequal(hasoffset,'yes')
  okflag = isfield(data, 'offset');

  if ~okflag && isfield(data, 'time') && isa(data.time, 'cell')
    if ~isfield(data, 'fsample')
      data.fsample = 1/(data.time{1}(2)-data.time{1}(1));
    end
    for i=1:length(data.time);
      data.offset(i) = time2offset(data.time{i}, data.fsample);
    end
    okflag = 1;
  end

  if ~okflag
    error('This function requires data with an ''offset'' field');
  end % if okflag

elseif isequal(hasoffset,'no') && isfield(data, 'offset')
  data = rmfield(data, 'offset');
end % if hasoffset

if isequal(hascumtapcnt,'yes') && ~isfield(data, 'cumtapcnt')
  error('This function requires data with a ''cumtapcnt'' field');
elseif isequal(hascumtapcnt,'no') && isfield(data, 'cumtapcnt')
  data = rmfield(data, 'cumtapcnt');
end % if hascumtapcnt

if isequal(hasdof,'yes') && ~isfield(data, 'hasdof')
  error('This function requires data with a ''dof'' field');
elseif isequal(hasdof,'no') && isfield(data, 'hasdof')
  data = rmfield(data, 'cumtapcnt');
end % if hasdof

if ~isempty(cmbrepresentation)
  if istimelock
    data = fixcov(data, cmbrepresentation);
  elseif isfreq
    data = fixcsd(data, cmbrepresentation);
  else
    error('This function requires data with a covariance or cross-spectrum');
  end
end % cmbrepresentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% represent the covariance matrix in a particular manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = fixcov(data, desired)
if isfield(data, 'cov') && ~isfield(data, 'labelcmb')
  current = 'full';
elseif isfield(data, 'cov') && isfield(data, 'labelcmb')
  current = 'sparse';
else
  error('Could not determine the current representation of the covariance matrix');
end
if isequal(current, desired)
  % nothing to do
elseif strcmp(current, 'full') && strcmp(desired, 'sparse')
  % FIXME should be implemented
elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
  % FIXME should be implemented
end
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% represent the cross-spectral-density matrix in a particular manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = fixcsd(data, desired, channelcmb)
if isfield(data, 'powspctrm')
  current = 'sparsewithpow';
elseif isfield(data, 'crsspctrm') && ~isfield(data, 'labelcmb')
  current = 'full';
elseif isfield(data, 'crsspctrm') && isfield(data, 'labelcmb')
  current = 'sparse';
elseif isfield(data, 'fourierspctrm') && ~isfield(data, 'labelcmb')
  current = 'fourier';
else
  error('Could not determine the current representation of the cross-spectrum matrix');
end
if isequal(current, desired)
  % nothing to do
elseif (strcmp(current, 'full')    && strcmp(desired, 'sparsewithpow')) || ...
       (strcmp(current, 'fourier') && strcmp(desired, 'sparsewithpow'))

elseif (strcmp(current, 'sparse')  && strcmp(desired, 'sparsewithpow'))
  % convert back to crsspctrm/powspctrm representation: useful for plotting functions etc


elseif (strcmp(current, 'full')          && strcmp(desired, 'fourier')) || ...
 (strcmp(current, 'sparse')        && strcmp(desired, 'fourier')) || ...
 (strcmp(current, 'sparsewithpow') && strcmp(desired, 'fourier'))
  % this is not possible

elseif strcmp(current, 'full')          && strcmp(desired, 'sparse')
  % why would you want this? FIXME give explicit error
elseif strcmp(current, 'fourier')       && strcmp(desired, 'sparse')

  if isempty(channelcmb), error('no channel combinations are specified'); end
  % this is what freqdescriptives currently does as an intermediate step
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)), 
    nrpt = length(data.cumtapcnt);
    flag = 0;
  else 
    nrpt = 1;
    flag = 1;
  end 
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end
  ncmb      = size(channelcmb,1);
  cmbindx   = zeros(ncmb,2);
  labelcmb  = cell(ncmb,2);
  for k = 1:ncmb
    ch1 = find(strcmp(data.label, channelcmb(k,1)));    
    ch2 = find(strcmp(data.label, channelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2),
      cmbindx(k,:)  = [ch1 ch2];
      labelcmb(k,:) = data.label([ch1 ch2])';
    end    
  end
  
  crsspctrm = zeros(nrpt,ncmb,nfrq,ntim)+i.*zeros(nrpt,ncmb,nfrq,ntim);
  sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
  for k = 1:ntim
    for m = 1:nfrq
      for p = 1:nrpt
        indx    = (sumtapcnt(p)+1):sumtapcnt(p+1);
        tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),m,k);
        tmpdat2 = data.fourierspctrm(indx,cmbindx(:,2),m,k);
	crsspctrm(p,:,m,k) = (sum(tmpdat1.*conj(tmpdat2),1))./data.cumtapcnt(p);
      end
    end
  end
  data.crsspctrm = crsspctrm;
  data.labelcmb  = labelcmb;
  data           = rmfield(data, 'fourierspctrm');
  if nrpt>1, 
    if ntim>1,
      data.dimord = 'rpt_chan_freq_time';
    else  
      data.dimord = 'rpt_chan_freq';
    end
  else
    if ntim>1,
      data.dimord = 'chan_freq_time';
    else
      data.dimord = 'chan_freq';
    end
  end
  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end
  

elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'sparse')

  % this is what freqdescriptives currently does as an intermediate step,
  % and will be the new default representation for sparse data i.e. autospectra 
  % as {'A' 'A'} in labelcmb
  if isfield(data, 'crsspctrm'),
    dimtok         = tokenize(data.dimord, '_');
    catdim         = match_str(dimtok, {'chan' 'chancmb'});
    data.crsspctrm = cat(catdim, data.powspctrm, data.crsspctrm);
    data.labelcmb  = [data.label(:) data.label(:); data.labelcmb];
    data           = rmfield(data, 'powspctrm');
  else
    data.crsspctrm = data.powspctrm;
    data.labelcmb  = [data.label(:) data.label(:)];
    data           = rmfield(data, 'powspctrm');
  end

elseif strcmp(current, 'fourier') && strcmp(desired, 'full')

  % this is how it is currently and the desired functionality of prepare_freq_matrices
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)), 
    nrpt = length(data.cumtapcnt);
    flag = 0;
  else 
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('rpttap',dimtok)), nrpt=length(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end
  nchan     = length(data.label);
  crsspctrm = zeros(nrpt,nchan,nchan,nfrq,ntim)+i.*zeros(nrpt,nchan,nchan,nfrq,ntim);
  sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
  for k = 1:ntim
    for m = 1:nfrq
      for p = 1:nrpt
        indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
        tmpdat = transpose(data.fourierspctrm(indx,:,m,k));
	crsspctrm(p,:,:,m,k) = (tmpdat*tmpdat')./data.cumtapcnt(p);
      end
    end
  end
  data.crsspctrm = crsspctrm;
  data           = rmfield(data, 'fourierspctrm');
  if nrpt>1, 
    if ntim>1,
      data.dimord = 'rpt_chan_chan_freq_time';
    else  
      data.dimord = 'rpt_chan_chan_freq';
    end
  else
    if ntim>1,
      data.dimord = 'chan_chan_freq_time';
    else
      data.dimord = 'chan_chan_freq';
    end
  end
  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end

elseif strcmp(current, 'sparse') && strcmp(desired, 'full')

  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), 
    nrpt = size(data.crsspctrm,1);
    flag = 0;
  else 
    nrpt = 1;
    data.crsspctrm = reshape(data.crsspctrm, [1 size(data.crsspctrm)]);
    flag = 1;
  end 
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end
  nchan     = length(data.label);
  cmbindx   = zeros(nchan);
  for k = 1:size(data.labelcmb,1)
    ch1 = find(strcmp(data.label, data.labelcmb(k,1)));    
    ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2), cmbindx(ch1,ch2) = k; end
  end
  crsspctrm = zeros(nrpt,nchan,nchan,nfrq,ntim)+i.*zeros(nrpt,nchan,nchan,nfrq,ntim);
  for k = 1:ntim
    for m = 1:nfrq
      for p = 1:nrpt
        tmpdat = nan+zeros(nchan);
        tmpdat(find(cmbindx)) = data.crsspctrm(p,cmbindx(find(cmbindx)),m,k);
	tmpdat                = ctranspose(tmpdat);
	tmpdat(find(cmbindx)) = data.crsspctrm(p,cmbindx(find(cmbindx)),m,k);
	crsspctrm(p,:,:,m,k)  = tmpdat;
      end
    end
  end
  data.crsspctrm = crsspctrm;
  data           = rmfield(data, 'labelcmb');
  if nrpt>1, 
    if ntim>1,
      data.dimord = 'rpt_chan_chan_freq_time';
    else  
      data.dimord = 'rpt_chan_chan_freq';
    end
  else
    if ntim>1,
      data.dimord = 'chan_chan_freq_time';
    else
      data.dimord = 'chan_chan_freq';
    end
  end
  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end
  
elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'full')
  
  % this is how is currently done in prepare_freq_matrices
  data = checkdata(data, 'cmbrepresentation', 'sparse');
  data = checkdata(data, 'cmbrepresentation', 'full');

end
keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = volume2source(data)
xgrid = 1:data.dim(1);
ygrid = 1:data.dim(2);
zgrid = 1:data.dim(3);
[x y z] = ndgrid(xgrid, ygrid, zgrid);
data.pos = warp_apply(data.transform, [x(:) y(:) z(:)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = source2volume(data)
% the volume representation requires a homogenous transformation matrix to convert voxel indices to head coordinates
xgrid = 1:data.dim(1);
ygrid = 1:data.dim(2);
zgrid = 1:data.dim(3);
[x y z] = ndgrid(xgrid, ygrid, zgrid);
ind =  [x(:) y(:) z(:)];    % these are the positions expressed in voxel indices along each of the three axes
pos = data.pos;             % these are the positions expressed in head coordinates	
% represent the positions in a manner that is compatible with the homogenous matrix multiplication, i.e. pos = H * ind
ind = ind'; ind(4,:) = 1;
pos = pos'; pos(4,:) = 1;
% recompute the homogenous transformation matrix
data.transform = pos / ind;
% remove the unwanted fields
if isfield(data, 'pos'),    data = rmfield(data, 'pos');    end
if isfield(data, 'xgrid'),  data = rmfield(data, 'xgrid');  end
if isfield(data, 'ygrid'),  data = rmfield(data, 'ygrid');  end
if isfield(data, 'zgrid'),  data = rmfield(data, 'zgrid');  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = freq2raw(freq)
if ~strcmp(freq.dimord, 'rpt_chan_freq_time')
  error('this only works for dimord=''rpt_chan_freq_time''');
end
nrpt  = size(freq.powspctrm,1);
nchan = size(freq.powspctrm,2);
nfreq = size(freq.powspctrm,3);
ntime = size(freq.powspctrm,4);
data = [];
% create the channel labels like "MLP11@12Hz""
k = 0;
for i=1:nfreq
  for j=1:nchan
    k = k+1;
    data.label{k} = sprintf('%s@%dHz', freq.label{j}, freq.freq(i));
  end
end
% reshape and copy the data as if it were timecourses only
for i=1:nrpt
  data.time{i}  = freq.time;
  data.trial{i} = reshape(freq.powspctrm(i,:,:,:), nchan*nfreq, ntime);
  if any(isnan(data.trial{i}(1,:))),
    tmp = data.trial{i}(1,:);
    begsmp = find(isfinite(tmp),1,'first');
    endsmp = find(isfinite(tmp),1,'last' );
    data.trial{i} = data.trial{i}(:, begsmp:endsmp);
    data.time{i}  = data.time{i}(begsmp:endsmp);
  end
end
data.fsample = 1/(data.time{1}(2)-data.time{1}(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = raw2timelock(data)
ntrial = numel(data.trial);
nchan  = numel(data.label);
if ntrial==1
  data.time   = data.time{1};
  data.avg    = data.trial{1};
  data        = rmfield(data, 'trial');
  data.dimord = 'chan_time';
else
  nsample  = size(data.trial{1},2);
  tmptrial = zeros(ntrial, nchan, nsample);
  for i=1:ntrial
    if ~isequal(data.time{i}, data.time{1})
      error('this conversion is only supported for trials of equal length and not shifted relative to each other');
    else
      tmptrial(i,:,:) = data.trial{i};
    end
  end
  data.avg     = reshape(mean(tmptrial, 1), nchan, nsample);
  data.trial   = tmptrial;
  data.time    = data.time{1};
  data.dimord = 'rpt_chan_time';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert between datatypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = timelock2raw(data)
switch data.dimord
  case 'chan_time'
    data.trial{1} = data.avg;
    data.time     = {data.time};
    data          = rmfield(data, 'avg');
  case 'rpt_chan_time'
    tmptrial = {};
    tmptime  = {};
    ntrial = size(data.trial,1);
    nchan  = size(data.trial,2);
    ntime  = size(data.trial,3);
    for i=1:ntrial
      tmptrial{i} = reshape(data.trial(i,:,:), [nchan, ntime]);
      tmptime{i}  = data.time;
    end
    data       = rmfield(data, 'trial');
    data.trial = tmptrial;
    data.time  = tmptime;
  case 'subj_chan_time'
    tmptrial = {};
    tmptime  = {};
    ntrial = size(data.individual,1);
    nchan  = size(data.individual,2);
    ntime  = size(data.individual,3);
    for i=1:ntrial
      tmptrial{i} = reshape(data.individual(i,:,:), [nchan, ntime]);
      tmptime{i}  = data.time;
    end
    data       = rmfield(data, 'individual');
    data.trial = tmptrial;
    data.time  = tmptime;
  otherwise
    error('unsupported dimord');
end
% remove the unwanted fields
if isfield(data, 'avg'),        data = rmfield(data, 'avg'); end
if isfield(data, 'var'),        data = rmfield(data, 'var'); end
if isfield(data, 'cov'),        data = rmfield(data, 'cov'); end
if isfield(data, 'dimord'),     data = rmfield(data, 'dimord'); end
if isfield(data, 'numsamples'), data = rmfield(data, 'numsamples'); end
if isfield(data, 'dof'),        data = rmfield(data, 'dof'); end

