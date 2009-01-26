function [data] = selectdata(varargin)

% this function serves to concatenate the input data-structures along the
% compatible dimensions and thus is a more general implementation of
% appenddata, which only raw data. Moreover, it can be used to equate the
% data of different conditions to match e.g. in channels time-axis etc
% Moreover, it can be used as a generalization to ...average with 'keepindividual'
%
% Finally, this function serves to subselect regions-of-interest from the input data,
% either or not averaging across the specified dimensions
% supported input: -freq
%                 -timelock (not yet)
%                 -source   (not yet)
%                 -volume   (not yet)
%
% supported options: -foilim
%                   -latency
%                   -roi
%                   -channel (FIXME this is also done by preprocessing?)
%
%                   -avgoverchan
%                   -avgoverfreq
%                   -avgovertime
%                   -avgoverroi
%                   -avgoverrpt

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% $Log: selectdata.m,v $
% Revision 1.4  2009/01/26 20:53:11  roboos
% ensure that some options are true|false
% changes some whitespace
%
% Revision 1.3  2009/01/12 17:05:58  roboos
% fixed some whitespace
%

% check the input data and options
isdata  = find(cellfun(@isstruct,varargin));
keyvals = setdiff(1:length(varargin),isdata);

data = varargin(isdata);
kvp  = varargin(keyvals);
for k = 1:length(data)
  data{k} = checkdata(data{k}, 'datatype', {'freq' 'timelock' 'source', 'volume'});
  [dtype{k}, dimord{k}]  = datatype(data{k});
end

if ~all(strcmp(dtype{1},dtype)),
  error('different types of input-data is not supported');
end

isfreq   = datatype(data{1},'freq');
istlck   = datatype(data{1},'timelock');
issource = datatype(data{1},'source');
isvolume = datatype(data{1},'volume');

selchan  = keyval('channel', kvp); selectchan = ~isempty(selchan);
selfoi   = keyval('foilim',  kvp); selectfoi  = ~isempty(selfoi);
seltime  = keyval('latency', kvp); selecttime = ~isempty(seltime);
selroi   = keyval('roi',     kvp); selectroi  = ~isempty(selroi);
selrpt   = keyval('rpt',     kvp); selectrpt  = ~isempty(selrpt);
param    = keyval('param',   kvp); if isempty(param), param = 'all'; end

avgoverchan  = keyval('avgoverchan',  kvp); if isempty(avgoverchan), avgoverchan = false; end
avgoverfreq  = keyval('avgoverfreq',  kvp); if isempty(avgoverfreq), avgoverchan = false; end
avgovertime  = keyval('avgovertime',  kvp); if isempty(avgovertime), avgoverchan = false; end
avgoverroi   = keyval('avgoverroi',   kvp); if isempty(avgoverroi),  avgoverroi  = false; end
avgoverrpt   = keyval('avgoverrpt',   kvp); if isempty(avgoverrpt),  avgoverrpt  = false; end

% create anonymous function and apply it to the boolean input arguments
istrue = @(x)(ischar(x) && (strcmpi(x, 'yes') || strcmpi(x, 'true')) || x==1);
avgoverchan = istrue(avgoverchan);
avgoverfreq = istrue(avgoverfreq);
avgovertime = istrue(avgovertime);
avgoverroi  = istrue(avgoverroi);
avgoverrpt  = istrue(avgoverrpt);

if length(data)>1 && selectrpt,
  error('multiple data-structure and subselection of trials is not supported');
end

% check consistency of input data
if any(~strmatch(dimord{1},dimord))
  error('inconsistent dimord');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from here on the data is concatenated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(data)>1,
  % determine the way to concatenate

  if issource || isvolume,
    param = parameterselection(param, data{1}); % FIXME check consistency across input data of presence of specific parameters
  else
    param = {param};
  end

  dimtok = tokenize(dimord{1}, '_');
  dimtok(strmatch('chan', dimtok)) = {'label'}; % data.chan does not exist

  dimmat = zeros(length(dimtok), length(data));
  dimmat(:,1) = 1;
  for k = 1:length(dimtok)
    try,
      dimdat = getfield(data{1}, dimtok{k});
    catch
      % dimtok is probably 'rpt' or so
      dimdat = getsubfield(data{1}, param{1});
    end
    for m = 2:length(data)
      try,
        dimdat2 = getfield(data{m},dimtok{k});
      catch
        % dimtok is probably 'rpt' or so
        dimdat2 = getsubfield(data{m}, param{1});
      end
      try, dimmat(k,m) = all(dimdat(:)==dimdat2(:));            catch end;
      try, dimmat(k,m) = all(cellfun(@isequal,dimdat,dimdat2)); catch end;
    end
  end
  catdim = find(sum(dimmat,2)<length(data));

  if length(catdim)>1,
    error('ambiguous dimensions for concatenation');
  end

  % concatenate the data
  for k = 1:length(param)
    tmp = cell(1,length(data));
    for m = 1:length(tmp)
      tmp{m} = getsubfield(data{m},param{k});
    end
    data{1} = setsubfield(data{1}, param{k}, cat(catdim,tmp{:}));
  end

  % concatenate the relevant descriptive fields in the data-structure
  if ~strcmp(dimtok{catdim},'rpt') && ~strcmp(dimtok{catdim},'rpttap'),
    for k = 1:length(data)
      if k==1,
        tmp = getsubfield(data{k}, dimtok{catdim});
      else
        if strcmp(dimtok{catdim},'pos'),
          tmp = [tmp; getsubfield(data{k}, dimtok{catdim})]; sortflag = 0;
        else
          tmp = [tmp  getsubfield(data{k}, dimtok{catdim})]; sortflag = 1;
        end
      end
    end
    data{1} = setsubfield(data{1}, dimtok{catdim}, tmp);
  else
    % no such field as {'label','time','freq','pos'} has to be concatenated
    sortflag = 0;
  end

  % concatenate the relevant descriptive fields in the data-structure (continued)
  tryfields = {'cumsumcnt' 'cumtapcnt' 'dof'};
  for k = 1:length(tryfields)
    try,
      for m = 1:length(data)
        if m==1,
          tmpfield = getfield(data{m}, tryfields{k});
        else
          tmpfield = [tmpfield; getfield(data{m}, tryfields{k})];
        end
      end
      data{1} = setfield(data{1}, tryfields{k}, tmpfield);
    catch
    end
  end

  % FIXME this is ugly: solve it
  if issource || isvolume,
    data{1}.dim(catdim) = max(size(tmp));
  end

  % sort concatenated data FIXME this is also ugly and depends on tmp
  if sortflag && ~iscell(tmp),
    [srt, ind] = sort(tmp, 2);
    data{1} = setsubfield(data{1}, dimtok{catdim}, tmp(ind));
    for k = 1:length(param)
      tmp     = getsubfield(data{1}, param{k});
      tmp     = permute(tmp, [catdim setdiff(1:length(size(tmp)), catdim)]);
      tmp     = ipermute(tmp(ind,:,:,:,:), [catdim setdiff(1:length(size(tmp)), catdim)]);
      data{1} = setsubfield(data{1}, param{k}, tmp);
    end
  elseif iscell(tmp)
    %in this case (ugly!) tmp is probably a cell-array containing functional data
  end

  % remove unspecified parameters
  rmparam = setdiff(parameterselection('all',data{1}),param);
  for k = 1:length(rmparam)
    data{1} = rmsubfield(data{1}, rmparam{k});
  end

  % keep the first structure only
  data = data{1};

else
  % nothing to do
  data = data{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from here on a subset is selected from the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if selectrpt,
  if ~isfreq,
    % do nothing
  else
    dimtok = tokenize(data.dimord, '_');
    if strcmp(dimtok{1}, 'rpttap'),
      error('here you have to ensure the correct handling of tapers');
    end
  end
end

if selectchan,
  selchan = match_str(data.label, channelselection(selchan, data.label));
end

if selectfoi,
  if length(selfoi)==1, selfoi(2) = selfoi; end;
  selfoi = nearest(data.freq, selfoi(1)):nearest(data.freq, selfoi(2));
end

if selecttime,
  if length(seltime)==1, seltime(2) = seltime; end;
  seltime = nearest(data.time, seltime(1)):nearest(data.time, seltime(2));
end

if isfreq,

  dimtok  = tokenize(data.dimord, '_');
  chandim = strmatch('chan', dimtok);

  % hard-coded assumption is that the order in dimord is always chan first (could be 1 or 2), then freq, then time
  if isfield(data, 'fourierspctrm')
    if chandim==1, % FIXME this is actually impossible
      if selectchan, data.fourierspctrm = data.fourierspctrm(selchan,:,:); end
      if selectfoi,  data.fourierspctrm = data.fourierspctrm(:,selfoi,:);  end
      if selecttime, data.fourierspctrm = data.fourierspctrm(:,:,seltime); end
    elseif chandim==2,
      if selectrpt,  data.fourierspctrm = data.fourierspctrm(selrpt,:,:,:);  end
      if selectchan, data.fourierspctrm = data.fourierspctrm(:,selchan,:,:); end
      if selectfoi,  data.fourierspctrm = data.fourierspctrm(:,:,selfoi,:);  end
      if selecttime, data.fourierspctrm = data.fourierspctrm(:,:,:,seltime); end
    end
  end
  if isfield(data, 'powspctrm')
    if chandim==1,
      if selectchan, data.powspctrm = data.powspctrm(selchan,:,:); end
      if selectfoi,  data.powspctrm = data.powspctrm(:,selfoi,:);  end
      if selecttime, data.powspctrm = data.powspctrm(:,:,seltime); end
    elseif chandim==2,
      if selectrpt,  data.powspctrm = data.powspctrm(selrpt,:,:,:);  end
      if selectchan, data.powspctrm = data.powspctrm(:,selchan,:,:); end
      if selectfoi,  data.powspctrm = data.powspctrm(:,:,selfoi,:);  end
      if selecttime, data.powspctrm = data.powspctrm(:,:,:,seltime); end
    end
  end
  if isfield(data, 'crsspctrm')
    if chandim==1,
      if selectchan, data.crsspctrm = data.crsspctrm(selchan,:,:); end
      if selectfoi,  data.crsspctrm = data.crsspctrm(:,selfoi,:);  end
      if selecttime, data.crsspctrm = data.crsspctrm(:,:,seltime); end
    elseif chandim==2,
      if selectrpt,  data.crsspctrm = data.crsspctrm(selrpt,:,:,:);  end
      if selectchan, data.crsspctrm = data.crsspctrm(:,selchan,:,:); end
      if selectfoi,  data.crsspctrm = data.crsspctrm(:,:,selfoi,:);  end
      if selecttime, data.crsspctrm = data.crsspctrm(:,:,:,seltime); end
    end
  end

  if selectrpt,  data.cumtapcnt = data.cumtapcnt(selrpt); end % FIXME mtconvol
  if selectchan, data.label     = data.label(selchan);    end
  if selectfoi,  data.freq      = data.freq(selfoi);      end
  if selecttime, data.time      = data.time(seltime);     end

  if avgoverrpt,  error('not yet implemented'); end
  if avgoverchan, error('not yet implemented'); end
  if avgoverfreq, error('not yet implemented'); end
  if avgovertime, error('not yet implemented'); end

elseif istlck,
  if selecttime && isfield(data, 'cov'),
    error('it is not possible to extract a latency window here, please re-run timelockanalysis');
  end
  if isfield(data, 'trial'),
    hastrl = 1;
    if selectrpt,  data.trial = data.trial(selrpt,:,:);  end;
    if selectchan, data.trial = data.trial(:,selchan,:); end;
    if selecttime, data.trial = data.trial(:,:,seltime); end;
  else
    hastrl = 0;
  end
  if isfield(data, 'avg'),
    if selectrpt,  data     = rmfield(data, 'avg'); end;
    if selectchan, data.avg = data.avg(selchan,:);  end;
    if selecttime, data.avg = data.avg(:,seltime);  end;
  end
  if isfield(data, 'var'),
    if selectrpt,  data     = rmfield(data, 'var'); end;
    if selectchan, data.var = data.var(selchan,:);  end;
    if selecttime, data.var = data.var(:,seltime);  end;
  end
  if isfield(data, 'cov'),
    if hastrl,
      if selectrpt,  data.cov = data.cov(selrpt,:,:);        end;
      if selectchan, data.cov = data.cov(:,selchan,selchan); end;
    else
      if selectchan, data.cov = data.cov(selchan,selchan);   end;
    end
  end
  if isfield(data, 'blcov'),
    if hastrl,
      if selectrpt,  data.blcov = data.blcov(selrpt,:,:);        end;
      if selectchan, data.blcov = data.blcov(:,selchan,selchan); end;
    else
      if selectchan, data.blcov = data.blcov(selchan,selchan);   end;
    end
  end
  if selectrpt,
    data.numsamples = data.numsamples(selrpt);
    data.dof(:)     = length(selrpt); % FIXME this only works for equal length trials
    try, data.numcovsamples = data.numcovsamples(selrpt); end;
  end
  if selectchan, data.label = data.label(selchan); end
  if selecttime,
    data.time = data.time(seltime);
    data.dof  = data.dof(:,seltime);
    data.numsamples(:) = length(seltime);

  end

elseif issource,
  error('this is not yet implemented');

elseif isvolume,
  error('this is not yet implemented');
end
