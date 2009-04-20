function [data] = selectdata(varargin)

% this function serves to concatenate the input data-structures along the
% compatible dimensions and thus is a more general implementation of
% appenddata, which only raw data. Moreover, it can be used to equate the
% data of different conditions to match e.g. in channels time-axis etc
% Moreover, it can be used as a generalization to ...average with 'keepindividual'
%
% Finally, this function serves to subselect regions-of-interest from the input data,
% either or not averaging across the specified dimensions.
%
% Supported input data:
%   freq
%   timelock
%   source   (not yet)
%   volume   (not yet)
%
% supported options:
%   foilim
%   toilim
%   roi
%   channel (FIXME this is also done by preprocessing?)
%   avgoverchan
%   avgoverfreq
%   avgovertime
%   avgoverroi
%   avgoverrpt

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% $Log: selectdata.m,v $
% Revision 1.6  2009/04/14 18:29:32  roboos
% deleted the subfunction istrue, since it now is a seperate function
%
% Revision 1.5  2009/03/18 19:49:54  roboos
% use the smart xxxdim functions
%
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

data   = varargin(isdata);
kvp    = varargin(keyvals);
dtype  = cell(1,length(data));
dimord = cell(1,length(data));

for k = 1:length(data)
  data{k} = checkdata(data{k}, 'datatype', {'freq' 'timelock' 'source', 'volume'});
  [dtype{k}, dimord{k}]  = datatype(data{k});
end

if any(~strmatch(dtype{1},dtype))
  error('different types of input data is not supported');
end

% check consistency of input data
if any(~strmatch(dimord{1},dimord))
  error('a different dimord in the input data is not supported');
end

isfreq   = datatype(data{1},'freq');
istlck   = datatype(data{1},'timelock');
issource = datatype(data{1},'source');
isvolume = datatype(data{1},'volume');

selchan  = keyval('channel', kvp); selectchan = ~isempty(selchan);
selfoi   = keyval('foilim',  kvp); selectfoi  = ~isempty(selfoi);
seltoi   = keyval('toilim',  kvp); selecttoi  = ~isempty(seltoi);
selroi   = keyval('roi',     kvp); selectroi  = ~isempty(selroi);
selrpt   = keyval('rpt',     kvp); selectrpt  = ~isempty(selrpt);
param    = keyval('param',   kvp); if isempty(param), param = 'all'; end

avgoverchan  = keyval('avgoverchan',  kvp); if isempty(avgoverchan), avgoverchan = false; end
avgoverfreq  = keyval('avgoverfreq',  kvp); if isempty(avgoverfreq), avgoverfreq = false; end
avgovertime  = keyval('avgovertime',  kvp); if isempty(avgovertime), avgovertime = false; end
avgoverroi   = keyval('avgoverroi',   kvp); if isempty(avgoverroi),  avgoverroi  = false; end
avgoverrpt   = keyval('avgoverrpt',   kvp); if isempty(avgoverrpt),  avgoverrpt  = false; end

% ensure that these are boolean arguments, optionally convert from "yes"/"no" to true/false
avgoverchan = istrue(avgoverchan);
avgoverfreq = istrue(avgoverfreq);
avgovertime = istrue(avgovertime);
avgoverroi  = istrue(avgoverroi);
avgoverrpt  = istrue(avgoverrpt);

if length(data)>1 && selectrpt,
  error('multiple data structures as input is not supported in combination with subselection of trials');
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
  dimord = dimord{1};

else
  % nothing to do
  data = data{1};
  dimord = dimord{1};
end

% determine the subselection in the data
if selectrpt,
  dimtok = tokenize(data.dimord, '_');
  if strcmp(dimtok{1}, 'rpttap'),
    error('here you have to ensure the correct handling of tapers');
  else
    % do nothing
  end
end

if selectchan,
  selchan = match_str(data.label, channelselection(selchan, data.label));
end

if selectfoi,
  if length(selfoi)==1, selfoi(2) = selfoi; end;
  selfoi = nearest(data.freq, selfoi(1)):nearest(data.freq, selfoi(2));
end

if selecttoi,
  if length(seltoi)==1, seltoi(2) = seltoi; end;
  seltoi = nearest(data.time, seltoi(1)):nearest(data.time, seltoi(2));
end

if selectroi,
  error('not yet implemented');
end

if isfreq,
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt);  end
  if selectchan, data = seloverdim(data, 'chan', selchan); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi);  end
  if selecttoi,  data = seloverdim(data, 'time', seltoi);  end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt');   end
  if avgoverchan, data = avgoverdim(data, 'chan');  end
  if avgoverfreq, data = avgoverdim(data, 'freq');  end
  if avgovertime, data = avgoverdim(data, 'time');  end

elseif istlck,
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt);  end
  if selectchan, data = seloverdim(data, 'chan', selchan); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi);  end
  if selecttoi,  data = seloverdim(data, 'time', seltoi);  end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt');   end
  if avgoverchan, data = avgoverdim(data, 'chan');  end
  if avgoverfreq, data = avgoverdim(data, 'freq');  end
  if avgovertime, data = avgoverdim(data, 'time');  end

elseif issource,
  error('this is not yet implemented');

elseif isvolume,
  error('this is not yet implemented');
end

