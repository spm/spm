function [trl] = trialfun_ctf_epoched(cfg);

% For epoched CTF data, the trialdef structure should contain
%   cfg.trialdef.includeTrigger    = array, e.g. [3 7]
%   cfg.trialdef.excludeTrigger    = array
%   cfg.trialdef.includeConditions = cell-array, e.g. {'C1', 'C2'}
%   cfg.trialdef.excludeConditions = cell-array, e.g. {'BAD'}
%   cfg.trialdef.prestim           = 0.300         latency in seconds
%   cfg.trialdef.poststim          = 0.700         latency in seconds
% Note that definitions according to the trigger values in the stimulus
% channel can be combined with contitions defined in the CTF DataEditor.
%
% Warning: this function is deprecated, since it is file-format specific.

warning('This function is deprecated, see http://www2.ru.nl/fcdonders/fieldtrip/doku.php?id=fieldtrip:development:deprecated for more details.');

% read the header information
hdr = read_fcdc_header(cfg.headerfile);

% set empty defaults for the fields that do not exist
if ~isfield(cfg.trialdef, 'includeTrigger'), cfg.trialdef.includeTrigger = []; end
if ~isfield(cfg.trialdef, 'excludeTrigger'), cfg.trialdef.excludeTrigger = []; end
if ~isfield(cfg.trialdef, 'includeConditions'), cfg.trialdef.includeConditions = {}; end
if ~isfield(cfg.trialdef, 'excludeConditions'), cfg.trialdef.excludeConditions = {}; end
% ensure that these are cell-arrays
if ischar(cfg.trialdef.includeConditions), cfg.trialdef.includeConditions = {cfg.trialdef.includeConditions}; end
if ischar(cfg.trialdef.excludeConditions), cfg.trialdef.excludeConditions = {cfg.trialdef.excludeConditions}; end
% Read trial definitions from the ClassFile.cls in the dataset
[condNumbers,condLabels] = read_ctf_cls(fullfile(cfg.dataset, 'ClassFile.cls'));
% Read trigger list from datafile
triggerChan = read_ctf_trigger(cfg.dataset);
% Shorten trials according to pre and poststim
if isfield(cfg.trialdef,'prestim')
  nSamplesPreTmp = round(cfg.trialdef.prestim*hdr.Fs);
  nPreCorr = hdr.nSamplesPre - nSamplesPreTmp;
  if nPreCorr < 0
    error('You have defined cfg.traildef.nsamplePre too large');
  end
else
  nPreCorr = 0;
end
if isfield(cfg.trialdef,'poststim')
  nSamplesPostTmp = round(cfg.trialdef.poststim*hdr.Fs) ;
  nPostCorr = (hdr.nSamples-hdr.nSamplesPre) - nSamplesPostTmp;
  if nPostCorr < 0
    error('You have defined cfg.traildef.nsamplePost too large');
  end
else
  nPostCorr = 0;
end
% Extract the relevant trials from CTF data
trl = [];
l = 0;
for k=1:hdr.nTrials
  if approveTrial(cfg,k,condNumbers,condLabels)
    % search 5 samples around where the trigger should be
    tbeg = (k-1)*hdr.nSamples+1 + hdr.nSamplesPre-2;
    tend = (k-1)*hdr.nSamples+1 + hdr.nSamplesPre+2;
    % be carefull not to to search outside of the data
    tbeg = max(tbeg, 1);
    tend = min(tend, hdr.nSamples*hdr.nTrials);
    % Pull out the trigger chan for trial k and find the trigger value,
    trigVal = max(triggerChan(tbeg:tend));
    if approveTrigger(cfg,trigVal)
      l = l + 1;
      trl(l,1) = 1+(k-1)*hdr.nSamples + nPreCorr;
      trl(l,2) = k*hdr.nSamples  - nPostCorr;
      trl(l,3) = -hdr.nSamplesPre + nPreCorr;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for single trial definition in case of CTF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok] = approveTrigger(cfg,trigVal)
ok = 1;
if isempty(cfg.trialdef.includeTrigger)
  ok = 1;
else
  ok = 0;
  for k=1:length(cfg.trialdef.includeTrigger)
    if trigVal == cfg.trialdef.includeTrigger(k)
      ok = 1;
    end
  end
end
if ~isempty(cfg.trialdef.excludeTrigger)
  for k=1:length(cfg.trialdef.excludeTrigger)
    if trigVal == cfg.trialdef.excludeTrigger(k)
      ok = 0;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for single trial definition in case of CTF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok] = approveTrial(cfg,Trial,condNumbers,condLabels)
if isempty(condLabels)
  ok = 1;
  return
end
ok = 0;
if isfield(cfg.trialdef,'includeConditions') && ~isempty(cfg.trialdef.includeConditions)
  for k=1:length(cfg.trialdef.includeConditions)
    for l=1:length(condLabels)
      if iscell(condNumbers{l})
        dum = cell2mat(condNumbers{l});
      else
        dum = condNumbers{l};
      end
      if strcmp(char(condLabels{l}),char(cfg.trialdef.includeConditions{k})) && ismember(Trial,dum)
        fprintf('Including trial %d; Label: %s\n',Trial,char(cfg.trialdef.includeConditions{k}));
        ok = 1;
      end
    end
  end
else
  ok = 1;
end
if isfield(cfg.trialdef,'excludeConditions') && ~isempty(cfg.trialdef.excludeConditions)
  for k=1:length(cfg.trialdef.excludeConditions)
    for l=1:length(condLabels)
      if iscell(condNumbers{l})
        dum = cell2mat(condNumbers{l});
      else
        dum = condNumbers{l};
      end
      if strcmp(char(condLabels{l}),char(cfg.trialdef.excludeConditions{k})) && ismember(Trial,dum)
        fprintf('Excluding trial %d; Label: %s\n',Trial,char(cfg.trialdef.excludeConditions{k}));
        ok = 0;
      end
    end
  end
end

