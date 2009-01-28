function trl = trialfun_realtime(cfg)
% TRIALFUN_REALTIME can be used in conjunction with rt_process to read realtime
% data. Trials are defined as [begsample endsample offset condition]
%
% options:
% cfg.minsample = the last sample number that was already considered (passed from rt_process)
% cfg.blocksize = [offset length] in seconds. In case of events, offset is
%                 wrt the trigger. In case of no events, offset is wrt
%                 prevSample. E.g., [-0.9 1] will read 1 second blocks with
%                 0.9 second overlap (default = [0 1]).
% cfg.bufferdata = {'first' 'last'}. If 'last' then only the last block of
%                 interest is read. Otherwise, all well-defined blocks are read (default = 'first')
%
% Copyright (C) 2009, Marcel van Gerven
%
% $Log: trialfun_realtime.m,v $
% Revision 1.2  2009/01/28 10:56:39  marvger
% we do not consider samples before minsample+1 in asynchronous mode
%
% Revision 1.1  2009/01/28 09:52:01  marvger
% generic trialfun for realtime processing
%

  if ~isfield(cfg,'minsample'),   cfg.minsample = 0;        end
  if ~isfield(cfg,'blocksize'),   cfg.blocksize = [0 1];    end
  if ~isfield(cfg,'bufferdata'),  cfg.bufferdata = 'first'; end
  if ~isfield(cfg,'triggers'),    cfg.triggers = [];        end
  
  % blocksize in terms of samples
  cfg.blocksize = round(cfg.blocksize * cfg.hdr.Fs);
      
  % retrieve trials of interest
  if isempty(cfg.event) % asynchronous mode
    trl = trialfun_asynchronous(cfg);
  else % synchronous mode
    trl = trialfun_synchronous(cfg);
  end
end

function trl = trialfun_asynchronous(cfg)
  
  prevSample = cfg.minsample;

  trl = [];

  if strcmp(cfg.bufferdata, 'last') % only get last block

    begsample  = cfg.hdr.nSamples*cfg.hdr.nTrials - cfg.blocksize(2) + 1;

    if begsample >= (prevSample + cfg.blocksize(1))
      endsample  = cfg.hdr.nSamples*cfg.hdr.nTrials;
      trl = [trl; [begsample endsample 0 nan]];
    end

  else % get all blocks

    while true

      % see whether new samples are available
      newsamples = (cfg.hdr.nSamples*cfg.hdr.nTrials-prevSample);

      if newsamples>=sum(cfg.blocksize)

        % we do not consider samples that have already been processed
        begsample  = max(cfg.minsample+1,prevSample+cfg.blocksize(1));
        endsample  = max(cfg.minsample+1,prevSample+sum(cfg.blocksize));        
        
        trl = [trl; [begsample endsample 0 nan]];
        prevSample = endsample;

      else
        break;
      end
    end
  end
end

function trl = trialfun_synchronous(cfg)

  trl = [];
  offset = cfg.hdr.Fs * cfg.blocksize(1);

  % process all events
  for j=1:length(cfg.event)

    if isempty(cfg.triggers)
      curtrig = cfg.event(j).value;
    else
      [m1,curtrig] = ismember(cfg.event(j).value,cfg.triggers);
    end
    
    if isempty(curtrig), curtrig = nan; end

    if isempty(cfg.triggers) || (~isempty(m1) && m1)
      % catched a trigger of interest

      begsample = cfg.event(j).sample + cfg.blocksize(1);
      endsample = begsample + cfg.blocksize(2);

      trl = [trl; [begsample endsample begsample + offset curtrig]];

    end
  end
end
