function [event] = read_trigger(filename, varargin)

% READ_TRIGGER extracts the events from a continuous trigger channel
% This function is a helper function to read_event and can be used for all
% dataformats that have one or multiple continuously sampled TTL channels
% in the data.
%
% This is a helper function for READ_EVENT
%
% TODO
%  - merge read_ctf_trigger into this function (requires trigshift and bitmasking option)
%  - merge biosemi code into this function (requires bitmasking option)

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: read_trigger.m,v $
% Revision 1.2  2008/04/29 14:54:39  roboos
% explicit specification of begsample and endsample, otherwise event.sample remains empty
%
% Revision 1.1  2008/04/29 13:53:50  roboos
% new implementation, works for ctf, bti and neuromag
%

event = [];

% get the optional input arguments
hdr         = keyval('header',      varargin);
begsample   = keyval('begsample',   varargin);
endsample   = keyval('endsample',   varargin);
chanindx    = keyval('chanindx',    varargin);
detectflank = keyval('detectflank', varargin);
fixctf      = keyval('fixctf',      varargin); if isempty(fixctf),      fixctf = 0;       end
fixneuromag = keyval('fixneuromag', varargin); if isempty(fixneuromag), fixneuromag = 0;  end

if isempty(begsample)
  begsample = 1;
end

if isempty(endsample)
  endsample = hdr.nSamples*hdr.nTrials;
end

% read the trigger channel as raw data, can safely assume that it is continuous
dat = read_data(filename, 'header', hdr, 'begsample', begsample, 'endsample', endsample, 'chainindx', chanindx);

if fixctf
  % correct for reading the data as signed 32-bit integer, whereas it should be interpreted as an unsigned int
  dat(dat<0) = dat(dat<0) + 2^32;
end

if fixneuromag
  % according to Joachim Gross, real events always have triggers > 5
  % this is probably to avoid the noisefloor
  dat(dat<6) = 0;
end

for i=1:length(chanindx)
  % process each trigger channel independently
  channel = hdr.label{chanindx(i)};
  trig    = dat(i,:);
  switch detectflank
    case 'up'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([0 trig(:)'])>0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j);                % assign the trigger value just _after_ going up
      end
    case 'down'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([0 trig(:)'])<0)
        event(end+1).type   = channel;
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j-1);              % assign the trigger value just _before_ going down
      end
    case 'both'
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([0 trig(:)'])>0)
        event(end+1).type   = [channel '_up'];   % distinguish between up and down flank
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j);                % assign the trigger value just _after_ going up
      end
      % convert the trigger into an event with a value at a specific sample
      for j=find(diff([0 trig(:)'])<0)
        event(end+1).type   = [channel '_down']; % distinguish between up and down flank
        event(end  ).sample = j + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(j-1);              % assign the trigger value just _before_ going down
      end
    otherwise
      error('incorrect specification of ''detectflank''');
  end
end
