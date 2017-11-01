function res = badsamples(this, chanind, sampind, trialind)
% Returns an array of 0/1 marking bad data based on artefact events and bad flags
% FORMAT res = badsamples(this, chanind, sampind, trialind)
% _______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: badsamples.m 7199 2017-11-01 16:42:12Z vladimir $

if ischar(chanind) && isequal(chanind, ':')
    chanind = 1:nchannels(this);
end

if ischar(sampind) && isequal(sampind, ':')
    sampind = 1:nsamples(this);
end

if ischar(trialind) && isequal(trialind, ':')
    trialind = 1:ntrials(this);
end

if ~isequal(type(this), 'continuous') && ~any(trialonset(this))
    error('Trial onset information is not available. Cannot map artefact events to samples.');
end

res = false(length(chanind), nsamples(this), length(trialind));
for i = 1:length(trialind)
    
    ev = events(this, trialind(i));
    if iscell(ev)
        ev = ev{1};
    end
    
    if ~isempty(ev)
        ev = ev(intersect(strmatch('artefact', {ev.type}),...
            find(cellfun(@ischar, {ev.value}))));
        for k = 1:numel(ev)
            [dum, chan] = intersect(chanind, selectchannels(this, ev(k).value));
            samples     = find((trialonset(this, trialind(i))+time(this))>=ev(k).time & ...
                (trialonset(this, trialind(i))+time(this))<=(ev(k).time+ev(k).duration));
            res(chan, samples, i) = true;
        end
    end
end

res = res(:, sampind, :);
res(badchannels(this, chanind), :, :) = true;
res(:, :, badtrials(this, trialind))  = true;
