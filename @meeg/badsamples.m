function res = badsamples(this, chanind, sampind, trialind)
% Returns an array of 0/1 marking bad data based on artefact events and bad flags
% FORMAT res = badsamples(this, chanind, sampind, trialind)
% _______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: badsamples.m 6311 2015-01-21 15:44:23Z vladimir $

if ischar(chanind) && isequal(chanind, ':')
    chanind = 1:nchannels(this);
end

if ischar(sampind) && isequal(sampind, ':')
    sampind = 1:nsamples(this);
end

if ischar(trialind) && isequal(trialind, ':')
    trialind = 1:ntrials(this);
end

res = false(length(chanind), nsamples(this), length(trialind));
for i = 1:length(trialind)
    ev = events(this, trialind(i), 'samples');
    if iscell(ev)
        ev = ev{1};
    end
    
    if ~isempty(ev)
        ev = ev(intersect(strmatch('artefact', {ev.type}),...
            find(cellfun(@ischar, {ev.value}))));
        for k = 1:numel(ev)
            [dum, chan] = intersect(chanind, selectchannels(this, ev(k).value));
            res(chan, ev(k).sample+(0:(ev(k).duration-1)), i) = true;
        end
    end
end

res = res(:, sampind, :);
res(badchannels(this, chanind), :, :) = true;
res(:, :, badtrials(this, trialind))  = true;
