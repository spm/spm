function res = badsamples(this, chanind, sampind, trialind)
% Returns an array of 0/1 marking bad data based on artefact events and bad flags
% FORMAT res = badsamples(this, chanind, sampind, trialind)
% _______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: badsamples.m 5610 2013-08-14 10:28:57Z vladimir $

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
        for j = 1:length(chanind)
            cev = ev(strmatch(char(chanlabels(this, chanind(j))), {ev.value}));
            for k = 1:numel(cev)
                res(j, cev(k).sample+(0:(cev(k).duration-1)), i) = true;
            end
        end
    end
end

res = res(:, sampind, :);
res(badchannels(this, chanind), :, :) = true;
res(:, :, badtrials(this, trialind))  = true;
