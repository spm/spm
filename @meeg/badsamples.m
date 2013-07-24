function res = badsamples(this, chanind, sampind, trialind)
% Returns an array of 0/1 marking bad data based on artefact events and bad flags
% FORMAT res = badsamples(this, chanind, sampind, trialind)
% _______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: badsamples.m 5592 2013-07-24 16:25:55Z vladimir $

if ischar(chanind) && isequal(chanind, ':')
    chanind = 1:nchannels(this);
end

if ischar(sampind) && isequal(sampind, ':')
    sampind = 1:nsamples(this);
end

if ischar(trialind) && isequal(trialind, ':')
    trialind = 1:ntrials(this);
end

res = false(length(chanind), length(sampind), length(trialind));
for i = 1:length(trialind)
    ev = events(this, trialind(i), 'samples');
    if iscell(ev)
        ev = ev{1};
    end
    
    if ~isempty(ev)
        ev = ev(strmatch('artefact', {ev.type}));
        for j = 1:length(chanind)
            cev = ev(strmatch(char(chanlabels(this, chanind(j))), {ev.value}));
            for k = 1:numel(cev)
                res(j, cev(k).sample+(0:(cev(k).duration-1)), i) = true;
            end
        end
    end
end

res(badchannels(this, chanind), :, :) = true;
res(:, :, badtrials(this, trialind))  = true;