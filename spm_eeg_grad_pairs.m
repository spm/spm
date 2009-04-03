function pairs = spm_eeg_grad_pairs(D)
% takes an MEEG object and returns a two-column matrix of pairs of
% gradiometers (for neuromag FIF labelling)
%______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Laurence Hunt
% $Id: spm_eeg_grad_pairs.m 3047 2009-04-03 08:28:59Z vladimir $

megplanars = D.meegchannels('MEGPLANAR');
channames = D.chanlabels;

pairs = [];
for i = 1:length(megplanars)
    if sum(strncmpi(channames(megplanars),channames{megplanars(i)},6))==2 %pair of planars exists
        pairs = [pairs; megplanars(find(strncmpi(channames(megplanars),channames{megplanars(i)},6)))];
        channames{megplanars(i)}='Ive been found'; %to stop us from duplicating this pair in list
    end
end

if size(pairs,2)~=2
    error('Failed to find pairs of planar gradiometers. File may have non-standard chanlabels.')
end

