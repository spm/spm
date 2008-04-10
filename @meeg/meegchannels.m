function ind = meegchannels(this)
% Method for getting index vector of m/eeg channels for display
% FORMAT ind = meegchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

type = chantype(this);
ind = union(find(strcmpi('EEG', type)), find(strcmpi('MEG', type)));
ind = ind(:)'; % must be row to allow to use it as loop indices

