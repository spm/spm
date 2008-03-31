function ind = meegchannels(this)
% Method for getting index vector of m/eeg channels for display
% FORMAT ind = meegchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

type = chantype(this);
ind = unique([strmatch('EEG', type); strmatch('MEG', type)]);

