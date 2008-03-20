function ind = meegchannels(obj)
% Method for getting index vector of m/eeg channels for display
% FORMAT ind = meegchannels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

type = cat(1,obj.channels(:).type);
ind = unique([find(strmatch('EEG', type)) find(strmatch('MEG', type))]);

