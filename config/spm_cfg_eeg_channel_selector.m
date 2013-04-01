function channels = spm_cfg_eeg_channel_selector(jobtree)
% generic M/EEG channel selector based on label and type
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_channel_selector.m 5375 2013-04-01 17:12:53Z vladimir $

if nargin == 0 || ischar(jobtree)
    chanall = cfg_const;
    chanall.tag = 'all';
    chanall.name = 'All';
    chanall.val = {'all'};
    
    type = cfg_menu;
    type.tag = 'type';
    type.name = 'Select channels by type';
    type.help = {'Select channels by type'};
    type.labels = {'MEG', 'MEGPLANAR', 'MEGMAG', 'MEGGRAD', 'MEGCOMB','EEG', 'EOG', 'ECG', 'EMG', 'LFP', 'PHYS', 'Other', 'REF', 'REFMAG', 'REFGRAD'};
    type.values = {'MEG', 'MEGPLANAR', 'MEGMAG', 'MEGGRAD', 'MEGCOMB','EEG', 'EOG', 'ECG', 'EMG', 'LFP', 'PHYS', 'Other', 'REF', 'REFMAG', 'REFGRAD'};
    
    chan = cfg_entry;
    chan.tag = 'chan';
    chan.name = 'Custom channel';
    chan.strtype = 's';
    chan.num = [1 Inf];
    chan.help = {'Enter a single channel name.'};
    
    chanfile = cfg_files;
    chanfile.tag = 'chanfile';
    chanfile.name = 'Channel file';
    chanfile.filter = 'mat';
    chanfile.num = [1 1];
    
    channels = cfg_repeat;
    channels.tag = 'channels';
    channels.name = 'Channel selection';
    % Sometimes it doesn't make sense to select by type
    if nargin == 0
        channels.values = {chanall, type, chan, chanfile};
    else
        channels.values = {chanall, chan, chanfile};
    end
    channels.num = [1 Inf];
    channels.val = {chanall};
else
    channels = {};
    for j = 1:numel(jobtree)
        if isfield(jobtree{j}, 'type')
            channels = [channels {jobtree{j}.type}];
        elseif isfield(jobtree{j}, 'all')
            channels = [channels {'all'}];
        elseif isfield(jobtree{j}, 'chan')
            channels = [channels {jobtree{j}.chan}];
        elseif isfield(jobtree{j}, 'chanfile')
            channels = [channels getfield(load(char(jobtree{j}.chanfile)), 'label')];
        end
    end
end