function D = spm_eeg_rereference(S);
% rereference EEG data to new reference channel(s)
% FORMAT D = spm_eeg_rereference(S)
%
% S		    - struct (optional)
% (optional) fields of S:
% D			- filename of EEG-file with continuous data
% Nnewref   - number of new reference channels
% newref    - new reference channels
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
% 
% spm_eeg_rereference references all data to a new channel or average of
% channels. If there is only one reference channel, the reference channel
% is deleted from the data and saved to a special reference channel. If
% the new reference is an average of channels, no channels are deleted. The
% new reference(s) are indicated by indices (in their order in the file).
% This order can be best looked at in D.channels.name.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id$

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG rereference setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

D.fnamedat = ['R' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');

if ~isfield(D, 'reref')
	D.reref = [];
end

try
    D.reref.newref = S.newref;
catch
	D.reref.newref = ...			
        spm_input('Indices of new reference channel(s)', '+1', 'n', '');
end

if ~all(ismember(D.reref.newref, D.channels.eeg))
    error(['New references can only have the following indices: ' sprintf('%d ', D.channels.eeg)]);
end

% get time course of new reference for all events
Tref = zeros(D.Nsamples, D.Nevents);
for i = 1:D.Nevents
	if length(D.reref.newref) > 1
		Tref(:,i) = mean(squeeze(D.data(D.reref.newref, :, i)))';
	else
		Tref(:,i) = squeeze(D.data(D.reref.newref, :, i))';
	end
end

Csetup = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));

% assume that the old reference channel doesn't need to be restored

% save the data
fpd = fopen(fullfile(P, D.fnamedat), 'w');

% delete new reference channel from channel struct
if length(D.reref.newref) == 1
    D.Nchannels = D.Nchannels - 1;
    D.channels.eeg = setdiff(D.channels.eeg, D.reref.newref);
    D.channels.order(D.reref.newref) = [];
    D.channels.name(D.reref.newref) = [];
    D.channels.base(D.reref.newref) = [];
    D.channels.sens(D.reref.newref) = [];
    D.channels.calib(D.reref.newref) = [];

end

d = zeros(D.Nchannels, D.Nsamples);
D.scale.dim = [1 3];
D.scale.values = zeros(D.Nchannels, D.Nevents);

% if there wasn't a reference channel before:
if ~isfield(D.channels, 'reference')
    D.channels.reference = length(D.channels.name)+1;
    D.channels.name{D.channels.reference} = 'reference';
    D.scale.values = [D.scale.values; zeros(1, D.Nevents)];
    D.Nchannels = D.Nchannels +1;
end

for i = 1:D.Nevents
    
    % concatenate data for each event
    d = squeeze(D.data(:, :, i));
    
    % Subtract new reference time series from other channels
    d(D.channels.eeg, :) = d(D.channels.eeg, :) - repmat(Tref(:,i)', length(D.channels.eeg), 1);
    
    % remove data of new (single) reference channel
    if length(D.reref.newref) == 1
        d(D.channels.eeg, :) = [];
    end

    % save reference time courses (and overwrite old reference)

    d(D.channels.reference, :) = Tref(:,i)';
    D.scale.values(:, i) = spm_eeg_write(fpd, d, 2, D.datatype);
							
end
fclose(fpd);

D.data = [];

D.fname = ['R' D.fname];

if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');


