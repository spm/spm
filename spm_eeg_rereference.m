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
% is deleted from the data and, if possible, replaced by the old reference. If
% the new reference is an average of channels, no channels are deleted.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_rereference.m 710 2006-12-21 14:59:04Z stefan $

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

if ~isfield(D, 'reref')
	D.reref = [];
end

% type all channel names and their indices in matlab window
for i = 1:D.Nchannels
   disp(sprintf('%d: %s', i, D.channels.name{i}))
end

try
    D.reref.newref = S.newref;
catch
	D.reref.newref = ...
        spm_input('New reference channel(s)', '+1', 'n', '');
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

if ~iscell(D.channels.ref_name) 
	D.channels.ref_name = {D.channels.ref_name};
end
ref_name_old = D.channels.ref_name;
reference_old = D.channels.reference;
DNchannels_old = D.Nchannels;


% Add new reference(s)... (D.channels.reference refers to position in template file)      
D.channels.reference = D.channels.order(D.reref.newref);
D.channels.ref_name=[];

for i = 1:length(D.reref.newref)
    	D.channels.ref_name{i} = D.channels.name{D.reref.newref(i)};
end

%    if reference_old < length(D.channels.eeg)
%		D.Nchannels = D.Nchannels - 1;
%    end
%    D.channels.order(D.reref.newref) = [];
%    D.channels.name(D.reref.newref) = [];
%    D.channels.eeg = setdiff([1:D.Nchannels], [D.channels.veog D.channels.heog D.channels.Bad]);
%    
%    % re-find EOG channels
%    D.channels.heog = find(strncmpi('heog', D.channels.name, 4));
%    if isempty(D.channels.heog)
%        D.channels.veog = 0;
%    end
%    D.channels.veog = find(strncmpi('veog', D.channels.name, 4));
%    if isempty(D.channels.veog)
%        D.channels.veog = 0;
%    end
%    
%    if isfield(D.channels, 'base')
%        D.channels.base(D.reref.newref) = [];
%        D.channels.sens(D.reref.newref) = [];
%        D.channels.calib(D.reref.newref) = [];
%    end
%else
%    D.channels.reference = 0;
%    D.channels.ref_name = 'NIL';
%end

% restore old reference(s), if possible
if reference_old ~= 0
    for i = 1:length(reference_old)
        if ~ismember(reference_old(i),D.channels.order)
            D.Nchannels = D.Nchannels + 1;
            D.channels.order(D.Nchannels) = reference_old(i);
            D.channels.name{D.Nchannels} = ref_name_old{i};
            D.channels.eeg = setdiff([1:D.Nchannels], [D.channels.veog D.channels.heog D.channels.Bad]);

            %Assumes reference channels are same as others...
            if isfield(D.channels, 'base')
                D.channels.base(D.Nchannels) = D.channels.base(D.Nchannels-1);
                D.channels.sens(D.Nchannels) = D.channels.sens(D.Nchannels-1);
                D.channels.calib(D.Nchannels) = D.channels.calib(D.Nchannels-1);
                D.channels.thresholded{D.Nchannels} = [];
            end
        end
    end
end

d = zeros(D.Nchannels, D.Nsamples);
D.scale = zeros(D.Nchannels, 1, D.Nevents);

fpd = fopen(fullfile(P, D.fnamedat), 'w');

spm_progress_bar('Init', D.Nevents, 'Events re-referenced'); drawnow;
if D.Nevents > 100, Ibar = floor(linspace(1, D.Nevents,100));
else, Ibar = [1:D.Nevents]; end

for i = 1:D.Nevents
    
    d = squeeze(D.data(:, :, i));
    
    d = [d; zeros(D.Nchannels - DNchannels_old, size(d, 2))];
    
    % Subtract new reference time series from other channels
    d(D.channels.eeg, :) = d(D.channels.eeg, :) - repmat(Tref(:,i)', length(D.channels.eeg), 1);
    
    D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);
    if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
    end
end

fclose(fpd);

D.data = [];

D.fname = ['R' D.fname];

if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');


