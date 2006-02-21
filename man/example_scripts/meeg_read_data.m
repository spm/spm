% This is an example script to convert epoched and averaged data (ERP/ERF) 
% to SPM format. Please adapt to your needs.

clear all

%-----------------------------------
% part 1: information about the data
%-----------------------------------
% change the following values according to your data

% sampling rate in Hz
D.Radc = 200;

% Number of time bins in peri-stimulus time
D.Nsamples = 101;

% Number of ERPs
D.Nevents = 2;

% channel template file
D.channels.ctf = 'montage1.mat';

% Names of channels in order of the data
D.channels.name = {'O1', 'O2', 'Cz'};
D.Nchannels = length(D.channels.name);

% index vector of bad channels
D.channels.Bad = [];

% nr of time bins before stimulus onset (this excludes the time bin at
% zero)
D.events.start = 20;

% nr of time bins after stimulus onset (this excludes the time bin at
% zero)
D.events.stop = 80;

% name of SPM mat file
D.fname = 'test.mat';

D.modality = 'EEG';
D.units = '\muV';
%-------------------------------------
% Don't change these
D.channels.reference = 0;
D.channels.ref_name = 'none';

D.events.types = [1:D.Nevents];
D.events.Ntypes = length(D.events.types);
D.events.code = [1:D.Nevents];
D.events.reject = zeros(1, D.events.Ntypes);

D.path = fileparts(D.fname);
D.fname = spm_str_manip(D.fname, 't');
D.fnamedat = [spm_str_manip(D.fname, 'r') '.dat'];

%----------------------------------------------------
% part 2: read data: adapt this to read your own data
%----------------------------------------------------

data = randn(D.Nchannels, D.Nsamples, D.Nevents);

%-----------------------------------
% part 3: write data, don't change
%-----------------------------------
Csetup = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));

% Map name of channels (EEG only) to channel order specified in channel
% template file. 
D.channels.heog = find(strncmpi('heog', D.channels.name, 4));
if isempty(D.channels.heog)
    warning(sprintf('No HEOG channel found.'));
    D.channels.heog = 0;
else
    D.channels.heog = D.channels.heog(1);
end

D.channels.veog = find(strncmpi('veog', D.channels.name, 4));
if isempty(D.channels.veog)
    warning(sprintf('No VEOG channel found.'));
    D.channels.veog = 0;
else
    D.channels.veog = D.channels.veog(1);
end

D.channels.eeg = setdiff([1:D.Nchannels], [D.channels.veog D.channels.heog]);

for i = D.channels.eeg
	index = [];
	for j = 1:Csetup.Nchannels
		if ~isempty(find(strcmpi(D.channels.name{i}, Csetup.Cnames{j})))
			index = [index j];
		end
	end
    
	if isempty(index)
		warning(sprintf('No channel named %s found in channel template file.', D.channels.name{i}));
	else
		% take only the first found channel descriptor
		D.channels.order(i) = index(1);
	end

end

D.scale = ones(D.Nchannels, 1, D.Nevents);
D.datatype  = 'float32';

fpd = fopen(fullfile(D.path, D.fnamedat), 'w');

for i = 1:D.Nevents
    d = squeeze(data(:, :, i));
    D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);        						
end


fclose(fpd);

if str2num(version('-release'))>=14 
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
    save(fullfile(D.path, D.fname), 'D');
end

spm('Pointer', 'Arrow');

