function D = spm_eeg_rdata_dat(S)
% converts EEG data from neuroscan ASCII format to SPM-format
% FORMAT D = spm_eeg_rdata_dat(S)
%
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of ASCII dat-file
% events      - struct with the following (optional) fields
%    start    - number of samples before stimulus onset
%    stop     - number of samples after stimulus onset
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_rdata_dat.m 539 2006-05-19 17:59:30Z Darren $


try
    Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.dat$', 'Select dat neuroscan-file');
end

% Read dat-file into cell vector of strings
txt = textread(Fdata, '%s', 'delimiter', '\n', 'whitespace', '');

% Fill in the values in struct D
D.channels = [];
D.events = [];
ind = strmatch('[Channels]', txt);
D.Nchannels = sscanf(txt{ind(1)}, '[Channels] %d');

D.channels.ctf = 'montage1.mat';

ind = strmatch('[Rate]', txt);
D.Radc = sscanf(txt{ind(1)}, '[Rate] %f');

try
    D.events.start = S.start;
catch
    D.events.start =...
        spm_input('Nr of samples before stimulus onset', '+1', 'i', '', 1);
end

try
    D.events.stop = S.stop
catch
    D.events.stop = ...
        spm_input('Nr of samples after stimulus onset', '+1', 'i', '', 1);
end

ind = strmatch('[Points]', txt);
D.Nsamples = sscanf(txt{ind(1)}, '[Points] %f');

ind = strmatch('[Sweeps]', txt);
D.Nevents = sscanf(txt{ind(1)}, '[Sweeps] %f');

ind = strmatch('[Electrode Labels]', txt);

ind2 = strfind(txt{ind(1)+1}, ']');
if length(ind2) ~= D.Nchannels, warning('Channel mismatch'), end

for i = 1 : D.Nchannels
    tmp = txt{ind(1)+1}(ind2(i)-8:ind2(i)-1);
    tmp = deblank(tmp(end:-1:1));
    D.channels.name{i} = tmp(end:-1:1);
end

D.channels.heog = find(strncmpi('heog', D.channels.name, 4));
if isempty(D.channels.heog)
    warning(sprintf('No HEOG channel found.'));
else
    D.channels.heog = D.channels.heog(1);
end

D.channels.veog = find(strncmpi('veog', D.channels.name, 4));
if isempty(D.channels.veog)
    warning(sprintf('No VEOG channel found.'));
else
    D.channels.veog = D.channels.veog(1);
end

Csetup = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));
% Map name of channels to channel order specified in CTF
for i = 1:length(D.channels.name)
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

% read data
ind_type = strmatch('[Trial Type]', txt);
ind_accept = strmatch('[Accept]', txt);
ind_data = strmatch('[Epoch Data]', txt);

if length(ind_type) ~= D.Nevents, warning('Data mismatch'), end
if length(ind_accept) ~= D.Nevents, warning('Data mismatch'), end
if length(ind_data) ~= D.Nevents, warning('Data mismatch'), end

% Assemble epochs
k = 1;

D.scale = zeros(D.Nchannels, D.Nevents);

D.events.code = zeros(1, D.Nevents);
D.events.reject = zeros(1, D.Nevents);

% open outputfile, take care not to overwrite input file, add '_spm'
P = spm_str_manip(Fdata, 'H');

Fout = [spm_str_manip(Fdata, 'rt') '_spm.dat'];
D.fnamedat = Fout;
Fmat = [spm_str_manip(Fdata, 'rt') '_spm.mat'];
D.fname = Fmat;

fpout = fopen(fullfile(P, Fout), 'w');
% D.datatype = 'int16';
D.datatype = 'float';

for i = 1 : D.Nevents
    D.events.reject(i) = 1 - sscanf(txt{ind_accept(i)}, '[Accept] %d');
    D.events.code(i) = sscanf(txt{ind_type(i)}, '[Trial Type] %d');
    
    d = sscanf(cat(2, txt{ind_data(i)+1 : ind_data(i) + D.Nsamples}), '%f', [D.Nchannels D.Nsamples]);
    
    D.scale(:, 1, k) = spm_eeg_write(fpout, d, 2, D.datatype);
       
    k = k +1;
    
end

D.events.types = unique(D.events.code);
D.events.Ntypes = length(D.events.types);

fclose(fpout);

if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

