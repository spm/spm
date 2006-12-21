function D = spm_eeg_rdata(S)
% converts EEG data from CNT- to SPM-format
% FORMAT D = spm_eeg_rdata(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of CNT-file
% Fchannels   - filename of channel template
% reference   - name of reference channel(s)
%_______________________________________________________________________
% 
% spm_eeg_rdata reads a continuous *.cnt Neuroscan file, stores everything
% in struct D and saves struct D to mat-file. The data is stored separately
% in a dat-file.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_rdata.m 710 2006-12-21 14:59:04Z stefan $

try
    Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.cnt$', 'Select cont. neuroscan-file');
end

try
    Fchannels = S.Fchannels;
catch
    Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
end

% There doesn't seem to be information in the CNT-file which channel(s) is the reference
% so ask for it now
try
    S.reference;
catch
    S.reference = spm_input('Input reference channel name', '+1', 's');
end

spm('Pointer','Watch'); drawnow;

D = struct('data', [], 'channels', [], 'scale', [], 'filter', [], 'events', [], 'reref', [], 'descrip', []);
D.channels = struct('name', [], 'base', [], 'sens', [], 'calib', [], 'order', []);
D.events = struct('time', [], 'code', []);

Csetup = load(Fchannels);

fp = fopen(Fdata, 'r');
F = spm_str_manip(Fdata, 'rt');
P = spm_str_manip(Fdata, 'H');

D.fnamedat = [F '.dat'];
fpd = fopen(fullfile(P, D.fnamedat), 'w');

[fullname, permission, mformat] = fopen(fp);

% file size
Dn = dir(fullname);
Sf = Dn.bytes;

% Read data and time strings
s = fseek(fp, 225, -1);
D.descrip.date = char(fread(fp, 10, 'uchar')');

s = fseek(fp, 235, -1);
D.descrip.time = char(fread(fp, 12, 'uchar')');

% Read number of channels
s = fseek(fp, 370, -1);
Nchannels = fread(fp, 1, 'ushort');

% Read AD rate
s = fseek(fp, 376, -1);
D.Radc = fread(fp, 1, 'ushort');

% Read the event table offset
s = fseek(fp, 886, -1);
Oevent = fread(fp, 1, 'int');

% Read the events
%-----------------------------------------------
s = fseek(fp, Oevent, -1);

% Type of event
Tevents = fread(fp, 1, 'uchar');

% Number of events
Nevents = fread(fp, 1, 'int');

% read events from cnt-file
if Tevents == 1
	Nevents = Nevents/8;
else 
	Nevents = Nevents/19;
end

Et = zeros(1, Nevents);
Ec = zeros(1, Nevents);
for i = 1:Nevents
	if Tevents == 1
		fseek(fp, Oevent + 9 + 8*(i-1), -1);
	else
		fseek(fp, Oevent + 9 + 19*(i-1), -1);
	end
	
	Ec(i) = fread(fp, 1, 'uchar');

	if Tevents == 1
		fseek(fp, Oevent + 9 + 8*(i-1) + 4, -1);
	else
		fseek(fp, Oevent + 9 + 19*(i-1) + 4, -1);
	end
	
	offset = fread(fp, 1, 'int');
	
	Et(i) = (offset - 900 - 75*Nchannels)/(2*Nchannels);
end

% Number of expected samples
Nsamples = (Oevent - 900 - 75*Nchannels)/(2*Nchannels); 

% store only relative name of channel template file
D.channels.ctf = spm_str_manip(Fchannels, 't');

%%%%%%%% RIK
% If necessary, insert here option to read in new gain file
%
% Eg: [cn,cg]=textread('/home/rhenson/matlab/SPMEEG/cal.txt','%s%f');
%
% match cn to D.channels.name, and then update D.channels.calib to cg
%%%%%%%%%

% Read the electrode names
for i = 1:Nchannels
	fseek(fp, 900 + 75*(i-1), -1);
	tmp = char(fread(fp, 10, 'uchar'))';
	
	D.channels.name{i} = deblank(tmp);
	
end

% Read HEOG channel number (but doesn't seem to be stored in header)
% s = fseek(fp, 406, -1);
% D.channels.heog = fread(fp, 1, 'short');
%
% % Read VEOG channel number (but doesn't seem to be stored in header)
% s = fseek(fp, 408, -1);
% D.channels.veog = fread(fp, 1, 'short');

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

% Map name of channels to channel order specified in channel template file
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

% Report about reference channel
index = [];
for j = 1:Csetup.Nchannels
    if ~isempty(find(strcmpi(S.reference, Csetup.Cnames{j})))
        index = j;
    end
end
if isempty(index)
    warning(sprintf('Reference %s not found in channel template file.', S.reference));
    D.channels.reference = 0;
else
    D.channels.reference = index;
end
D.channels.ref_name = S.reference;

% initialize Bad channels vector
D.channels.Bad = [];

% indices of EEG channels (without EOGs)
D.channels.eeg = setdiff(1:length(D.channels.order), [D.channels.heog D.channels.veog]);

% Read baseline, sensitivity and calibration to convert to microvolts.
for i = 1:Nchannels	
	fseek(fp, 900 + 75*(i-1) + 47, -1);
	D.channels.base(i) = fread(fp, 1, 'short')';

	fseek(fp, 900 + 75*(i-1) + 59, -1);
	D.channels.sens(i) = fread(fp, 1, 'float')';
		
	fseek(fp, 900 + 75*(i-1) + 71, -1);
	D.channels.calib(i) = fread(fp, 1, 'float')';
end

% Read data
%-----------------------------------------------
fseek(fp, 900 + 75*Nchannels, -1);

D.scale = zeros(Nchannels, 1);

% if any(D.channels.base~=0)
% 	error('base is unequal zero');
% end

% Conversion to microvolt saved as scalefactor
D.scale = (1/204.8 .* D.channels.sens .* D.channels.calib)';

% progress bar
spm_progress_bar('Init', Nsamples, 'Time bins converted'); drawnow;
if Nsamples > 100, Ibar = floor(linspace(1, Nsamples,100));
else, Ibar = [1:Nsamples]; end

for i = 1:Nsamples
	d = int16(fread(fp, Nchannels, 'short'));
	fwrite(fpd, d, 'int16');

    if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
    end

end

fclose(fp);
fclose(fpd);

D.Nchannels = Nchannels;
D.Nsamples = Nsamples;
D.Nevents = 1; % continuous data
D.events.time = Et;
D.events.code = Ec;
D.fname = [F '.mat'];

D.datatype = 'int16';
D.units = '\muV';

D.modality = 'EEG';

if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm_progress_bar('Clear');

spm('Pointer','Arrow');
