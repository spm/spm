function D = spm_eeg_rdata_bdf(S)
% converts EEG data from BDF- to SPM-format
% FORMAT D = spm_eeg_rdata_bdf(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of bdf-file
% Fchannels   - filename of channel template
% exg_name    - cell vector that code type of exogeneous channels. Allowed
%               types are 'heog', 'veog', 'reference' and 'other' 
%_______________________________________________________________________
% 
% spm_eeg_rdata_bdf reads a continuous *.bdf file, stores everything
% in struct D and saves struct D to mat-file. The data is stored separately
% in a dat-file.
% There are calls to openbdf and readbdf, which are distributed under the
% GNU-license.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id$


try
    Fdata = S.Fdata;
catch
	Fdata = spm_get(1, '.bdf', 'Select bdf-file');
end

try
    Fchannels = S.Fchannels;
catch
    Fchannels = spm_get(1, '.mat', 'Select channel template file', fullfile(spm('dir'), 'EEGtemplates'));
end


spm('Pointer','Watch'); drawnow;

Csetup = load(Fchannels);
F = spm_str_manip(Fdata, 'rt');
P = spm_str_manip(Fdata, 'H');

D.fnamedat = [F '.dat'];

% read header
Hdata = openbdf(Fdata);

% build header for SPM format
%--------------------------------------------------
% find out how many channels
Cind = strmatch('Active Electrode', Hdata.Head.Transducer);
D.Nchannels = length(Cind);

% read channel names
for i = 1:D.Nchannels
    D.channels.name{i} = deblank(Hdata.Head.Label(i,:));
end

% go through 'external' channels and ask what they are
Cexg = find(strncmpi('exg', D.channels.name, 3));

for i = 1:length(Cexg)
    try
        D.channels.name{D.Nchannels-length(Cexg) + i} = S.exg_name{i};
    catch        
        D.channels.name{D.Nchannels-length(Cexg) + i} = spm_input(sprintf('Type of exg%d-channel', D.Nchannels-length(Cexg) + i), 'i+1','heog|veog|reference|other');
    end
end

% find heog, veog and reference channels
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

D.channels.reference = find(strncmpi('reference', D.channels.name, 9));
if isempty(D.channels.reference)
    warning(sprintf('No reference channel found.'));
    D.channels.reference = 0;
end

% Map name of channels to channel order specified in channel template file
for i = setdiff(Cind, Cexg)
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

% Assume that each channel has the same sampling rate
D.Radc = Hdata.Head.SampleRate(Cind(1));

% open dat-file for writing
%---------------------------------------------------
Fout = [F '.dat'];
D.fnamedat = Fout;
fpout = fopen(Fout, 'w');

% No scaling
D.scale.dim = 1;
D.scale.values = ones(D.Nchannels, 1);

D.events.code = [];
D.events.time = [];
D.Nsamples = 0;

% loop over seconds
% Nblocks = 500;
Nblocks = Hdata.Head.NRec;
for t = 1: Nblocks
    disp(sprintf('Second %d/%d', t, Hdata.Head.NRec))
    
    data = readbdf(Hdata,t,0);
    
    % calibrate the data
    % data.Record = data.Record.*repmat(data.Head.Cal, 1, size(data.Record, 2));

    Nsamples = size(data.Record, 2);
    D.Nsamples = D.Nsamples + Nsamples;
    
    % store data
    fwrite(fpout, data.Record(1:D.Nchannels,:), 'float');
    
    % identify events
    
    if size(data.Record, 1) >= 144
        diffs=diff(data.Record(144,:));
        idiffs=find(diffs~=0);
        idiffs=idiffs+1;
        if ~isempty(idiffs)
            for i = idiffs
                bytes=dec2bin(data.Record(144,i));
                bytes=bytes(end-7:end);
                bytes=flipdim(bytes,2);
                event=bin2dec(bytes);
                D.events.code = [D.events.code event];
                
            end
            D.events.time = [D.events.time idiffs+(t-1)*Nsamples];
        end
    end    
end
D.Nevents = 1;
D.events.types = unique(D.events.code);
D.datatype = 'float';

fclose(fpout);

D.fname = [F '.mat'];

if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer','Arrow');
