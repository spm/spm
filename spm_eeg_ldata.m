function D = spm_eeg_ldata(P)
% read an EEG file in SPM format. 
% FORMAT D = spm_eeg_ldata(P)
%
% P 		- filename of EEG-data file
% D			- EEG data struct 
%_______________________________________________________________________
% 
% spm_eeg_ldata loads an EEG file that is in SPM format. Importantly, the
% data is memory mapped and made accessible under D.data.
%_______________________________________________________________________
% Stefan Kiebel $Id$

try
    P;
catch
    P = spm_get(1, '.mat', 'Select EEG mat file');
end

Ppath = spm_str_manip(P, 'H');
if strcmp('.', Ppath)
    Ppath = pwd;
end

try
	load(P);
catch    
	error(sprintf('Trouble reading file %s', P));
end

spm('Pointer', 'Watch');
if ~isfield(D, 'datatype')
    dtype = spm_type('int16');
else
    dtype = spm_type(D.datatype);
end

% save path temporarily in structure
D.path = Ppath;

% memory map the data
if isfield(D, 'Nfrequencies')
    % time-frequency data
	D.data = spm_file_array(fullfile(Ppath, D.fnamedat), [D.Nchannels D.Nfrequencies D.Nsamples D.Nevents],...
		dtype, D.scale);
else
	D.data = spm_file_array(fullfile(Ppath, D.fnamedat), [D.Nchannels D.Nsamples D.Nevents],...
	dtype, D.scale);
end

spm('Pointer', 'Arrow');
