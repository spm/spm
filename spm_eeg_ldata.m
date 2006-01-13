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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_ldata.m 404 2006-01-13 18:42:21Z stefan $


try
    P;
catch
    P = spm_select(1, '\.mat$', 'Select EEG mat file');
end

Ppath = spm_str_manip(P, 'H');
if strcmp('.', Ppath) | strcmp('..', Ppath)
    Ppath = pwd;
end

try
	load(P);
catch    
	error(sprintf('Trouble reading file %s', P));
end

spm('Pointer', 'Watch');

% check whether there is a struct D
if exist('D') ~= 1
    error('%s doesn''t contain SPM M/EEG data', P);
end

if ~isfield(D, 'datatype')
    dtype = spm_type('int16');
else
    % compatablity with old spm_file_array
    if strcmp(D.datatype, 'float')
        D.datatype = 'float32';
    end
    
    dtype = spm_type(D.datatype);
end

% save path temporarily in structure
D.path = Ppath;

% convert scale factors from old spm_file_array if necessary
if isstruct(D.scale)
    % convert
    if isfield(D, 'Nfrequencies')
        dim = [D.Nchannels, D.Nfrequencies, D.Nsamples, D.Nevents];
        dim(setdiff(1:4, D.scale.dim)) = 1;
        tmp = zeros(dim);
        tmp(:,:,:,:) = D.scale.values;
        D.scale = tmp;
    else
        dim = [D.Nchannels, D.Nsamples, D.Nevents];
        dim(setdiff(1:3, D.scale.dim)) = 1;
        tmp = zeros(dim);
        tmp(:,:,:) = D.scale.values;
        D.scale = tmp;
    end
end

% memory map the data
if isfield(D, 'Nfrequencies')
    % time-frequency data
	D.data = file_array(fullfile(Ppath, D.fnamedat), [D.Nchannels D.Nfrequencies D.Nsamples D.Nevents],...
		dtype, 0, D.scale);
else
	D.data = file_array(fullfile(Ppath, D.fnamedat), [D.Nchannels D.Nsamples D.Nevents],...
	dtype, 0, D.scale);
end

spm('Pointer', 'Arrow');
