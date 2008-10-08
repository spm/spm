function D = spm_eeg_load(P)
% read an EEG file in SPM format. 
% FORMAT D = spm_eeg_load(P)
%
% P         - filename of EEG-data file
% D         - MEEG object 
%__________________________________________________________________________
% 
% spm_eeg_load loads an MEEG file that is in SPM8 format. Importantly, the
% data is memory mapped and the struct is converted to meeg object.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_load.m 2315 2008-10-08 14:43:18Z jean $


% get filename
%--------------------------------------------------------------------------
try
    P = deblank(P);
catch
    P = spm_select(1, '\.mat$', 'Select M/EEG mat file');
end

Ppath = spm_str_manip(P, 'H');
if strcmp('.', Ppath) | strcmp('..', Ppath)
    Ppath = pwd;
end

% get filename
%--------------------------------------------------------------------------
try
    load(P);
catch    
    error(sprintf('Trouble reading file %s', P));
end

% check whether there is a struct D
%--------------------------------------------------------------------------
if exist('D') ~= 1
    error('%s doesn''t contain SPM M/EEG data', P);
end

% save path temporarily in structure
%--------------------------------------------------------------------------
D.path = Ppath;

try
    D.type;
catch
    D.type = 'other';
end

% check modality
%--------------------------------------------------------------------------
switch D.type
    case{'continuous', 'single', 'evoked','grandmean'}
        D = meeg(D);
    otherwise
end

