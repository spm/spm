function D = spm_eeg_load(P)
% Load an M/EEG file in SPM format
% FORMAT D = spm_eeg_load(P)
%
% P        - filename of M/EEG file
% D        - MEEG object 
%__________________________________________________________________________
% 
% spm_eeg_load loads an M/EEG file using the SPM MEEG format. Importantly,
% the data array is memory-mapped and the struct is converted to MEEG object.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_load.m 5052 2012-11-14 14:31:08Z guillaume $


%-Bypass if the input is already an MEEG object
%--------------------------------------------------------------------------
if nargin && isa(P, 'meeg')
    D = P;
    return;
end

%-Get filename
%--------------------------------------------------------------------------
if ~exist(spm_file(P, 'ext', '.mat'), 'file')
    error('Cannot find file "%s".',P);
end
if ~nargin
    [P, sts] = spm_select(1, 'mat', 'Select SPM M/EEG file');
    if ~sts, D = []; return; end
end

P      = spm_file(P, 'ext', '.mat');
[p, f] = fileparts(P);
if isempty(p)
    p  = pwd;
end

%-Load MAT file
%--------------------------------------------------------------------------
try
    load(P);
catch    
    error('Trouble reading file "%s".', P);
end

%-Check whether there is a struct D
%--------------------------------------------------------------------------
if ~exist('D','var')
    error('File "%s" doesn''t contain SPM M/EEG data.', P);
end

%-Handle situations where the object has been directly saved in file
%--------------------------------------------------------------------------
if ~isa(D, 'struct')
    try
        D = struct(D);
    catch
        error('The file should contain an SPM M/EEG struct named D.');
    end
end

%-Save path and fname in structure
%--------------------------------------------------------------------------
D.path  = p;
D.fname = [f '.mat'];

%-And return an MEEG object
%--------------------------------------------------------------------------
D = meeg(D);
