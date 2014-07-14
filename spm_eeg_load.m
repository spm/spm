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
% $Id: spm_eeg_load.m 6102 2014-07-14 09:19:09Z vladimir $


%-Bypass if the input is already an MEEG object
%--------------------------------------------------------------------------
if nargin && isa(P, 'meeg')
    D = P;
    return;
end

%-Get filename
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select(1, 'mat', 'Select SPM M/EEG file');
    if ~sts, D = []; return; end
end

P = spm_file(spm_file(P, 'ext', '.mat'), 'cpath');

%-Load MAT file
%--------------------------------------------------------------------------
try
    load(P);
catch
    
    try
        load(spm_file(P,'filename'))
        
        warning('Ignoring path. The dataset was loaded from the current directory.');
    catch
        error('Cannot load file "%s".', P);
    end
    
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
D.path  = spm_file(P, 'path');
D.fname = spm_file(P, 'filename');

%-And return an MEEG object
%--------------------------------------------------------------------------
D = meeg(D);
