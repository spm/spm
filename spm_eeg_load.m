function D = spm_eeg_load(P)
% Load an M/EEG file in SPM format
% FORMAT D = spm_eeg_load(P)
%
% P         - filename of M/EEG file
% D         - MEEG object 
%__________________________________________________________________________
% 
% spm_eeg_load loads an M/EEG file using the SPM MEEG format. Importantly, the
% data array is memory-mapped and the struct is converted to MEEG object.
%__________________________________________________________________________
% Copyright (C) 2008-2011 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_load.m 5025 2012-10-31 14:44:13Z vladimir $

% bypass if the input is already an MEEG object
%--------------------------------------------------------------------------
if nargin && isa(P, 'meeg')
    D = P;
    return;
end

% get filename
%--------------------------------------------------------------------------

if nargin ==0 || ~exist(P, 'file')    
    [P, sts] = spm_select(1, 'mat', 'Select SPM M/EEG file');
    if ~sts, D = []; return; end
end

P = deblank(P);

[p, f] = fileparts(P);
if isempty(p)
    p = pwd;
end

% load MAT file
%--------------------------------------------------------------------------
try
    load(P);
catch    
    error('Trouble reading file %s', P);
end

% check whether there is a struct D
%--------------------------------------------------------------------------
if ~exist('D','var')
    error('%s doesn''t contain SPM M/EEG data', P);
end

% This is for the case when people save the object in a file
%--------------------------------------------------------------------------
if ~isa(D, 'struct')
    try
        D = struct(D);
    catch
        error('The file should contain an SPM M/EEG struct named D');
    end
end

% save path and fname in structure
%--------------------------------------------------------------------------
D.path = p;
D.fname = [f '.mat'];

% return an MEEG object
%--------------------------------------------------------------------------
D = meeg(D);
