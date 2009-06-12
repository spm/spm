function D = spm_eeg_copy(S)
% Copy EEG/MEG data to new files
% FORMAT D = spm_eeg_copy(S)
% S           - input struct (optional)
% (optional) fields of S:
%   S.D       - MEEG object or filename of MEEG mat-file
%   S.newname - filename for the new dataset
%
% D           - MEEG object of the new dataset
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_copy.m 3199 2009-06-12 16:18:28Z guillaume $

SVNrev = '$Rev: 3199 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);

% get MEEG object
%--------------------------------------------------------------------------
try
    D   = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D       = spm_eeg_load(D);

% get filename for the new dataset
%--------------------------------------------------------------------------
if ~isfield(S, 'newname')
    S.newname = spm_input('New file name:', '+1', 's');
end

S.newname = [spm_str_manip(S.newname, 'rt') '.dat'];

% copy dataset (.mat and .dat)
%--------------------------------------------------------------------------
Dnew = clone(D, S.newname);
copyfile(D.fnamedat, Dnew.fnamedat, 'f');

D = Dnew;
D = D.history('spm_eeg_copy', S); % maybe not?
save(D);