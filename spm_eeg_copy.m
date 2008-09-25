function newD = spm_eeg_copy(S)
% function used for epoching continuous EEG/MEG data
% FORMAT D = spm_eeg_epochs(S)
%
% S  - filename or input struct (optional)
% (optional) fields of S:
% S.D         - MEEG object or filename of MEEG mat-file
% S.newname   - name for the new dataset
% 
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_copy.m 2198 2008-09-25 17:54:45Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','MEEG copy',0);

if nargin == 0
    S =[];
end

if ~isfield(S, 'D')
    S.D = spm_select(1, 'mat', 'Select M/EEG mat file');
end

D = S.D;

if ~isa(D, 'meeg')
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

if ~isfield(S, 'newname')
    S.newname = spm_input('New file name:', '+1', 's');
end

S.newname = [spm_str_manip(S.newname, 'rt') '.dat'];

newD = clone(D, S.newname);

save(newD);

copyfile(D.fnamedat, newD.fnamedat, 'f');