function D = spm_eeg_copy(S)
% Copy EEG/MEG data to new files
% FORMAT D = spm_eeg_copy(S)
% S           - input struct (optional)
%  fields of S:
%   S.D       - MEEG object or filename of MEEG mat-file
%   S.outfile - filename for the new dataset
%   S.updatehistory - update history information [default: true]
%
% D           - MEEG object of the new dataset
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename);
spm('FigName','M/EEG copy');

if ~isfield(S, 'updatehistory'), S.updatehistory = 1; end

D = spm_eeg_load(S.D);


D = copy(D, S.outfile);

if ~isfield(S, 'updatehistory') || S.updatehistory
    D = D.history('spm_eeg_copy', S); 
end

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG copy: done');
