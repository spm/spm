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
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_copy.m 5079 2012-11-25 18:38:18Z vladimir $

SVNrev = '$Rev: 5079 $';


%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG copy');

if ~isfield(S, 'updatehistory'),       S.updatehistory = 1;           end

D       = spm_eeg_load(S.D);


D       = copy(D, S.outfile);

if ~isfield(S, 'updatehistory') || S.updatehistory
    D = D.history('spm_eeg_copy', S); 
end

save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG copy: done');