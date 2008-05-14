% Recompiles the Biosemi mex file
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_fixbiosemi.m 1633 2008-05-14 11:23:23Z vladimir $

cd(fullfile(spm('dir'), 'external\fileio\private\mex'));
eval(['mex read_24bit.c -outdir ' fullfile(spm('dir'), 'external\fileio\private')]);