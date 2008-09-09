% Recompiles the Biosemi mex file
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_fixbiosemi.m 2060 2008-09-09 17:34:21Z guillaume $

swd = pwd;
cd(fullfile(spm('dir'),'external','fileio','private','mex'));
eval(['mex read_24bit.c -outdir ' fullfile(spm('dir'),'external','fileio','private')]);
cd(swd);