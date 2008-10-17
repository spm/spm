% Copies CTF grad and fiducials from the ooriginal dataset to an SPM file
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_copygrad.m 2350 2008-10-17 11:45:42Z vladimir $

spmfile = spm_select(1, '\.mat$', 'Select an SPM8 EEG file');

ctffile = spm_select(1, '\.*', 'Select a raw MEG data file');


hdr = fileio_read_header(ctffile);

D = spm_eeg_load(spmfile);
D = sensors(D, 'MEG', forwinv_convert_units(hdr.grad, 'mm'));
D = fiducials(D, forwinv_convert_units(fileio_read_headshape(ctffile), 'mm'));

save(D);
