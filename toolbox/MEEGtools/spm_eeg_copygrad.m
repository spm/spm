% Copies CTF grad and fiducials from the ooriginal dataset to an SPM file
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_copygrad.m 2351 2008-10-17 11:50:35Z vladimir $

spmfile = spm_select(1, '\.mat$', 'Select an SPM8 EEG file');

ctffile = spm_select(1, '\.*', 'Select a raw MEG data file');


hdr = fileio_read_header(ctffile);

D = spm_eeg_load(spmfile);
D = sensors(D, 'MEG', forwinv_convert_units(hdr.grad, 'mm'));
D = fiducials(D, forwinv_convert_units(fileio_read_headshape(ctffile), 'mm'));

% Create 2D positions for MEG (when there are no EEG sensors)
% by projecting the 3D positions to 2D
if ~isempty(strmatch('MEG', D.chantype, 'exact')) &&...
    ~isempty(D.sensors('MEG')) && isempty(D.sensors('EEG'))
    S = [];
    S.task = 'project3D';
    S.modality = 'MEG';
    S.updatehistory = 1;
    S.D = D;
    
    D = spm_eeg_prep(S);
end

save(D);
