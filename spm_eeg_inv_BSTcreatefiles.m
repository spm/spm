function D = spm_eeg_inv_BSTcreatefiles(S)

%=======================================================================
% FORMAT D = spm_eeg_inv_BSTcreatefiles(S)
% 
% This function defines the required fields for running the BrainStorm
% forward solution (bst_headmodeler.m)
%
% S		    - input struct
% (optional) fields of S:
% D			- filename of EEG/MEG mat-file
%
% Output:
% D			- EEG/MEG data struct (also written to files)
% 
% See help lines in bst_headmodeler.m for default field values
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_BSTcreatefiles.m 308 2005-11-23 19:21:56Z jeremie $

spm_defaults

try
    D = S; 
    clear S
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
end

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end



save(fullfile(P,D.fname),'D');


