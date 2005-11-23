function D = spm_eeg_inv_forward_ui(S)

%=======================================================================
% Forward Solution user-interface routine
% commands the forward computation for either EEG or MEG data
% and calls for various types of solutions using BrainStorm functions
% as well as a realistic sphere solution (for EEG).
%
% FORMAT D = spm_eeg_inv_forward_ui(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the forward solution files and variables
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id$

spm_defaults

try
    D = S; 
    clear S
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
end

val = length(D.inv);
if val > 1
    D = spm_eeg_inv_copyfields(D,[0 0 1 0]);
end

D = spm_eeg_inv_BSTfwdsol(D);

save(fullfile(D.path,D.fname),'D');

return