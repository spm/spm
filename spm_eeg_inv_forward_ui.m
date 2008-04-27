function D = spm_eeg_inv_forward_ui(varargin)
% Forward Solution user-interface routine
% commands the forward computation for either EEG or MEG data
% and calls for various types of solutions using BrainStorm functions
% as well as a realistic sphere solution (for EEG).
%
% FORMAT D = spm_eeg_inv_forward_ui(D,val)
% Input:
% D         - input data struct (optional)
% Output:
% D         - same data struct including the forward solution files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_forward_ui.m 1488 2008-04-27 14:11:48Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% get method
%--------------------------------------------------------------------------
D.inv{D.val}.method = 'Imaging';

% compute forward model
%==========================================================================
D = spm_eeg_inv_forward(D);

fprintf('Foward model complete - thank you\n')
