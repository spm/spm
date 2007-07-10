function D = spm_eeg_inv_forward_ui(varargin)

%==========================================================================
% Forward Solution user-interface routine
% commands the forward computation for either EEG or MEG data
% and calls for various types of solutions using BrainStorm functions
% as well as a realistic sphere solution (for EEG).
%
% FORMAT D = spm_eeg_inv_forward_ui(D,val)
% Input:
% D		    - input data struct (optional)
% Output:
% D			- same data struct including the forward solution files and variables
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_forward_ui.m 849 2007-07-10 15:30:31Z rik $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% get method
%--------------------------------------------------------------------------
if strcmp(D.modality,'MEG')
    method = 'Imaging';
else
    try 
	method = D.inv{val}.method;
    catch
	method = questdlg('recontruction','Please select','Imaging','ECD','Imaging');
    end
end
D.inv{D.val}.method = method;

% compute forward model
%==========================================================================
if strcmp(method,'Imaging')
    D = spm_eeg_inv_BSTfwdsol(D);
else
    D = spm_eeg_inv_elec_Rsph_ui(D);
end

fprintf('Foward model complete - thank you\n')
