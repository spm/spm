function D = spm_eeg_weight_epochs_TF(varargin)
%This function has been deprecated
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% $Id: spm_eeg_weight_epochs_TF.m 3209 2009-06-17 11:07:47Z vladimir $

persistent runonce
if isempty(runonce)
   warning('spm_eeg_weight_epochs_TF is deprecated. Use spm_eeg_weight_epochs instead.');
   runonce = 1;
end

D = spm_eeg_weight_epochs(varargin{:});