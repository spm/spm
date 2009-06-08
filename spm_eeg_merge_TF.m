function Dout = spm_eeg_merge_TF(S)
%This function has been deprecated
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% $Id: spm_eeg_merge_TF.m 3188 2009-06-08 08:47:46Z vladimir $

persistent runonce
if isempty(runonce)
   warning('spm_eeg_merge_TF is deprecated. Use spm_eeg_merge instead.');
   runonce = 1;
end

Dout = spm_eeg_merge(S);