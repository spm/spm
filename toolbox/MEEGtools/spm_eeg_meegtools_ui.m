function spm_eeg_meegtools_ui
% GUI gateway to MEEGtools toolbox
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_meegtools_ui.m 1633 2008-05-14 11:23:23Z vladimir $

fun = spm_input('MEEG tools',1,'m', 'Copy MEG sensors|Fix Biosemi conversion', strvcat('spm_eeg_copygrad', 'spm_eeg_fixbiosemi'));
eval(fun);