function spm_MEEGtools
% GUI gateway to MEEGtools toolbox
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_MEEGtools.m 1637 2008-05-14 14:55:21Z vladimir $

fun = spm_input('MEEG tools',1,'m', 'Copy MEG sensors|Fix Biosemi conversion', strvcat('spm_eeg_copygrad', 'spm_eeg_fixbiosemi'));
eval(fun);