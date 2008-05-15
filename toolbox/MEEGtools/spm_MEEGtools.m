function spm_MEEGtools
% GUI gateway to MEEGtools toolbox
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_MEEGtools.m 1669 2008-05-15 18:34:49Z vladimir $

fun = spm_input('MEEG tools',1,'m', 'Copy MEG sensors|Fix Biosemi conversion|Interactive plotting',...
    strvcat('spm_eeg_copygrad', 'spm_eeg_fixbiosemi', 'spm_eeg_plot_interactive'));
eval(fun);