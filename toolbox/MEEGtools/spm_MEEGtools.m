function spm_MEEGtools
% GUI gateway to MEEGtools toolbox
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_MEEGtools.m 2333 2008-10-13 13:19:22Z vladimir $


funlist = {
    'Copy MEG sensors', 'spm_eeg_copygrad';
    'Fix Biosemi conversion', 'spm_eeg_fixbiosemi';
    'Fieldtrip interactive plotting', 'spm_eeg_plot_interactive';
    'Fieldtrip visual artefact rejection', 'spm_eeg_ft_artefact_visual';
    'Fieldtrip dipole fitting', 'spm_eeg_ft_dipolefitting';
    'Vector-AR connectivity measures', 'spm_eeg_var_measures';
    'Fieldtrip-SPM MEG head modelling' , 'spm_eeg_ft_indiv_meg_model'
    'Define spatial confounds' , 'spm_eeg_spatial_confounds'
    };

str = sprintf('%s|', funlist{:, 1});
str = str(1:(end-1));

fun = spm_input('MEEG tools',1,'m', str, strvcat(funlist(:, 2)));
  
eval(fun);