function spm_MEEGtools
% GUI gateway to MEEGtools toolbox
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_MEEGtools.m 6029 2014-05-30 18:52:03Z vladimir $


funlist = {
    'Transform EEG cap', 'spm_eeg_transform_cap';
    'Re-reference EEG', 'spm_eeg_reref_eeg';
    'Split conditions into separate datasets', 'spm_eeg_split_conditions';   
    'Fieldtrip interactive plotting', 'spm_eeg_plot_interactive';
    'Fieldtrip visual artefact rejection', 'spm_eeg_ft_artefact_visual';
    'Fieldtrip dipole fitting', 'spm_eeg_ft_dipolefitting';
    'Vector-AR connectivity measures', 'spm_eeg_var_measures';  
    'Use CTF head localization' , 'spm_eeg_megheadloc'
    'Fix CTF head position data' ,'spm_eeg_fix_ctf_headloc'
    'Fieldtrip manual coregistration' , 'spm_eeg_ft_datareg_manual'
    'Remove spikes from EEG' , 'spm_eeg_remove_spikes'
    'Reduce jumps in MEG data' , 'spm_eeg_remove_jumps'
    'Detrending and Hanning for ERPs', 'spm_eeg_erp_correction'
    'Extract dipole waveforms', 'spm_eeg_dipole_waveforms'
    'Fieldtrip-SPM robust multitaper coherence', 'spm_eeg_ft_multitaper_coherence'
    'Fieldtrip multitaper power map', 'spm_eeg_ft_multitaper_powermap'
    'Interpolate artefact segment', 'spm_eeg_interpolate_artefact'     
    'Relabel trials for epoched CTF datasets', 'spm_eeg_recode_epoched_ctf'
    'Correct TMS artefact', 'spm_eeg_tms_correct'
    'Plot scalp maps from M/EEG image', 'spm_eeg_img2maps'
    'Continuous data power', 'spm_eeg_cont_power'
    };

str = sprintf('%s|', funlist{:, 1});
str = str(1:(end-1));

fun = spm_input('MEEG tools',1,'m', str, strvcat(funlist(:, 2)));
  
eval(fun);