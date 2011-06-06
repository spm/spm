function spm_MEEGtools
% GUI gateway to MEEGtools toolbox
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_MEEGtools.m 4342 2011-06-06 12:44:56Z vladimir $


funlist = {
    'Copy MEG sensors', 'spm_eeg_copygrad';
    'Transform EEG cap', 'spm_eeg_transform_cap';
    'Re-reference EEG', 'spm_eeg_reref_eeg';
    'Split conditions into separate datasets', 'spm_eeg_split_conditions';
    'Estimate multiple DCMs', 'spm_dcm_estimate_group';
    'Fieldtrip interactive plotting', 'spm_eeg_plot_interactive';
    'Fieldtrip visual artefact rejection', 'spm_eeg_ft_artefact_visual';
    'Fieldtrip dipole fitting', 'spm_eeg_ft_dipolefitting';
    'Vector-AR connectivity measures', 'spm_eeg_var_measures';
    'Define spatial confounds' , 'spm_eeg_spatial_confounds'
    'Correct sensor data',        'spm_eeg_correct_sensor_data'
    'Use CTF head localization' , 'spm_eeg_megheadloc'
    'Fix CTF head position data' ,'spm_eeg_fix_ctf_headloc'
    'Fieldtrip manual coregistration' , 'spm_eeg_ft_datareg_manual'
    'Remove spikes from EEG' , 'spm_eeg_remove_spikes'
    'Reduce jumps in MEG data' , 'spm_eeg_remove_jumps'
    'Detrending and Hanning for ERPs', 'spm_eeg_erp_correction'
    'Extract dipole waveforms', 'spm_eeg_dipole_waveforms'
    'Fieldtrip multitaper TF', 'spm_eeg_ft_multitaper_tf'
    'Fieldtrip-SPM robust multitaper coherence', 'spm_eeg_ft_multitaper_coherence'
    'Fieldtrip multitaper power map', 'spm_eeg_ft_multitaper_powermap'
    'Interpolate artefact segment', 'spm_eeg_interpolate_artefact'
    'FMRIB Detect ECG peaks',   'spm_eeg_fmrib_qrsdetect'     
    'Detect eyeblinks',  'spm_eeg_detect_eyeblinks'
    'Relabel trials for epoched CTF datasets', 'spm_eeg_recode_epoched_ctf'
    'Correct TMS artefact', 'spm_eeg_tms_correct'
    'Plot scalp maps from M/EEG image', 'spm_eeg_img2maps'
    'Continuous data power', 'spm_eeg_cont_power'
    };

str = sprintf('%s|', funlist{:, 1});
str = str(1:(end-1));

fun = spm_input('MEEG tools',1,'m', str, strvcat(funlist(:, 2)));
  
eval(fun);