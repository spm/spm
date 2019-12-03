function meeg = spm_cfg_eeg
% SPM M/EEG Configuration file for MATLABBATCH
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_eeg.m 7745 2019-12-03 14:55:56Z gareth $

%--------------------------------------------------------------------------
% M/EEG preprocessing
%--------------------------------------------------------------------------
meegprep        = cfg_choice;
meegprep.tag    = 'preproc';
meegprep.name   = 'Preprocessing';
meegprep.help   = {'M/EEG preprocessing.'};
meegprep.values = {spm_cfg_eeg_epochs spm_cfg_eeg_prepare spm_cfg_eeg_montage spm_cfg_eeg_filter...
    spm_cfg_eeg_bc spm_cfg_eeg_artefact spm_cfg_eeg_downsample spm_cfg_eeg_merge...
    spm_cfg_eeg_fuse spm_cfg_eeg_combineplanar spm_cfg_eeg_reduce spm_cfg_eeg_crop...
    spm_cfg_eeg_remove_bad_trials spm_cfg_eeg_spatial_confounds spm_cfg_eeg_correct_sensor_data}; 

%--------------------------------------------------------------------------
% M/EEG averaging
%--------------------------------------------------------------------------
meegavg        = cfg_choice;
meegavg.tag    = 'averaging';
meegavg.name   = 'Averaging';
meegavg.help   = {'M/EEG Averaging'};
meegavg.values = {spm_cfg_eeg_average spm_cfg_eeg_grandmean spm_cfg_eeg_contrast}; 

%--------------------------------------------------------------------------
% M/EEG images
%--------------------------------------------------------------------------
meegimg        = cfg_choice;
meegimg.tag    = 'images';
meegimg.name   = 'Images';
meegimg.help   = {'M/EEG Images'};
meegimg.values = {spm_cfg_eeg_convert2images spm_cfg_eeg_collapse_timefreq}; 

%--------------------------------------------------------------------------
% M/EEG time-frequency
%--------------------------------------------------------------------------
meegtf        = cfg_choice;
meegtf.tag    = 'tf';
meegtf.name   = 'Time-frequency';
meegtf.help   = {'M/EEG time-frequency.'};
meegtf.values = {spm_cfg_eeg_tf spm_cfg_eeg_tf_rescale spm_cfg_eeg_avgfreq spm_cfg_eeg_avgtime, spm_cfg_eeg_cfc}; 

%--------------------------------------------------------------------------
% M/EEG source reconstruction
%--------------------------------------------------------------------------
source        = cfg_choice;
source.tag    = 'source';
source.name   = 'Source reconstruction';
source.help   = {'M/EEG source reconstruction.'};
%source.values = { spm_cfg_eeg_inv_headmodel, spm_cfg_eeg_inv_headmodelhelmet, spm_cfg_eeg_inv_invert, spm_cfg_eeg_inv_invertiter ,spm_cfg_eeg_inv_simulate,spm_cfg_eeg_inv_mix, spm_cfg_eeg_inv_results, spm_cfg_eeg_inv_extract,spm_cfg_eeg_inv_coregshift,spm_cfg_eeg_inv_sensorshift, spm_cfg_eeg_inv_post, spm_cfg_eeg_inv_patchdef, spm_cfg_eeg_inv_prepro, spm_cfg_eeg_inv_priors,spm_cfg_eeg_inv_optimize}; 
source.values = { spm_cfg_eeg_inv_headmodel, spm_cfg_eeg_inv_headmodelhelmet, spm_cfg_eeg_inv_invert, spm_cfg_eeg_inv_invertiter ,spm_cfg_eeg_inv_simulate,spm_cfg_eeg_inv_mix, spm_cfg_eeg_inv_results, spm_cfg_eeg_inv_extract,spm_cfg_eeg_inv_coregshift,spm_cfg_eeg_inv_sensorshift, spm_cfg_eeg_dipfit, spm_cfg_eeg_momentfit}; 

%--------------------------------------------------------------------------
% M/EEG Modelling
%--------------------------------------------------------------------------
meegmodel        = cfg_choice;
meegmodel.tag    = 'modelling';
meegmodel.name   = 'Modelling';
meegmodel.help   = {'M/EEG Modelling'};
meegmodel.values = {spm_cfg_eeg_firstlevel, spm_cfg_eeg_regressors}; 
%--------------------------------------------------------------------------
% M/EEG other
%--------------------------------------------------------------------------
meegothr        = cfg_choice;
meegothr.tag    = 'other';
meegothr.name   = 'Other';
meegothr.help   = {'M/EEG Other'};
meegothr.values = {spm_cfg_eeg_review, spm_cfg_eeg_copy, spm_cfg_eeg_delete}; 

%--------------------------------------------------------------------------
% OPM
%--------------------------------------------------------------------------
meegopm        = cfg_choice;
meegopm.tag    = 'OPM';
meegopm.name   = 'OPM Preprocessing';
meegopm.help   = {'OPM Preprocessing'};
meegopm.values = {spm_cfg_opm_create,spm_cfg_opm_synth_gradiometer,spm_cfg_opm_epoch_trigger}; 

%--------------------------------------------------------------------------
% M/EEG
%--------------------------------------------------------------------------
meeg         = cfg_choice;
meeg.tag     = 'meeg';
meeg.name    = 'M/EEG';
meeg.help    = {'M/EEG functions.'};
meeg.values  = {spm_cfg_eeg_convert meegprep meegavg meegimg meegtf source meegmodel  meegopm meegothr};
