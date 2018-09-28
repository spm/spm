function denoise = spm_cfg_opm_synth_gradiometer
% configuration file for performing synthetic gradiometery on OPM data
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_cfg_opm_synth_gradiometer.m 7429 2018-09-28 09:29:20Z spm $

%--------------------------------------------------------------------------
% Input Dataset
%--------------------------------------------------------------------------
D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the M/EEG mat file.'};

%--------------------------------------------------------------------------
% Confounds
%--------------------------------------------------------------------------
confounds        = cfg_entry;
confounds.tag     = 'confounds';
confounds.name    = 'Confounds';
confounds.help    = {'Labels of channel types to use for denoising. Will default to REF to use reference sensors for denoising'};
confounds.strtype = 's+';
confounds.num     = [1,100];
confounds.val     = {{'REF'}};
%--------------------------------------------------------------------------
% derivatives
%--------------------------------------------------------------------------
derivative       = cfg_menu;
derivative.tag     = 'derivative';
derivative.name    = 'Derivatives';
derivative.help    = {'Boolean to specify whether Derivatices should be used for denoising. Default is true'};
derivative.labels = {'TRUE', 'FALSE'};
derivative.values = {1, 0};
derivative.val     = {1};


%--------------------------------------------------------------------------
% simulation parameters
%--------------------------------------------------------------------------
denoise          = cfg_exbranch;
denoise.tag      = 'denoise';
denoise.name     = 'Synthetic Gradiometery';
denoise.val      = {D,confounds,derivative};
denoise.help     = {'Denoise will regress all channels of the selected type(s) from the input dataset. Optionally the derivatives of the selected type(s) can be used as well. This funciton wil automatically regress on a trial by trial basis or accross the whole sesison based on whether or not the dataset has been epoched.'}';
denoise.prog     = @synth_gradiometer;
denoise.vout     = @vout_synth_gradiometer;
denoise.modality = {'EEG'};


%==========================================================================
function out = synth_gradiometer(job)
% construct the S struct

% datset parameters
S=[];
S.D= spm_eeg_load(job.D{1});
S.confounds= job.confounds;
S.derivative = job.derivative;

% run the main function 
out.D= spm_opm_synth_gradiometer(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};


%==========================================================================
function dep = vout_synth_gradiometer(job)
% return dependencies
dep = cfg_dep;
dep.sname = 'Denoised Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Denoised Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});
