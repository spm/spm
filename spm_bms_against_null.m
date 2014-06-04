function spm_bms_against_null(logbf_file)
% Plot PPM showing evidence against null
% FORMAT spm_bms_against_null(logbf_file)
%
% logbf_file  -  Log Bayes Factor file providing evidence against null
%
% Call this function when SPM is already running
% or set SPM to appropriate modality eg. spm('defaults','FMRI');
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_bms_against_null.m 6038 2014-06-04 15:22:42Z will $

if nargin < 1 || isempty(logbf_file)
    logbf_file = spm_select(1,'image',...
        'Select Log Bayes factor image providing evidence against null');
end

p = spm_file(logbf_file,'fpath');

% Create BMS mat file and PPM file
j2.dir         = {spm_file(logbf_file,'fpath')};
j2.sess_map{1}.mod_map = {logbf_file};
j2.mod_name    = {};
j2.method_maps = 'FFX';
j2.out_file    = 0;
j2.mask        = {''};
j2.nsamp       = '1e6';
spm_run_bms_map(j2);

% Display Map
j3.file  = {fullfile(p,'BMS.mat')};
j3.img   = {fullfile(p,['m1_model_ppm' spm_file_ext])};
j3.thres = 0.7;
j3.k     = 0;
j3.scale = 0;
spm_run_bms_vis(j3);

disp('Use the ''threshold'' button to change probability threshold');
disp(' ');
disp('Use the ''scale'' button to change map from Probability to Log-Odds');
