function spm_bms_test_null(logbf_file)
% Plot PPM showing evidence for null
% FORMAT spm_bms_test_null(logbf_file)
%
% logbf_file  -  Log Bayes Factor file providing evidence against null
%
% Call this function when SPM is already running
% or set SPM to appropriate modality eg. spm('defaults','FMRI');
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2013-2022 Wellcome Centre for Human Neuroimaging


if nargin < 1 || isempty(logbf_file)
    logbf_file = spm_select(1,'image',...
        'Select Log Bayes factor image providing evidence against null');
end

p = spm_file(logbf_file,'fpath');

% Create negative logBF image
j1.input      = {logbf_file};
j1.output     = {spm_file(logbf_file,'basename','logBF_for_null')};
j1.expression = '-i1';
j1.flags      = {0,0,0,spm_type('float32')};
spm_imcalc(char(j1.input), char(j1.output), j1.expression, j1.flags);

% Create BMS mat file and PPM file
j2.dir         = {spm_file(logbf_file,'fpath')};
j2.sess_map{1}.mod_map = j1.output;
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
