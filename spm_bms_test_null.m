function [] = spm_bms_test_null (logbf_file)
% Plot PPM showing evidence for null
% FORMAT [] = spm_bms_test_null (logbf_file)
%
% logbf_file    Log Bayes Factor file providing evidence against null
%
% Call this function when SPM is already running
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_bms_test_null.m 5739 2013-11-12 17:14:13Z will $

if nargin < 1 | isempty(logbf_file)
    display('Select Log Bayes factor image providing evidence against null');
    logbf_file=spm_select;
end

[p,nam,ext,num] = spm_fileparts(logbf_file);

% Create negative logBF image
j1.input = {logbf_file};
nam='\logBF_for_null';
out_file=fullfile(p,[nam ext num]);
j1.output={out_file};
j1.expression = '-i1';
spm_imcalc(char(j1.input), char(j1.output), j1.expression);

% Create BMS mat file and PPM file
j2.dir = {[p,'\']};
j2.sess_map{1}.mod_map = j1.output;
j2.mod_name = {};
j2.method_maps = 'FFX';
j2.out_file = 0;
j2.mask = {''};
j2.nsamp = '1e6';
spm_run_bms_map(j2);

% Display Map
j3.file={[p,'\BMS.mat']};
j3.img={[p,'\m1_model_ppm.nii,1']};
j3.thres = 0.8;
j3.k = 0;
j3.scale = 0;
spm_run_bms_vis(j3);

disp('Use the threshold button to change probability threshold');
disp(' ');
disp('Use the scale button to change map from Probability to Log-Odds');