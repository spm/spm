function tests = test_spm_cfg_dcm_fmri
% Unit Tests for spm_cfg_dcm_fmri (DCM fMRI spec batch)
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_cfg_dcm_fmri.m 7479 2018-11-09 14:17:33Z peter $

tests = functiontests(localfunctions);


% -------------------------------------------------------------------------
function setup(testCase)
% Delete artefacts before each test

% Delete GCM from the tmp directory
outpath = fullfile(get_data_path(), 'tmp');
if exist(outpath,'file')
    delete(fullfile(outpath,'*.mat'));
end

% Delete DCMs from the GLM directory
glmdir = fullfile(get_data_path(), 'GLM');
dcms   = cellstr(spm_select('FPListRec',glmdir,'^DCM_.*mat$'));
for i = 1:length(dcms)
    spm_unlink(dcms{i});
end

% -------------------------------------------------------------------------
function test_specify_group(testCase)

data_path = get_data_path();

% Output path & directory
outpath = fullfile(data_path, 'tmp');
dir_out = outpath;
if ~exist(dir_out,'file')
    mkdir(dir_out);
end

% Template DCMs
fullm   = fullfile(data_path,'models','DCM_s1_m1.mat');
altm    = {fullfile(data_path,'models','DCM_s1_m2.mat');
           fullfile(data_path,'models','DCM_s1_m3.mat');
           fullfile(data_path,'models','DCM_s1_m4.mat')};

% All subjects' SPMs and VOIs
spms    = cellstr(spm_select('FPListRec',data_path,'SPM.mat'));
voi1    = cellstr(spm_select('FPListRec',data_path,'VOI_R1.mat'));
voi2    = cellstr(spm_select('FPListRec',data_path,'VOI_R2.mat'));

% Run
matlabbatch{1}.spm.dcm.spec.fmri.group.output.dir  = cellstr(dir_out);
matlabbatch{1}.spm.dcm.spec.fmri.group.output.name = 'specify_group_test';
matlabbatch{1}.spm.dcm.spec.fmri.group.template.fulldcm = cellstr(fullm);
matlabbatch{1}.spm.dcm.spec.fmri.group.template.altdcm  = altm;
matlabbatch{1}.spm.dcm.spec.fmri.group.data.spmmats      = spms;
matlabbatch{1}.spm.dcm.spec.fmri.group.data.session      = 1;
matlabbatch{1}.spm.dcm.spec.fmri.group.data.region = {voi1
                                                      voi2}';
                                                  
out = spm_jobman('run',matlabbatch);

% Load output
GCM = load(out{1}.gcmmat{1});
GCM = GCM.GCM;
GCM = spm_dcm_load(GCM);

% Check each subject's DCM timeseries match their VOIs
for s = 1:size(GCM,1)
    
    % Load VOI files
    xY1 = load(voi1{s});
    xY2 = load(voi2{s});
    Y_expected = [xY1.Y xY2.Y];
    
    % Check VOI match actual timeseries
    for m = 1:size(GCM,2)     
        Y_actual   = GCM{s,m}.Y.y;
        testCase.assertTrue(all(Y_expected(:)==Y_actual(:)), ...
            sprintf('Timeseries for subject %d model %m wrong',s,m));
    end
end

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');