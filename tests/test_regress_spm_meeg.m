classdef test_regress_spm_meeg < matlab.unittest.TestCase
    % regresion test for multimodal MEEG tutorial
    %__________________________________________________________________________

    % Copyright (C) 2026 Department of Imaging Neuroscience

    methods (TestClassSetup)
        function setupSPM(testCase)
            spm('defaults','eeg');
        end
    end % methods (TestClassSetup)

    methods (Test)

        function test_regress_spm_meeg_1(testCase)

            rawpth = fullfile(spm('Dir'),'tests','data','multimodal_faces');
                       

            nrun = 2; % enter the number of runs here

            jobfile = {fullfile(spm('Dir'), 'tests', 'batch_preproc_meeg_pipeline_job.m')};
            jobs = repmat(jobfile, 1, nrun);
            inputs = cell(3, nrun);
            inputs2 = {};

            for crun = 1:nrun
                inputs{1, crun} = cellstr(fullfile(rawpth, 'meg',sprintf('run_%02d_sss.fif',crun)));
                inputs{2, crun} = cellstr(fullfile(rawpth, 'meg','bad_channels.mat'));
                inputs{3, crun} = cellstr(fullfile(rawpth, 'meg','Trials',sprintf('run_%02d_trldef.mat',crun)));
                inputs2{crun, 1} = sprintf('befdspmeeg_run_%02d_sss.mat', crun);
            end

            spm_jobman('run', jobs, inputs{:});
            jobfile = {fullfile(spm('Dir'), 'tests', 'batch_preproc_meeg_pipeline_merged_job.m')};
            spm_jobman('run', jobfile, inputs2);
            jobfile = {fullfile(spm('Dir'), 'tests', 'batch_evoked_meeg.m')};
            spm_jobman('run', jobfile, {'Mcbefdspmeeg_run_01_sss.mat'});
            jobfile = {fullfile(spm('Dir'), 'tests', 'batch_preproc_meeg_tf_job.m')};
            spm_jobman('run', jobfile, {'Mcbefdspmeeg_run_01_sss.mat'});

            jobfile = {fullfile(spm('Dir'), 'tests', 'batch_localise_job.m')};
            inputs = cell(5, 1);
            inputs{1, 1} = {'apMcbefdspmeeg_run_01_sss.mat'};
            inputs{2, 1} = {fullfile(rawpth, 'anat', 'sub-15-T1w.nii')};
            fid = load(fullfile(rawpth, 'anat', 'mri_fids.mat'));
            inputs{3, 1} = fid.mri_fids(1, :);
            inputs{4, 1} = fid.mri_fids(2, :);
            inputs{5, 1} = fid.mri_fids(3, :);
            spm_jobman('run', jobfile, inputs{:});
            
            D = spm_eeg_load('apMcbefdspmeeg_run_01_sss.mat');          

            outcome = ~isnan(D.inv{1}.inverse.R2) & D.inv{1}.inverse.R2>98 &...
                ~isnan(D.inv{2}.inverse.R2) & D.inv{2}.inverse.R2>97;

            testCase.verifyTrue(outcome);
           
            % clean up
            delete('*.mat');
            delete('*.dat');
            delete('*.nii');
            delete(fullfile(rawpth, 'anat', '*.gii'));
            delete(fullfile(rawpth, 'anat', '*BEM*'));
            delete(fullfile(rawpth, 'anat', '*seg8*'));
            delete(fullfile(rawpth, 'anat', 'y_*.nii'));
            rmdir('*_sss', 's');

        end
    end % methods (Test)

end % classdef