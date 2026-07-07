classdef test_spm_cfg_getfile < matlab.unittest.TestCase
    % Unit Tests for spm_cfg_getfile
    %__________________________________________________________________________

    % Copyright (C) 2026 Department of Imaging Neuroscience


    methods (Test)

        % -------------------------------------------------------------------------
        function test_file_filter(testCase)
            % Test selecting an image file

            % Init SPM
            spm_jobman('initcfg');

            fn=fullfile(spm('dir'),'toolbox','OldSeg','grey.nii');
            out = cfg_getfile('filter',{fn},'image','.*');

            testCase.assertTrue(exist(out{1},'file') > 0);

        end
    end
end
