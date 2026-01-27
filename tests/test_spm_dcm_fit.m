classdef test_spm_dcm_fit < matlab.unittest.TestCase
% Unit Tests for spm_dcm_fit
%__________________________________________________________________________

% Copyright (C) 2016-2023 Wellcome Centre for Human Neuroimaging


methods (TestClassSetup)
    function setupSPM(testCase)
        spm_get_defaults('dcm.verbose',false);
    end
end % methods (TestClassSetup)


methods (Test)

function test_fmri_noparfor(testCase)

data_path = test_spm_dcm_fit.get_data_path();

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;
GCM = GCM(1:1,1);

% Remove old posteriors
for i = 1:length(GCM)
    GCM{i} = rmfield(GCM{i},'Ep');
    GCM{i} = rmfield(GCM{i},'F');
    GCM{i}.options.nograph = true;
end

% Estimate
GCM = spm_dcm_fit(GCM);

% Check DCM
testCase.assertTrue(length(GCM) == 1);
for i = 1:length(GCM)
    testCase.assertTrue(isfield(GCM{i},'Ep'));
    testCase.assertTrue(isfield(GCM{i},'F'));
end
end

% -------------------------------------------------------------------------
function test_fmri_parfor(testCase)

data_path = test_spm_dcm_fit.get_data_path();

% Load DCM (one will be enough)
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;
GCM = GCM(1:1,1);

% Remove old posteriors
for i = 1:length(GCM)
    GCM{i} = rmfield(GCM{i},'Ep');
    GCM{i} = rmfield(GCM{i},'F');
    GCM{i}.options.nograph = true;
end

% Estimate
GCM = spm_dcm_fit(GCM, true);

% Check DCM
testCase.assertTrue(length(GCM) == 1);
for i = 1:length(GCM)
    testCase.assertTrue(isfield(GCM{i},'Ep'));
    testCase.assertTrue(isfield(GCM{i},'F'));
end
end

end % methods (Test)

methods (Static, Access = private)
% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');
end
end % methods (Static, Access = private)

end % classdef