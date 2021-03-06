function tests = test_spm_dcm_bmr
% Unit Tests for test_spm_dcm_bmr
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_bmr(testCase)

data_path = get_data_path();

% Load PEB of full DCM
PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB;
PEB = PEB(1);

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Limit to 3 subjects for performance
GCM = GCM(1:3,:);

% Prune connections on the B-matrix only
[RCM,BMC,BMA] = spm_dcm_bmr(GCM,{'B'});

% Check outputs have expected sizes
[ns,nm] = size(GCM);
assert(all(size(RCM) == size(GCM)));
assert(all(size(BMC) == [1 ns]));
assert(all(size(BMA) == [1 ns]));

% Again with only two outputs
[RCM,BMC] = spm_dcm_bmr(GCM,{'B'});

% Check outputs have expected sizes
[ns,nm] = size(GCM);
assert(all(size(RCM) == size(GCM)));
assert(all(size(BMC) == [1 ns]));

% -------------------------------------------------------------------------
function test_bmr_single_subject(testCase)

data_path = get_data_path();

% Load PEB of full DCM
PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB;
PEB = PEB(1);

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Limit to 1 subject
GCM = GCM(1,:);

% Prune connections on the B-matrix only
[RCM,BMC,BMA] = spm_dcm_bmr(GCM,{'B'});

% Check outputs have expected sizes
[ns,nm] = size(GCM);
assert(all(size(RCM) == size(GCM)));
assert(all(size(BMC) == [1 ns]));
assert(all(size(BMA) == [1 ns]));

% Again with only two outputs
[RCM,BMC] = spm_dcm_bmr(GCM,{'B'});

% Check outputs have expected sizes
[ns,nm] = size(GCM);
assert(all(size(RCM) == size(GCM)));
assert(all(size(BMC) == [1 ns]));


% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');
