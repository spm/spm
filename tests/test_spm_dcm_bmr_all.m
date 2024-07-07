function tests = test_spm_dcm_bmr_all
% Unit Tests for test_spm_dcm_bmr_all
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_bmr_all(testCase)

data_path = get_data_path();

% Load PEB of full DCM
PEB = load(fullfile(data_path,'PEB_test.mat'));
PEB = PEB.PEB;
PEB = PEB(1);

% Load DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Prune connections on the B-matrix only
[rPEB,BMR,BMA] = spm_dcm_bmr_all(PEB,{'B'},'BMA');

% Check the only surviving B-matrix (group difference) parameter is R1->R2
b_reduced = BMA.Cp(:,2);
expected  = logical([0 0 0 0 1 0]);
testCase.assertTrue(all(b_reduced(~expected) < 1e-4));
testCase.assertTrue(all(b_reduced(expected) > 1e-4));

Pp = full(spm_unvec(rPEB.Pp,PEB.Ep));
testCase.assertTrue(Pp(5,2) > 0.9);

% -------------------------------------------------------------------------
% function test_divide(testCase)
% % This test illustrates a suspected Matlab bug that occurs specifically  
% % on Mac with MATLAB 2023b
% Pk(1,:) =  [0.0264    0.0312    0.0312    0.0243    0.0312    0.0312];
% Pk(2,:) =  [0.0049    0.0000    0.0000    0.0069    0.0000    0.0001];
%    
% Pk = Pk(1,:)./sum(Pk);
% disp(Pk);
% 
% testCase.assertTrue(~all(abs(Pk) < 0.001));

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');
