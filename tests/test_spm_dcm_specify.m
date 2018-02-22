function tests = test_spm_dcm_specify
% Unit Tests for spm_dcm_specify_ui
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_specify.m 7264 2018-02-22 14:43:47Z peter $

tests = functiontests(localfunctions);


function test_ui(testCase)

data_path = get_data_path();

% Identify SPM & VOIs
glm_dir = fullfile(data_path,'GLM','S001');
SPM = load(fullfile(glm_dir,'SPM.mat'));
SPM = SPM.SPM;

xY = {fullfile(glm_dir,'VOI_R1.mat');
      fullfile(glm_dir,'VOI_R2.mat')};

% Specify connectivity
a = ones(2,2);
b = zeros(2,2,2);
c = [1 0; 1 0];
d = ones(2,2,2);

b(:,:,2) = eye(2);

% DCM settings
s = struct();
s.name       = 'tmp';
s.u          = [1 1]';
s.delays     = [1.2 1.2];
s.TE         = 0.05;
s.nonlinear  = true;
s.two_state  = true;
s.stochastic = true;
s.centre     = true;
s.induced    = 0;
s.a          = a;
s.b          = b;
s.c          = c;
s.d          = d;

% Test
DCM = spm_dcm_specify(SPM,xY,s);

% Check
testCase.assertEqual(length(DCM.U.name),2);
testCase.assertEqual(DCM.delays,s.delays);
testCase.assertEqual(DCM.TE,s.TE);
testCase.assertEqual(DCM.options.nonlinear,s.nonlinear);
testCase.assertEqual(DCM.options.two_state,s.two_state);
testCase.assertEqual(DCM.options.stochastic,s.stochastic);
testCase.assertEqual(DCM.options.centre,s.centre);
testCase.assertEqual(DCM.options.induced,s.induced);
testCase.assertEqual(DCM.a,s.a);
testCase.assertEqual(DCM.b,s.b);
testCase.assertEqual(DCM.c,s.c);
testCase.assertEqual(DCM.d,s.d);

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');