function tests = test_spm_dcm_peb
% Unit Tests for test_spm_dcm_peb
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_peb.m 7479 2018-11-09 14:17:33Z peter $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_peb(testCase)

data_path = get_data_path();

% Load first level DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;
nm  = size(GCM,2);

% Prepare group level design matrix
X = load(fullfile(data_path,'design_matrix.mat'));
X = X.X;

ns = size(X,1);
X  = [ones(ns,1) X];
nx = size(X,2);

% Number of free DCM parameters expected in full model
np = 6;

% Estimate PEB
PEB = spm_dcm_peb(GCM, X);

% Check output sizes
testCase.assertEqual(size(PEB,2),           nm);
testCase.assertEqual(size(PEB(1).Xnames,2), nx);
testCase.assertEqual(size(PEB(1).Ep,2),     nx);
testCase.assertEqual(size(PEB(1).Ep,1),     np);
testCase.assertEqual(size(PEB(1).Pnames,1), np);
testCase.assertEqual(size(PEB(1).Pind,1),   np);
testCase.assertEqual(size(PEB(1).Ce),       [np np]);
testCase.assertEqual(size(PEB(1).Cp),       [np*nx np*nx]);

% % Save test PEB for other tests
% %save(fullfile(data_path,'PEB_test.mat'),'PEB');
% -------------------------------------------------------------------------
function test_precision_def(testCase)
% Tests different configurations of precision components

data_path = get_data_path();

% Load first level DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Get the indices of free parameters
q   = spm_find_pC(GCM{1},{'A','B'});

% Prepare group level design matrix
X  = load(fullfile(data_path,'design_matrix.mat'));
X  = X.X;
ns = size(X,1);
X  = [ones(ns,1) X];
M.X = X;

% Estimate PEB (single precision component)
M.Q = 'single';
PEB = spm_dcm_peb(GCM(:,1), M);
testCase.assertEqual(length(PEB.Eh), 1);
testCase.assertEqual(length(PEB.Ch), 1);

% Estimate PEB (one precision component per parameter)
M.Q = 'all';
PEB = spm_dcm_peb(GCM(:,1), M);
testCase.assertEqual(length(PEB.Eh), length(q));
testCase.assertEqual(length(PEB.Ch), length(q));

% Estimate PEB (one precision component per field)
M.Q = 'fields';
PEB = spm_dcm_peb(GCM(:,1), M);
testCase.assertEqual(length(PEB.Eh), 2);
testCase.assertEqual(length(PEB.Ch), 2);

% Estimate PEB (manual precision component definition)
np   = size(GCM{1}.M.pC,1);
Q{1} = sparse(q(1:2),q(1:2),[1 1],np,np);
Q{2} = sparse(q(3:6),q(3:6),[1 1 1 1],np,np);
M.Q = Q;
PEB = spm_dcm_peb(GCM(:,1), M);
testCase.assertEqual(length(PEB.Eh), 2);
testCase.assertEqual(length(PEB.Ch), 2);

% -------------------------------------------------------------------------
function test_peb_of_pebs(testCase)
% Tests hierarchical models with PEB as the input to other PEBs

data_path = get_data_path();

% Load first level DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM;

% Prepare group level design matrix
X = load(fullfile(data_path,'design_matrix.mat'));
X = X.X;
ns = size(X,1);
X  = [ones(ns,1) X];

% Assign subjects to groups
g1 = (X(:,2) == -1);
g2 = (X(:,2) == 1);

% Estimate PEB
M = struct();
M.Q = 'all';
PEB1 = spm_dcm_peb(GCM(g1,1),M);
PEB2 = spm_dcm_peb(GCM(g2,1),M);

% Group model
M = struct();
M.X = [1 1; 1 -1];
M.Q = 'none';
PEB  = spm_dcm_peb({PEB1;PEB2},M);

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region');