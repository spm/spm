classdef test_spm_dcm_loo < matlab.unittest.TestCase
% Unit Tests for test_spm_dcm_peb
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


methods (Test)

% -------------------------------------------------------------------------
function test_loo_group(testCase)
% Tests LOO cross-validation in the presence of a binary group effect

data_path = test_spm_dcm_loo.get_data_path();

% Reduce number of subjects to increase performance
s = [1:6 24:30];

% Load first level DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM(s,:);

% Prepare group level design matrix
X  = load(fullfile(data_path,'design_matrix.mat'));
X  = X.X(s,:);
ns = size(X,1);
X  = [ones(ns,1) X];

% Run
M.X = X;
[qE,qC,Q] = spm_dcm_loo(GCM(:,2),M,{'B'});

% Test that classical p-value is significant
[T,df] = spm_ancova(M.X(:,1:2),[],qE(:),[0;1]);
p = spm_Tcdf(-T,df(2));
testCase.assertTrue(p < 0.05);
end

% -------------------------------------------------------------------------
function test_loo_group_null(testCase)
% Tests LOO cross-validation in the absence of a binary group effect

data_path = test_spm_dcm_loo.get_data_path();

% Reduce number of subjects to increase performance
s = [1:6 24:30];

% Load first level DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM(s,:);

% Prepare group level design matrix
X  = load(fullfile(data_path,'design_matrix.mat'));
X  = X.X(s,:);
ns = size(X,1);
X  = [ones(ns,1) X];

% Run
M.X = X;
[qE,qC,Q] = spm_dcm_loo(GCM(:,2),M,{'A'});

% Test that classical p-value is significant
[T,df] = spm_ancova(M.X(:,1:2),[],qE(:),[0;1]);
p = spm_Tcdf(-T,df(2));
testCase.assertTrue(p > 0.05);
end

% -------------------------------------------------------------------------
function test_loo_continuous_null(testCase)
% Tests LOO cross-validation in the absence of a continuous group effect

data_path = test_spm_dcm_loo.get_data_path();

% Reduce number of subjects to increase performance
s = [1:6 24:30];

% Load first level DCMs
GCM = load(fullfile(data_path,'models','GCM_simulated.mat'));
GCM = GCM.GCM(s,:);

% Prepare group level design matrix
X  = load(fullfile(data_path,'design_matrix.mat'));
X  = X.X(s,:);
ns = size(X,1);
X  = [ones(ns,1) X];

% Re-order to test random age covariate
X = X(:,[1 3 2]);

% Run
M.X = X;
[qE,qC,Q] = spm_dcm_loo(GCM(:,2),M,{'B'});

% Test that classical p-value is significant
[T,df] = spm_ancova(M.X(:,1:2),[],qE(:),[0;1]);
p = spm_Tcdf(-T,df(2));
testCase.assertTrue(p >= 0.05);
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