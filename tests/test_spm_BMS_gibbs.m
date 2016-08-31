function tests = test_spm_BMS_gibbs
% Unit Tests for spm_BMS_gibbs
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_BMS_gibbs.m 6865 2016-08-31 16:49:57Z guillaume $

tests = functiontests(localfunctions);


function test_strong_evidence(testCase)

% Setup
import matlab.unittest.constraints.*
rng('default');
rng(1);

n = 20;      % Subjects
models = 4;  % Models per subject

% Frequency of each model in the group
r = [0.60 0.20 0.15 0.05];

num_subjects_per_model = n .* r;

% Build subjects x models matrix of model assignments
m = [];
for model = 1:models  
    row = zeros(1,models); 
    row(model) = 1;
    
    m = [m; repmat(row, num_subjects_per_model(model), 1)];
end

% Free energy of the generative model in each subject
mean_F = -10; 
std_F  = 0.5;
F_gen  = mean_F + std_F .* randn(sum(m(:) > 0),1);

% Build free energy matrix (non-generative models: mean -13, STD 1)
F    = -13 + randn(n,4);
F(m > 0) = F_gen;

% Run
[exp_r,xp] = spm_BMS_gibbs(F);

% Check expected frequencies are within 10%
actual = exp_r;
testCase.verifyThat(actual, IsEqualTo(r, 'Within', AbsoluteTolerance(0.1) ) );

% Check exceedance probability
actual = xp(1);
testCase.verifyThat(actual, IsGreaterThanOrEqualTo(0.95) );
