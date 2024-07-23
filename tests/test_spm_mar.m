function tests = test_spm_mar
% Unit Tests for spm_mar
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_estimation(testCase)

% Adapted from spm_mar_demo.m by Will Penny
if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox', 'spectral'));
end

% Generate MAR(2) and fit MAR model
p=2;     % model order
T=200;   % time points
w=[0;0]; % intercept

% Coeffs at lag 1
A1 = [ 0.4   1.2;   
       0.3   0.7 ];

% Coeffs at lag 2
A2 = [ 0.35 -0.3;  
      -0.40 -0.5 ];

A = [ A1 A2 ];

% Noise covaariance
C = [ 1.00  0.50;   
      0.50  1.50 ];

% Generate random timeseries x
rng('default')
x = spm_mar_gen(w, A, C, T);

% Estimate models for orders 1 to 5 and get free energy (F)
F = nan(5,1);
mar   = cell(1,5);
for m=1:5
    mar{m}=spm_mar(x,m);
    F(m) = mar{m}.fm;
end

% Log Bayes factor
F=F-min(F);

% Identify winning model
[~,winning_model_idx] = max(F);

% Test winning model
testCase.verifyEqual(winning_model_idx,2);

% Test coefficients (lag 1)
testCase.verifyEqual(A1(:),-mar{p}.lag(1).a(:),'AbsTol',0.2);

% Test coefficients (lag 2)
testCase.verifyEqual(A2(:),-mar{p}.lag(2).a(:),'AbsTol',0.2);

% Test noise covariance
testCase.verifyEqual(C(:),mar{p}.noise_cov(:),'AbsTol',0.2);