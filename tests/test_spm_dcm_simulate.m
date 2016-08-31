function tests = test_spm_dcm_simulate
% Unit Tests for test_spm_dcm_simulate
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_simulate.m 6865 2016-08-31 16:49:57Z guillaume $

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_simulate_snr_var_fmri(testCase)
% Test simulation based on SNR variance
data_path = get_data_path();

% Load DCMs
GCM = load(fullfile(data_path,'GCM_simulated.mat'));
GCM = GCM.GCM;

% Simulate data
SNR     = 1;
gen_idx = 1;
rng('default');rng(1);
%st = randn('state');
%randn('state',100);
[GCM,gen] = spm_dcm_simulate(GCM, 'SNR_var', SNR, gen_idx);
%randn('state',st);

ns = size(GCM,1);

% Gather actual SNR
actual_snr = [];
for s = 1:ns
    DCM_gen = gen{s};
    
    signal  = DCM_gen.y;
    noise   = DCM_gen.Y.y - signal; 
    
    actual_snr(s,:) = var(signal) ./ var(noise);
end

nr = size(actual_snr,2); % Number of regions

% Check that there isn't significant evidence that the SNR deviated from
% the expected value (one-sample ttest)
for n = 1:nr    
    data = actual_snr(:,n);
    p    = one_sample_tt(data,SNR);
    testCase.assertTrue(p >= 0.05);
end

% -------------------------------------------------------------------------
function test_simulate_snr_std_fmri(testCase)
% Test simulation based on SNR standard deviatoin
data_path = get_data_path();

% Load DCMs
GCM = load(fullfile(data_path,'GCM_simulated.mat'));
GCM = GCM.GCM;

% Simulate data
SNR     = 1;
gen_idx = 1;
rng('default');rng(1);
%st = randn('state');
%randn('state',100);
[GCM,gen] = spm_dcm_simulate(GCM, 'SNR_std', SNR, gen_idx);
%randn('state',st);

ns = size(GCM,1);

% Gather actual SNR
actual_snr = [];
for s = 1:ns
    DCM_gen = gen{s};
    
    signal  = DCM_gen.y;
    noise   = DCM_gen.Y.y - signal; 
    actual_snr(s,:) = std(signal) ./ std(noise);
end

nr = size(actual_snr,2); % Number of regions

% Check that there isn't significant evidence that the SNR deviated from
% the expected value (one-sample ttest)
for n = 1:nr    
    data = actual_snr(:,n);
    p    = one_sample_tt(data,SNR);
    testCase.assertTrue(p >= 0.05);
end

% -------------------------------------------------------------------------
function test_simulate_snr_var(testCase)
% Test simulation based on fixed noise variance
data_path = get_data_path();

% Load DCMs
GCM = load(fullfile(data_path,'GCM_simulated.mat'));
GCM = GCM.GCM;

% Simulate data
noise_var = 1;
gen_idx   = 1;
rng('default');rng(1);
%st = randn('state');
%randn('state',100);
[GCM,gen] = spm_dcm_simulate(GCM, 'var', noise_var, gen_idx);
%randn('state',st);

ns = size(GCM,1);

% Gather actual SNR
actual_var = [];
for s = 1:ns
    DCM_gen = gen{s};
    
    signal  = DCM_gen.y;
    noise   = DCM_gen.Y.y - signal; 
    
    actual_var(s,:) = var(noise);
end

nr = size(actual_var,2); % Number of regions

% Check that there isn't significant evidence that the SNR deviated from
% the expected value (one-sample ttest)
for n = 1:nr    
    data = actual_var(:,n);
    p    = one_sample_tt(data,noise_var);
    testCase.assertTrue(p >= 0.05);
end

% -------------------------------------------------------------------------
function p = one_sample_tt(data,mu)
% One sample t-test
t = (mean(data)-mu) / (std(data) / sqrt(length(data)));
p = (1-spm_Tcdf(t,length(data)-1)) * 2;

% -------------------------------------------------------------------------
function data_path = get_data_path()

data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'fMRI', 'simulated_2region', 'models');
