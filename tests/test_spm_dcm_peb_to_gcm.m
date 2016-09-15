function tests = test_spm_dcm_peb_to_gcm
% Unit Tests for spm_dcm_peb_to_gcm
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_dcm_peb_to_gcm.m 6875 2016-09-15 09:30:18Z peter $

tests = functiontests(localfunctions);

%--------------------------------------------------------------------------
function data_path = get_data_path()
data_path = fullfile( spm('Dir'), 'tests', ...
    'data', 'test_spm_dcm_peb_to_gcm');

%--------------------------------------------------------------------------
function test_2groups(testCase)
% Tests creating GCMs for two groups of subjects who differ in a single
% parameter

% Settings
n     = 100;        % Subjects
b_Hz  = 0.6;        % Strength of group difference
param = 'B(2,1,3)'; % Parameter to vary between groups
beta  = 16;         % Within:between variance ratio for simulated parameters

rng(1);

% Load template DCM
DCM_template = fullfile(get_data_path(), 'DCM_attention.mat');
DCM_template = load(DCM_template);
DCM_template = DCM_template.DCM;

% Indices of subjects in each group
group    = cell(1,2);
group{1} = 1:ceil(n/2);
group{2} = ceil(n/2)+1:n;

% PEB design matrix
X             = ones(n,2);
X(group{2},2) = -1;

% Create dummy PEB
PEB = create_peb(DCM_template, X, b_Hz, param, beta);

% Create GCM array
options = struct('nsubjects',n);
[GCM,PEB] = spm_dcm_peb_to_gcm(PEB,DCM_template,options);

% Check size
testCase.assertTrue( all(size(GCM)==[n 1]) );

% Get group arithmetic average of parameters
Eps = [];
for s = 1:n
    Eps(:,s) = spm_vec(GCM{s,1}.Ep);
end
Eps = Eps(PEB.Pind,:);

% Check subject' parameters are similar to PEB.Ep(:,1)
testCase.assertGreaterThanOrEqual(...
    one_sample_tt(mean(Eps,2)-PEB.Ep(:,1)),...
    0.05);

% Check group difference in B-parameter is similar to PEB.Ep(idx,2)
idx       = spm_fieldindices(GCM{1}.M.pE,param);
common    = PEB.Ep(PEB.Pind == idx,1);
attention = PEB.Ep(PEB.Pind == idx,2);
expected  = [common+attention common-attention]';

for g = 1:length(group)    
    % Gather group's parameters
    Eps = [];
    for s = 1:length(group{g})
        Eps(:,s) = spm_vec(GCM{group{g}(s),1}.Ep);
    end
    
    % Check group mean is similar to PEB parameter
    testCase.assertGreaterThanOrEqual(...
        one_sample_tt(Eps(idx,:), expected(g)),...
        0.05);
    
    % Check that between-subjects variance of parameters is as expected
    testCase.assertGreaterThanOrEqual(...
        one_sample_tt(diag(PEB.Ce) - var(Eps(PEB.Pind,:),[],2)),...
        0.05);    
    
end

%--------------------------------------------------------------------------
function PEB = create_peb(DCM, X, b_Hz, param, beta)
% Generate dummy PEB for simulation
%
% DCM   - Template DCM
% X     - PEB design matrix [parameters x regressors]
% b_Hz  - Group difference on a connection in Hz
% param - Parameter name to vary between groups
% beta  - Within:between variance ratio for simulated parameters

% Fields to include in the PEB
fields = {'A','B','C','transit','decay','epsilon'};

% Get within-subject priors from template DCM
[pE,pC] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d);

% Sample once from priors to get 'true' group-level parameter means
pE_struct = pE;
pE        = full(spm_vec(pE));
gP        = spm_normrnd(pE,pC,1);

% Force C-matrix (driving input) to be positive
idx = spm_fieldindices(pE_struct,'C');
gP(idx) = abs(gP(idx));

% Get indices of parameters to include in the PEB
Pind = spm_find_pC(pC,pE_struct,fields);

% Create group difference regressor
att = zeros(length(pE),1);
idx = spm_fieldindices(pE_struct,param);
att(idx) = b_Hz;

% Collate PEB
clear PEB;
PEB.Ep   = [gP(Pind) att(Pind)];
PEB.beta = beta;
PEB.Pind = Pind;
PEB.M.X  = X;

% -------------------------------------------------------------------------
function p = one_sample_tt(data,mu)
% One sample t-test

if nargin < 2, mu = 0; end
t = (mean(data)-mu) / (std(data) / sqrt(length(data)));
p = (1-spm_Tcdf(t,length(data)-1)) * 2;