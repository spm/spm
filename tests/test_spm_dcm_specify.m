function tests = test_spm_dcm_specify
% Unit Tests for spm_dcm_specify_ui
%__________________________________________________________________________

% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_inputs_by_name(testCase)

% Identify SPM
glm_dir = get_data_path();
SPM = load(fullfile(glm_dir,'SPM.mat'));
SPM = SPM.SPM;

% Select experimental conditions to import from the SPM.mat
cond = struct();
cond(1).name    = 'Faces';                  % Name for use in the DCM
cond(1).spmname = {'N1', 'N2', 'F1', 'F2'}; % SPM conditions (SPM.Sess(s).U)

cond(2).name    = 'Fame';
cond(2).spmname = {'F1','F2'};

cond(3).name    = 'Repetition';
cond(3).spmname = {'N2','F2'};

cond(4).name    = 'Lag';
cond(4).spmname = {'N2xLag^1','F2xLag^1'};

cond(5).name    = 'Lag quadratic';
cond(5).spmname = {'N2xLag^2','F2xLag^2'};

% Specify DCM
DCM = run_specify(cond,SPM,testCase);

% Combine conditions manually and check they match those in the DCM
k = 33:(size(SPM.Sess.U(1).u,1));

testCase.assertEqual(DCM.U.u(:,1), ...
    double(SPM.Sess.U(1).u(k,1) | SPM.Sess.U(2).u(k,1) | ...
    SPM.Sess.U(3).u(k,1) | SPM.Sess.U(4).u(k,1)));

testCase.assertEqual(DCM.U.u(:,2), ...
    double(SPM.Sess.U(3).u(k,1) | SPM.Sess.U(4).u(k,1)));

testCase.assertEqual(DCM.U.u(:,3), ...
    double(SPM.Sess.U(2).u(k,1) | SPM.Sess.U(4).u(k,1)));

testCase.assertEqual(DCM.U.u(:,4), ...
    double(SPM.Sess.U(2).u(k,2) | SPM.Sess.U(4).u(k,2)));

testCase.assertEqual(DCM.U.u(:,5), ...
    double(SPM.Sess.U(2).u(k,3) | SPM.Sess.U(4).u(k,3)));


% -------------------------------------------------------------------------
function test_inputs_by_ID(testCase)

% Identify SPM
glm_dir = get_data_path();
SPM = load(fullfile(glm_dir,'SPM.mat'));
SPM = SPM.SPM;

u = [1 0 0
     0 1 0
     1 0 0
     0 0 1];
DCM = run_specify(u,SPM,testCase);

% Get conditions manually and check they match those in the DCM
k = 33:(size(SPM.Sess.U(1).u,1));
testCase.assertEqual(DCM.U.u(:,1), SPM.Sess.U(1).u(k,1));
testCase.assertEqual(DCM.U.u(:,2), SPM.Sess.U(2).u(k,2));
testCase.assertEqual(DCM.U.u(:,3), SPM.Sess.U(3).u(k,1));
testCase.assertEqual(DCM.U.u(:,4), SPM.Sess.U(4).u(k,3));

% -------------------------------------------------------------------------
function DCM = run_specify(spec,SPM,testCase)
% input 'spec' determines how the DCM input is specified. Can be a 
% structure or matrix of condition indices.

glm_dir = get_data_path();
xY = {fullfile(glm_dir,'VOI_lFus_1.mat');
      fullfile(glm_dir,'VOI_rFus_1.mat')};

% Number of regions
n = 2;

% Number of conditions in the DCM
ncond = length(spec);

% Connect both regions
a = ones(n,n);

% All conditions except the first modulate the self-connections
b = zeros(n,n,ncond);
for i = 2:ncond
    b(:,:,i) = eye(n);
end

% Faces drive both regions
c = [1 0 0 0 0; 
     1 0 0 0 0];
 
% No non-linear coupling
d = zeros(n,n,n);

% DCM settings
s = struct();
s.name       = 'tmp';
s.delays     = repmat(SPM.xY.RT/2,1,n);
s.TE         = 0.05;
s.nonlinear  = false;
s.two_state  = false;
s.stochastic = false;
s.centre     = true;
s.induced    = 0;
s.a          = a;
s.b          = b;
s.c          = c;
s.d          = d;

if isstruct(spec)
    s.cond = spec;
elseif isnumeric(spec)
    s.u = spec;
end

% Test
DCM = spm_dcm_specify(SPM,xY,s);

% Check
testCase.assertEqual(length(DCM.U.name),ncond);
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
    'data', 'test_spm_dcm_specify');
