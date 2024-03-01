function tests = test_spm_bireduce
% Unit Tests for test_spm_bireduce
%__________________________________________________________________________

% Copyright (C) 2016-2024 Wellcome Centre for Human Neuroimaging

tests = functiontests(localfunctions);

% -------------------------------------------------------------------------
function test_bireduce(testCase)

% Consider the system:
%
%              dx/dt  = 1./(1 + exp(-x*P)) + u
%                y    = x + e
%

% Specify a model structure
M.f  = inline('1./(1 + exp(-P*x)) + [u; 0]','x','u','P','M');
M.g  = inline('x','x','u','P','M');
M.pE = [-1 .3;.5 -1];           % Prior expectation of parameters
M.pC = speye(4,4);              % Prior covariance for parameters
M.x  = zeros(2,1);               % initial state x(0)
M.m  = 1;                       % number of inputs
M.n  = 2;                       % number of states
M.l  = 2;                       % number of outputs

%y = spm_int(M.pE,M,zeros(50,1));
[M0,M1,L1,L2] = spm_bireduce(M,M.pE);

