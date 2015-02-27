function tests = test_spm_filter
% Unit Tests for spm_filter
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_filter.m 6352 2015-02-27 18:30:35Z guillaume $

tests = functiontests(localfunctions);


function test_spm_filter_1(testCase)
import matlab.unittest.constraints.*
K = struct('RT',2.4,'row',1:64,'HParam',128);
K = spm_filter(K);

testCase.verifyThat(K, HasField('X0'));
testCase.verifyThat(K, HasLength(1));
testCase.verifyThat(K.X0, HasSize([64 2]));


function test_spm_filter_2(testCase)
K  = struct('RT',2.4,'row',1:64,'HParam',128);
Y  = rand(64,1);
Yf = spm_filter(K,Y);
K  = spm_filter(K);

exp = Yf;
act = spm_filter(eye(64)- K.X0*K.X0',Y);
tol = 1e-10;
testCase.verifyEqual(exp, act,'AbsTol',tol);


function test_spm_filter_3(testCase)
import matlab.unittest.constraints.*
K  = struct('RT',2.4,'row',{1:64,65:128},'HParam',128);
Y  = rand(128,2);
K  = spm_filter(K);
Yf = spm_filter(K,Y);

testCase.verifyThat(Yf, HasSize(size(Y)));
