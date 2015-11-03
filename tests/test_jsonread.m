function tests = test_jsonread
% Unit Tests for jsonread
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_jsonread.m 6589 2015-11-03 16:01:08Z guillaume $

tests = functiontests(localfunctions);


function test_jsonread_from_string(testCase)
exp = struct('name','Karl','age',56);
act = spm_jsonread('{ "name" : "Karl", "age" : 56 }');
testCase.verifyTrue(isequal(exp, act));
