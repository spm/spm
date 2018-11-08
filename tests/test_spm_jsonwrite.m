function tests = test_spm_jsonwrite
% Unit Tests for spm_jsonwrite
%__________________________________________________________________________
% Copyright (C) 2016-2017 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_jsonwrite.m 7478 2018-11-08 14:51:54Z guillaume $

tests = functiontests(localfunctions);


function test_jsonwrite_array(testCase)
exp = {'one';'two';'three'};
act = spm_jsonread(spm_jsonwrite(exp));
testCase.verifyTrue(isequal(exp, act));

exp = 2;
act = nnz(spm_jsonwrite(1:3) == ',');
testCase.verifyTrue(isequal(exp, act));

function test_jsonwrite_object(testCase)
exp = struct('Width',800,'Height',600,'Title','View from the 15th Floor','Animated',false,'IDs',[116;943;234;38793]);
act = spm_jsonread(spm_jsonwrite(exp));
testCase.verifyTrue(isequal(exp, act));

function test_jsonwrite_all_types(testCase)
exp = [];
act = spm_jsonread(spm_jsonwrite(exp));
testCase.verifyTrue(isequal(exp, act));

exp = [true;false];
act = spm_jsonread(spm_jsonwrite(exp));
testCase.verifyTrue(isequal(exp, act));

exp = struct('a','');
act = spm_jsonread(spm_jsonwrite(exp));
testCase.verifyTrue(isequal(exp, act));

str = struct('str',reshape(1:9,3,3));
exp = spm_jsonread('{"str":[[1,4,7],[2,5,8],[3,6,9]]}');
act = spm_jsonread(spm_jsonwrite(str));
testCase.verifyTrue(isequal(act, exp));

str = [1,2,NaN,3,Inf];
exp = spm_jsonread('[1,2,null,3,null]');
act = spm_jsonread(spm_jsonwrite(str));
testCase.verifyTrue(isequaln(act, exp));
