function tests = test_spm_cat_struct
% Unit Tests for spm_cat_struct
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_cat_struct.m 7074 2017-05-05 10:51:08Z guillaume $

tests = functiontests(localfunctions);


function test_spm_cat_struct_1(testCase)
s1 = struct('a',1,'b',2);
s2 = struct('a',3,'b',4);

exp = struct('a',{1,3},'b',{2,4});
act = spm_cat_struct(s1,s2);
testCase.verifyEqual(act, exp);

function test_spm_cat_struct_2(testCase)
s1 = struct('a',1);
s2 = struct('b',2);

exp = struct('a',{1,[]},'b',{[],2});
act = spm_cat_struct(s1,s2);
testCase.verifyEqual(act, exp);

function test_spm_cat_struct_3(testCase)
s1 = struct('a',{1,2});
s2 = struct('b',3);

exp = struct('a',{1,2,[]},'b',{[],[],3});
act = spm_cat_struct(s1,s2);
testCase.verifyEqual(act, exp);

exp = struct('a',{[],1,2},'b',{3,[],[]});
act = spm_cat_struct(s2,s1);
testCase.verifyEqual(act, exp);

function test_spm_cat_struct_4(testCase)
s1 = struct('a',{1;2});
s2 = struct('b',3);

exp = struct('a',{1;2;[]},'b',{[];[];3});
act = spm_cat_struct(s1,s2);
testCase.verifyEqual(act, exp);

exp = struct('a',{[],1,2},'b',{3,[],[]});
act = spm_cat_struct(s2,s1);
testCase.verifyEqual(act, exp);

function test_spm_cat_struct_5(testCase)
s1 = struct('a',{1,2});
s2 = struct('b',3);
s3 = struct('a',4,'c',5);

exp = struct('a',{1,2,[],4},'b',{[],[],3,[]},'c',{[],[],[],5});
act = spm_cat_struct(s1,s2,s3);
testCase.verifyEqual(act, exp);

function test_spm_cat_struct_6(testCase)
s1 = struct('a',{1,2;3 4},'b',{5,6;7,8});
s2 = struct('b',{0,1},'c',{3,2});

exp = struct('a',{1,3,2,4,[],[]},'b',{5,7,6,8,0,1},'c',{[],[],[],[],3,2})';
act = spm_cat_struct(s1,s2);
testCase.verifyEqual(act, exp);
