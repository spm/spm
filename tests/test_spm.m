function tests = test_spm
% Unit Tests for spm
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% $Id: test_spm.m 8124 2021-07-12 16:26:10Z guillaume $

tests = functiontests(localfunctions);


function test_spm_multi(testCase)
import matlab.unittest.constraints.*

d = spm('Dir');
testCase.verifyThat(d,IsOfClass('char'));

v = spm('Ver');
testCase.verifyThat(v,IsOfClass('char'));

v = spm('Version');
testCase.verifyThat(v,IsOfClass('char'));

xTB = spm('TBs');
testCase.verifyThat(xTB,IsOfClass('struct'));

u = spm('GetUser');
testCase.verifyThat(u,IsOfClass('char'));
u = spm('GetUser','hello %s!');
testCase.verifyThat(u,IsOfClass('char'));

t = spm('Time');
testCase.verifyThat(t,IsOfClass('char'));

mem = spm('Memory');
testCase.verifyThat(mem, IsOfClass('double'));
mem_avail = spm('Memory','available');
testCase.verifyThat(mem_avail, IsOfClass('double'));
mem_total = spm('Memory','total');
testCase.verifyThat(mem_total, IsOfClass('double'));
