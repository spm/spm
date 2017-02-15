function tests = test_spm_Ce
% Unit Tests for spm_Ce
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_Ce.m 7018 2017-02-15 13:36:48Z guillaume $

tests = functiontests(localfunctions);


function test_spm_Ce_1(testCase)
import matlab.unittest.constraints.*
N = [16 32];
C = spm_Ce(N);
testCase.verifyThat(C, IsOfClass('cell'));
testCase.verifyThat(C, HasLength(numel(N)));
for i=1:numel(N)
    testCase.verifyThat(C{i}, HasSize([sum(N) sum(N)]));
end

C = spm_Ce(N,[]);
testCase.verifyThat(C, IsEqualTo(spm_Ce(N)));

C = spm_Ce(N,0.2);
testCase.verifyThat(C, HasLength(2*numel(N)));


function test_spm_Ce_2(testCase)
import matlab.unittest.constraints.*
N = [16 32 64];
C = spm_Ce('ar',N);
testCase.verifyThat(C, IsEqualTo(spm_Ce(N)));

C = spm_Ce('ar',N,[]);
testCase.verifyThat(C, IsEqualTo(spm_Ce(N,[])));

C = spm_Ce('ar',N,0.4);
testCase.verifyThat(C, IsEqualTo(spm_Ce(N,0.4)));


function test_spm_Ce_3(testCase)
import matlab.unittest.constraints.*
C = spm_Ce('fast',16,3);
testCase.verifyThat(C, IsOfClass('cell'));

C = spm_Ce('fast',[16 32],3);
testCase.verifyThat(C, IsOfClass('cell'));
