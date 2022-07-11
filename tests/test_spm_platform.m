function tests = test_spm_platform
% Unit Tests for spm_platform
%__________________________________________________________________________

% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_platform_init(testCase)
import matlab.unittest.constraints.*

spm_platform;
p = spm_platform;
testCase.verifyThat(p, IsOfClass('struct'));
p = spm_platform('init');
testCase.verifyThat(p, IsOfClass('struct'));

function test_spm_platform_query(testCase)
import matlab.unittest.constraints.*

val = spm_platform('bigend');
testCase.verifyThat(val, IsOfClass('logical'));
val = spm_platform('user');
testCase.verifyThat(val, IsOfClass('char'));
val = spm_platform('host');
testCase.verifyThat(val, IsOfClass('char'));
val = spm_platform('tempdir');
testCase.verifyThat(val, IsOfClass('char'));
val = spm_platform('desktop');
testCase.verifyThat(val, IsOfClass('logical'));

function test_spm_platform_fonts(testCase)
import matlab.unittest.constraints.*

fonts = spm_platform('fonts');
testCase.verifyThat(fonts, IsOfClass('struct'));
font = spm_platform('font','helvetica');
testCase.verifyThat(font, IsOfClass('char'));
font = spm_platform('font','times');
testCase.verifyThat(font, IsOfClass('char'));
font = spm_platform('font','courier');
testCase.verifyThat(font, IsOfClass('char'));
font = spm_platform('font','symbol');
testCase.verifyThat(font, IsOfClass('char'));

function test_spm_platform_memory(testCase)
import matlab.unittest.constraints.*

meminfo = spm_platform('memory');
testCase.verifyThat(meminfo, IsOfClass('struct'));
meminfo = spm_platform('memory','available');
testCase.verifyThat(meminfo, IsOfClass('double'));
meminfo = spm_platform('memory','total');
testCase.verifyThat(meminfo, IsOfClass('double'));
