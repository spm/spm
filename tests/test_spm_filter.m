classdef test_spm_filter < matlab.unittest.TestCase
% Unit Tests for spm_filter
%__________________________________________________________________________

% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging



methods (Test)


function test_spm_filter_1(testCase)
import matlab.unittest.constraints.*
K = struct('RT',2.4,'row',1:64,'HParam',128);
K = spm_filter(K);

testCase.verifyThat(K, HasField('X0'));
testCase.verifyThat(K, HasLength(1));
testCase.verifyThat(K.X0, HasSize([64 2]));
end

function test_spm_filter_2(testCase)
K  = struct('RT',2.4,'row',1:64,'HParam',128);
Y  = rand(64,1);
Yf = spm_filter(K,Y);
K  = spm_filter(K);

exp = Yf;
act = spm_filter(eye(64)- K.X0*K.X0',Y);
tol = 1e-10;
testCase.verifyEqual(exp, act,'AbsTol',tol);
end

function test_spm_filter_3(testCase)
import matlab.unittest.constraints.*
K  = struct('RT',2.4,'row',{1:64,65:128},'HParam',128);
Y  = rand(128,2);
K  = spm_filter(K);
Yf = spm_filter(K,Y);

testCase.verifyThat(Yf, HasSize(size(Y)));
end
end % methods (Test)

end % classdef