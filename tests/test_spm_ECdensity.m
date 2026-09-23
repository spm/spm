classdef test_spm_ECdensity < matlab.unittest.TestCase
% Unit Tests for spm_ECdensity
%
% The EC densities are checked three ways:
%  * against closed-form reference values (Worsley et al. 1996),
%  * chi^2_1 == Z^2, so its EC densities must equal twice the Gaussian
%    densities at sqrt(t),
%  * chi^2_v == v*F(v,inf), so its EC densities must match the F-field
%    densities at t/v in the large denominator-df limit.
%__________________________________________________________________________

% Copyright (C) 2026 Wellcome Centre for Human Neuroimaging


methods (Test)


function test_chi2_reference_values(testCase)
% closed-form chi^2 EC densities: rho_d(t) for d = 1,2,3
tol = 1e-12;

EC  = spm_ECdensity('X',15,[1 10]);
testCase.verifyEqual(EC(2), 1.875960195166651e-01,'RelTol',tol);
testCase.verifyEqual(EC(3), 1.930554347139083e-01,'RelTol',tol);
testCase.verifyEqual(EC(4), 6.622457656369271e-02,'RelTol',tol);

EC  = spm_ECdensity('X',8,[1 3]);
testCase.verifyEqual(EC(2), 7.766134591989517e-02,'RelTol',tol);
testCase.verifyEqual(EC(3), 1.094369452067079e-01,'RelTol',tol);
testCase.verifyEqual(EC(4), 1.113765748270419e-01,'RelTol',tol);
end

function test_chi2_df1_matches_gaussian(testCase)
% chi^2_1 == Z^2: excursion {chi^2_1 >= t} == {|Z| >= sqrt(t)}, so
% every chi^2_1 EC density equals twice the Gaussian one at sqrt(t)
% (the d=0 row compares gammainc- against erf-based tail probabilities,
% which agree less tightly than the closed-form density rows)
t   = [8 15 25];

ECX = spm_ECdensity('X',t,[1 1]);
ECZ = spm_ECdensity('Z',sqrt(t),[]);
testCase.verifyEqual(ECX(1,:), 2*ECZ(1,:),'RelTol',1e-8);
for d = 2:4
    testCase.verifyEqual(ECX(d,:), 2*ECZ(d,:),'RelTol',1e-10);
end
end

function test_chi2_matches_F_limit(testCase)
% chi^2_v == v*F(v,inf): EC densities must match the F-field densities
% at t/v as the denominator df grows (reference error ~O(1/V))
tol = 1e-5;
t   = [8 15 25];
V   = 1e8;

for v = [3 6 10]
    ECX = spm_ECdensity('X',t,[1 v]);
    ECF = spm_ECdensity('F',t/v,[v V]);
    for d = 2:4
        testCase.verifyEqual(ECX(d,:), ECF(d,:),'RelTol',tol);
    end
end
end

function test_t_matches_gaussian_limit(testCase)
% T-field densities converge to the Gaussian ones as df -> inf
tol = 1e-5;
t   = [2 3 4 5];

ECT = spm_ECdensity('T',t,[1 1e8]);
ECZ = spm_ECdensity('Z',t,[]);
for d = 1:4
    testCase.verifyEqual(ECT(d,:), ECZ(d,:),'RelTol',tol);
end
end


end
end
