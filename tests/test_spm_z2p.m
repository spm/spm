% Unit Tests for spm_z2p
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_z2p.m 6352 2015-02-27 18:30:35Z guillaume $

%% test 1: Gaussian
exp = 0.05;
act = spm_z2p(1.644853626951472,[],'Z');
tol = 1e-8;
assert(abs(act - exp) < tol);

%% test 2: T
exp = 0.05;
act = spm_z2p(1.745883689098006,[1 16],'T');
tol = 1e-8;
assert(abs(act - exp) < tol);

%% test 3: Chi squared
exp = 0.05;
act = spm_z2p(26.296227607876062,[1 16],'X');
tol = 1e-8;
assert(abs(act - exp) < tol);

%% test 4: F
exp = 0.05;
act = spm_z2p(3.633723519202212,[2 16],'F');
tol = 1e-8;
assert(abs(act - exp) < tol);

%% test 5: P
exp = 0.5;
act = spm_z2p(0.5,[],'P');
tol = 1e-8;
assert(abs(act - exp) < tol);
