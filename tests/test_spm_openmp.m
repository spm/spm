function tests = test_spm_openmp
% Unit Tests for OpenMP
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_openmp.m 7690 2019-11-07 17:02:55Z guillaume $

tests = functiontests(localfunctions);
tests(1:end) = []; % ALL TESTS ARE DISABLED

% Push returns different results because of the order of the summation.
% No order is more right than the other, so it shouldn't be a problem.
% We could use something like a -DSPM_DETERMINISTIC mode under which these
% kind of sections are not parallelised.

% There is a slight issue in the checkerboard scheme with circulant
% boundary conditions when one of the dimensions is not a multiple of 3.
% This makes fmg return slighlty different results in mono vs multi
% threads. Setting up a tolerance for this step would prevent the test to
% fail.

function test_spm_openmp_diffeo_bsplinc(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('bsplinc\n');
fprintf('----------\n');
dim     = [100 100 100];
i       = randn(dim, 'single');
fprintf('Number of threads: %d\n', spm_set_num_threads(1));
tic
o = spm_diffeo('bsplinc', i, [7 7 7]);
t1 = toc;
fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
tic
oo = spm_diffeo('bsplinc', i, [7 7 7]);
t2 = toc;
fprintf('Speedup: %g\n', t1/t2);
fprintf('Same output: %d\n', isequaln(o,oo));
testCase.verifyEqual(o, oo);


function test_spm_openmp_diffeo_bsplins(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('bsplins\n');
fprintf('----------\n');
dimi      = [100 100 100];
dimo      = 2*dimi;
i         = randn(dimi, 'single');
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dimo(1)),single(1:dimo(2)),single(1:dimo(3)));
id        = cat(4,id{:});
y         = (id -1)*0.5 + 1;
fprintf('Number of threads: %d\n', spm_set_num_threads(1));
tic
o = spm_diffeo('bsplins', i, y, [7 7 7]);
t1 = toc;
fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
tic
oo = spm_diffeo('bsplins', i, y, [7 7 7]);
t2 = toc;
fprintf('Speedup: %g\n', t1/t2);
fprintf('Same output: %d\n', isequaln(o,oo));
testCase.verifyEqual(o, oo);


function test_spm_openmp_diffeo_pull(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('pull\n');
fprintf('----------\n');
dimo      = [200 200 200];
dimi      = [100 100 100];
i         = randn(dimi, 'single');
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dimo(1)),single(1:dimo(2)),single(1:dimo(3)));
id        = cat(4,id{:});
y         = (id -1)*0.5 + 1;
fprintf('Number of threads: %d\n', spm_set_num_threads(1));
tic
o = spm_diffeo('pullc', i, y);
t1 = toc;
fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
tic
oo = spm_diffeo('pullc', i, y);
t2 = toc;
fprintf('Speedup: %g\n', t1/t2);
fprintf('Same output: %d\n', isequaln(o,oo));
testCase.verifyEqual(o, oo);


function test_spm_openmp_diffeo_push(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('push\n');
fprintf('----------\n');
dimo      = [200 200 200];
dimi      = [100 100 100];
i         = randn(dimi, 'single');
id        = cell(3,1);
[id{1:3}] = ndgrid(single(1:dimi(1)),single(1:dimi(2)),single(1:dimi(3)));
id        = cat(4,id{:});
y         = (id -1)*0.5 + 1;
fprintf('Number of threads: %d\n', spm_set_num_threads(1));
tic
[o,c] = spm_diffeo('pushc', i, y, dimo);
t1 = toc;
fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
tic
[oo,cc] = spm_diffeo('pushc', i, y, dimo);
t2 = toc;
fprintf('Speedup: %g\n', t1/t2);
fprintf('Same output: %d\n', isequaln(o,oo) && isequal(c,cc));
tol = single(1E-5);
testCase.verifyEqual(o, oo,'AbsTol',tol);
testCase.verifyEqual(c, cc,'AbsTol',tol);


function test_spm_openmp_field_vel2mom(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('vel2mom (field)\n');
fprintf('----------\n');
dim = [200 200 200];
i   = ones([dim 3], 'single');
for bnd=[0 1]
    spm_field('boundary', bnd);
    fprintf('* boundary: %d\n', bnd);
    fprintf('Number of threads: %d\n', spm_set_num_threads(1));
    tic
    o = spm_field('vel2mom',i,[1 1 1 10 100 10000]);
    t1 = toc;
    fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
    tic
    oo = spm_field('vel2mom',i,[1 1 1 10 100 10000]);
    t2 = toc;
    fprintf('Speedup: %g\n', t1/t2);
    fprintf('Same output: %d\n', isequaln(o,oo));
    testCase.verifyEqual(o, oo);
end


function test_spm_openmp_diffeo_vel2mom(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('vel2mom (diffeo)\n');
fprintf('----------\n');
dim = [200 200 200];
i   = ones([dim 3], 'single');
for bnd=[0 1]
    spm_diffeo('boundary', bnd);
    fprintf('* boundary: %d\n', bnd);
    fprintf('Number of threads: %d\n', spm_set_num_threads(1));
    tic
    o = spm_diffeo('vel2mom',i,[1 1 1 10 100 10000 10 10]);
    t1 = toc;
    fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
    tic
    oo = spm_diffeo('vel2mom',i,[1 1 1 10 100 10000 10 10]);
    t2 = toc;
    fprintf('Speedup: %g\n', t1/t2);
    fprintf('Same output: %d\n', isequaln(o,oo));
    testCase.verifyEqual(o, oo);
end


function test_spm_openmp_field_fmg(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('fmg (field)\n');
fprintf('----------\n');
dim = [200 200 200];
H   = cat(4, ones([dim 3], 'single'), 1E-5*ones([dim 3], 'single'));
g   = randn([dim 3], 'single');
for bnd=[0 1]
    spm_field('boundary', bnd);
    fprintf('* boundary: %d\n', bnd);
    fprintf('Number of threads: %d\n', spm_set_num_threads(1));
    tic
    o = spm_field(H,g,[1 1 1 10 100 10000 2 2]);
    t1 = toc;
    fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
    tic
    oo = spm_field(H,g,[1 1 1 10 100 10000 2 2]);
    t2 = toc;
    fprintf('Speedup: %g\n', t1/t2);
    fprintf('Same output: %d\n', isequaln(o,oo));
    testCase.verifyEqual(o, oo);
end


function test_spm_openmp_diffeo_fmg(testCase)
fprintf('\n');
fprintf('----------\n');
fprintf('fmg (diffeo)\n');
fprintf('----------\n');
dim = [200 200 200];
H   = cat(4, ones([dim 3], 'single'), 1E-5*ones([dim 3], 'single'));
g   = randn([dim 3], 'single');
for bnd=[0 1]
    spm_diffeo('boundary', bnd);
    fprintf('* boundary: %d\n', bnd);
    fprintf('Number of threads: %d\n', spm_set_num_threads(1));
    tic
    o = spm_diffeo('fmg',H,g,[1 1 1 10 100 10000 10 10 2 2]);
    t1 = toc;
    fprintf('Number of threads: %d\n', spm_set_num_threads(-1));
    tic
    oo = spm_diffeo('fmg',H,g,[1 1 1 10 100 10000 10 10 2 2]);
    t2 = toc;
    fprintf('Speedup: %g\n', t1/t2);
    fprintf('Same output: %d\n', isequaln(o,oo));
    testCase.verifyEqual(o, oo);
end


function n = spm_set_num_threads(n)
setenv('SPM_NUM_THREADS',sprintf('%d',n));
try, spm_diffeo; end
n = sscanf(getenv('SPM_NUM_THREADS'),'%d');
