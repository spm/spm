function tests = test_spm_mesh_area
% Unit Tests for spm_mesh_area
%__________________________________________________________________________

% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_area_polyhedron(testCase)

M = spm_mesh_polyhedron('icosahedron');
a = 2;

exp = 5*sqrt(3)*a^2;
act = spm_mesh_area(M);
testCase.verifyEqual(act, exp, 'AbsTol',1e-6);

M = spm_mesh_polyhedron('octahedron');
a = sqrt(2);

exp = 2*sqrt(3)*a^2;
act = spm_mesh_area(M);
testCase.verifyEqual(act, exp, 'AbsTol',1e-6);

M = spm_mesh_polyhedron('tetrahedron');
a = 2;

exp = sqrt(3)*a^2;
act = spm_mesh_area(M);
testCase.verifyEqual(act, exp, 'AbsTol',1e-6);


function test_spm_mesh_area_sphere(testCase)

M = spm_mesh_sphere;

exp = 4*pi;
act = spm_mesh_area(M);
testCase.verifyEqual(act, exp, 'AbsTol',1e-2);

exp = 4*pi;
act = sum(spm_mesh_area(M,'face'));
testCase.verifyEqual(act, exp, 'AbsTol',1e-2);

exp = 4*pi;
act = sum(spm_mesh_area(M,'vertex'));
testCase.verifyEqual(act, exp, 'AbsTol',1e-2);
