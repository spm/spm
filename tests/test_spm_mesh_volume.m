function tests = test_spm_mesh_volume
% Unit Tests for spm_mesh_volume
%__________________________________________________________________________

% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_volume_polyhedron(testCase)

M = spm_mesh_polyhedron('icosahedron');
a = 2;

exp = 5/12*(3+sqrt(5))*a^3;
act = spm_mesh_volume(M);
testCase.verifyEqual(act, exp, 'AbsTol',1e-6);

M = spm_mesh_polyhedron('octahedron');
a = sqrt(2);

exp = sqrt(2)/3*a^3;
act = spm_mesh_volume(M);
testCase.verifyEqual(act, exp, 'AbsTol',1e-6);

% M = spm_mesh_polyhedron('tetrahedron');
% a = 2;
% 
% exp = sqrt(2)/12*a^3;
% act = spm_mesh_volume(M);
% testCase.verifyEqual(act, exp, 'AbsTol',1e-6);


function test_spm_mesh_volume_sphere(testCase)

M = spm_mesh_sphere;
r = 1;

exp = 4/3*pi*r^3;
act = spm_mesh_volume(M);
testCase.verifyEqual(act, exp, 'AbsTol',1e-2);
