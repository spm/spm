function tests = test_spm_mesh_borders
% Unit Tests for spm_mesh_borders
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_mesh_borders.m 7413 2018-09-07 09:53:34Z guillaume $

tests = functiontests(localfunctions);


function test_spm_mesh_borders_sphere(testCase)
import matlab.unittest.constraints.*
M = spm_mesh_polyhedron('icosahedron');
[B,C] = spm_mesh_borders(M);

testCase.verifyThat(B, IsEmpty);
testCase.verifyThat(C, IsEmpty);

exp = sort(M.faces(end,:))';
M.faces(end,:) = [];
[act,C] = spm_mesh_borders(M);
testCase.verifyEqual(act, exp);
testCase.verifyEqual(numel(C), 1);

function test_spm_mesh_borders_sphere_4(testCase)
M = spm_mesh_sphere(4);
T = M.vertices(:,1) < 0;
M = spm_mesh_split(M,T);
[B,C] = spm_mesh_borders(M);

% figure, plot(gifti(M)); hold on
% plot3(M.vertices(B,1),M.vertices(B,2),M.vertices(B,3),'*r');

exp = 94;
act = numel(B);
testCase.verifyEqual(act, exp);
testCase.verifyEqual(numel(C), 1);

function test_spm_mesh_borders_sphere_5(testCase)
M = spm_mesh_sphere(5);
T = M.vertices(:,1) < 0;
M = spm_mesh_split(M,T);
[B,C] = spm_mesh_borders(M);

exp = 190;
act = numel(B);
testCase.verifyEqual(act, exp);
testCase.verifyEqual(numel(C), 1);

function test_spm_mesh_borders_two_spheres(testCase)
M1 = spm_mesh_sphere(4);
M2 = spm_mesh_sphere(5);
M2.vertices(:,2) = M2.vertices(:,2) + 3;
M = spm_mesh_join([M1,M2]);

T = M.vertices(:,1) < 0;
M = spm_mesh_split(M,T);
[B,C] = spm_mesh_borders(M);

exp = 190 + 94;
act = numel(B);
testCase.verifyEqual(act, exp);
testCase.verifyEqual(numel(C), 2);
testCase.verifyEqual(sort(cellfun(@numel,C)), [94 190]);
