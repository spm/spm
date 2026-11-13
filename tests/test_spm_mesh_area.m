classdef test_spm_mesh_area < matlab.unittest.TestCase
    % Unit Tests for spm_mesh_area
    %__________________________________________________________________________

    % Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging



    methods (Test)

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
        end

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
        end

        function test_spm_mesh_area_inputs(testCase)
            meshname = fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii');
            M = gifti(meshname);
            M = export(M,'patch');

            Nv = size(M.vertices,1);
            Nf = size(M.faces,1);

            A = spm_mesh_area(M,'face');
            testCase.verifyEqual(size(A), [1 Nf]);
            testCase.verifyEqual(round(median(A)), 14);

            B = spm_mesh_area(M,'vertex');
            testCase.verifyEqual(size(B), [Nv 1]);
            testCase.verifyEqual(round(median(B)), 29);

            C = spm_mesh_area(M,'sum');

            testCase.verifyEqual(size(C), [1 1]);
            testCase.verifyEqual(round(C), 151070);

            D = spm_mesh_area(M, false);
            testCase.verifyEqual(D, C);

            E = spm_mesh_area(M, true);
            testCase.verifyEqual(E, A);
        end

    end % methods (Test)

end % classdef