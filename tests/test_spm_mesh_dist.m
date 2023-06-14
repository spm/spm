function tests = test_spm_mesh_dist
% Unit Tests for spm_mesh_dist
%__________________________________________________________________________

% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);


function test_spm_mesh_dist_(testCase)
M.vertices = [...
     1    -1    -1
     1     0    -1
     1     1    -1
     1    -1     0
     1     0     0
     1     1     0
     1    -1     1
     1     0     1
     1     1     1
    -1    -1    -1
    -1     0    -1
    -1     1    -1
    -1    -1     0
    -1     0     0
    -1     1     0
    -1    -1     1
    -1     0     1
    -1     1     1
     0     1    -1
     0     1     0
     0     1     1
     0    -1    -1
     0    -1     0
     0    -1     1
     0     0     1
     0     0    -1];

M.faces = [...
     1     2     4
     2     5     4
     2     3     5
     3     6     5
     4     5     7
     5     8     7
     5     6     8
     6     9     8
    13    11    10
    13    14    11
    14    12    11
    14    15    12
    16    14    13
    16    17    14
    17    15    14
    17    18    15
    15    19    12
    15    20    19
    20     3    19
    20     6     3
    18    20    15
    18    21    20
    21     6    20
    21     9     6
    10    22    13
    22    23    13
    22     1    23
     1     4    23
    13    23    16
    23    24    16
    23     4    24
     4     7    24
    16    24    17
    24    25    17
    24     7    25
     7     8    25
    17    25    18
    25    21    18
    25     8    21
     8     9    21
    11    22    10
    11    26    22
    26     1    22
    26     2     1
    12    26    11
    12    19    26
    19     2    26
    19     3     2];

act = zeros(1,numel(-2:0.13:2)^3);
exp = act;
i = 1;
for x=-2:0.13:2
    for y=-2:0.13:2
        for z=-2:0.13:2
            act(i) = spm_mesh_dist(M,[x y z]);
            exp(i) = dist_cube([x y z]);
            i = i + 1;
        end
    end
end
testCase.verifyEqual(act, exp, 'AbsTol',1e-10);

[X,Y,Z] = meshgrid(-2:0.13:2,-2:0.13:2,-2:0.13:2);
XYZ = [X(:) Y(:) Z(:)];
act = spm_mesh_dist(M,XYZ);
exp = zeros(size(act));
for i=1:size(XYZ,1)
    exp(i) = dist_cube(XYZ(i,:));
end
testCase.verifyEqual(act, exp, 'AbsTol',1e-10);


function d = dist_cube(XYZ)
dx = max([[-1 -1 -1] - XYZ; XYZ - [1 1 1]], [], 1);
if max(dx) < 0
    d = max(dx);
else
    d = sqrt(sum(dx(dx>0).^2));
end
