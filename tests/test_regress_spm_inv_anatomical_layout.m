classdef test_regress_spm_inv_anatomical_layout < matlab.unittest.TestCase
% Unit Tests for spm_inv_anatomical_layout
%__________________________________________________________________________

% Copyright (C) 2026 Department of Imaging Neuroscience



methods (Test)


function test_spm_inv_anatomical_layout_1(testCase)
% Three rows per position (radial, primary and secondary tangential),
% orthonormal axes, radial pointing into the head

spm('defaults','eeg');

[head, fid] = canonical_head();
theta = linspace(0, 2*pi, 9)';
pos2d = [0 0; 0.15 * [cos(theta(1:end-1)), sin(theta(1:end-1))]];
T = spm_inv_anatomical_layout(pos2d, head, fid, false, 'mm');

exp = 3 * size(pos2d,1);
act = height(T);
testCase.verifyTrue(isequal(exp, act));

P = [T.Px, T.Py, T.Pz];
O = [T.Ox, T.Oy, T.Oz];
testCase.verifyLessThan(max(abs(vecnorm(O, 2, 2) - 1)), 1e-6);
for k = 1:size(pos2d,1)
    axes_k = O(3*k-2:3*k,:);
    testCase.verifyLessThan(max(abs(axes_k * axes_k' - eye(3)), [], 'all'), 1e-3);
end

% Radial into the head, towards the middle of the scalp
radial = O(1:3:end,:);
inward = mean(head) - P(1:3:end,:);
testCase.verifyTrue(all(sum(radial .* inward, 2) > 0));
end


function test_spm_inv_anatomical_layout_2(testCase)
% The precalculated grid gives the same positions

spm('defaults','eeg');

[head, fid] = canonical_head();
[T, precalculated] = spm_inv_anatomical_layout([0 0; 0.1 0.1], head, fid, false, 'mm');
T2 = spm_inv_anatomical_layout([0 0; 0.1 0.1], head, fid, false, 'mm', precalculated);

cz = [T.Px(1), T.Py(1), T.Pz(1)];
testCase.verifyLessThan(max(head(:,3)) - cz(3), 15);

act = max(abs([T.Px - T2.Px; T.Py - T2.Py; T.Pz - T2.Pz]));
testCase.verifyLessThan(act, 1e-9);
end


function test_spm_inv_anatomical_layout_3(testCase)
% Round trip: 2Dto 3D with spm_inv_anatomical_layout to 2D with
% spm_get_anatomical_layout gives back approximately the layout positions.

spm('defaults','eeg');

[head, fid] = canonical_head();
[gx, gy] = meshgrid(-0.25:0.05:0.25);
pos2d = [gx(:), gy(:)];
pos2d = pos2d(vecnorm(pos2d, 2, 2) <= 0.25,:);
T = spm_inv_anatomical_layout(pos2d, head, fid, false, 'mm');
P = [T.Px(1:3:end), T.Py(1:3:end), T.Pz(1:3:end)];
labels = arrayfun(@num2str, (1:size(P,1))', 'UniformOutput', false);
lay = spm_get_anatomical_layout(P, labels, head, fid, false, 'mm', [], inf);

act = max(vecnorm(lay.pos - pos2d, 2, 2));
testCase.verifyLessThan(act, 0.01);
end
end % methods (Test)

end % classdef


function [head, fid] = canonical_head()
% Canonical scalp vertices and fiducials (mm)
mesh = spm_eeg_inv_mesh([], 2);
scalp = gifti(mesh.tess_scalp);
head = double(scalp.vertices);
labels = upper(mesh.fid.fid.label);
for k = 1:numel(labels)
    if any(strcmp(labels{k}, {'NAS', 'LPA', 'RPA'}))
        fid.(labels{k}) = mesh.fid.fid.pnt(k,:);
    end
end
end
