classdef test_regress_spm_get_anatomical_layout < matlab.unittest.TestCase
% regression tests for spm_get_anatomical_layout
%__________________________________________________________________________

% Copyright (C) 2026 Department of Imaging Neuroscience



methods (Test)


function test_spm_get_anatomical_layout_1(testCase)
% Layout structure for points on the canonical scalp

spm('defaults','eeg');

[head, fid] = canonical_head();
sens = head(1:200:end,:) + 0.1;
labels = arrayfun(@num2str, (1:size(sens,1))', 'UniformOutput', false);
lay = spm_get_anatomical_layout(sens, labels, head, fid, false, 'mm');

exp = [size(sens,1), 2];
act = size(lay.pos);
testCase.verifyTrue(isequal(exp, act));

exp = labels;
act = lay.label(:);
testCase.verifyTrue(isequal(exp, act));

testCase.verifyTrue(iscell(lay.outline) && iscell(lay.mask));
testCase.verifyTrue(all(isfinite(lay.pos(:))));
end


function test_spm_get_anatomical_layout_2(testCase)
% Fiducials land on the head outline (radius 0.3125, through the
% fiducials) at the front, left and right

spm('defaults','eeg');

[head, fid] = canonical_head();
lay = spm_get_anatomical_layout([fid.NAS; fid.LPA; fid.RPA], {'NAS'; 'LPA'; 'RPA'}, ...
    head, fid, false, 'mm', [], inf);

act = vecnorm(lay.pos, 2, 2);
testCase.verifyLessThan(abs(act - 0.3125), 0.03);

% NAS to the front
testCase.verifyGreaterThan(lay.pos(1,2), 0);
% LPA to the left
testCase.verifyLessThan(lay.pos(2,1), 0);
% RPA to the right
testCase.verifyGreaterThan(lay.pos(3,1), 0);
end


function test_spm_get_anatomical_layout_3(testCase)
% Head points far from the sensors are removed by default, but kept with
% max_distance = inf, and the layout of points near the top is the same

spm('defaults','eeg');

[head, fid] = canonical_head();
[~, top] = maxk(head(:,3), 5);
sens = head(top,:);
labels = arrayfun(@num2str, (1:size(sens,1))', 'UniformOutput', false);
lay_default = spm_get_anatomical_layout(sens, labels, head, fid, false, 'mm');
lay_all = spm_get_anatomical_layout(sens, labels, head, fid, false, 'mm', [], inf);

% Near vertex
testCase.verifyLessThan(max(vecnorm(lay_all.pos, 2, 2)), 0.1);
testCase.verifyLessThan(max(abs(lay_default.pos(:) - lay_all.pos(:))), 0.02);
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
