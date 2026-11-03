classdef test_regress_spm_opm_simulate_array < matlab.unittest.TestCase
% regresion tests for spm_opm_simulate_array
%__________________________________________________________________________

% Copyright (C) 2026 Department of Imaging Neuroscience



methods (Test)


function test_spm_opm_simulate_array_1(testCase)
% Symmetric equidistant layout without a head: right side mirrors left,
% everything inside the default mask, no positions

spm('defaults','eeg');

S = [];
S.n_sensors = 8;
S.symmetry = 'left';
S.midline = true;
S.plot = false;
[positions, lay] = spm_opm_simulate_array(S);

testCase.verifyTrue(isempty(positions));

exp = [S.n_sensors, 2];
act = size(lay.pos);
testCase.verifyTrue(isequal(exp, act));

exp = 8;
act = numel(unique(lay.label));
testCase.verifyTrue(isequal(exp, act));

left = find(startsWith(lay.label, 'L'));
[~, right] = ismember(strrep(lay.label(left), 'L', 'R'), lay.label);
act = max(abs(lay.pos(right,:) - [-lay.pos(left,1), lay.pos(left,2)]), [], 'all');
testCase.verifyLessThan(act, 1e-9);
testCase.verifyTrue(all(lay.pos(startsWith(lay.label, 'M'),1) == 0));

mask = lay.mask{1};
testCase.verifyTrue(all(inpolygon(lay.pos(:,1), lay.pos(:,2), mask(:,1), mask(:,2))));
end


function test_spm_opm_simulate_array_2(testCase)
% A symmetric layout off the midline needs an even number of sensors

spm('defaults','eeg');

S = [];
S.n_sensors = 9;
S.midline = false;
S.plot = false;
failed = false;
try
    spm_opm_simulate_array(S);
catch
    failed = true;
end
testCase.verifyTrue(failed);
end


function test_spm_opm_simulate_array_3(testCase)
% Triaxial array on the canonical head: three channels per sensor,
% orthonormal axes, sensors offset from the scalp, positions written for
% spm_opm_sim

spm('defaults','eeg');

out_file = [tempname, '.tsv'];
S = [];
S.mesh = spm_eeg_inv_mesh([], 2);
S.n_sensors = 16;
S.axes = 'triaxial';
S.offset = 11;
S.output_file = out_file;
S.plot = false;
[positions, lay, info] = spm_opm_simulate_array(S);

exp = 3 * size(lay.pos,1);
act = height(positions);
testCase.verifyTrue(isequal(exp, act));

testCase.verifyTrue(all(isfield(info.fid, {'NAS', 'LPA', 'RPA'})));

O = [positions.Ox, positions.Oy, positions.Oz];
for k = 1:size(lay.pos,1)
    axes_k = O(3*k-2:3*k,:);
    testCase.verifyLessThan(max(abs(axes_k * axes_k' - eye(3)), [], 'all'), 1e-3);
end

act = vecnorm(info.pos - info.scalp_pos, 2, 2);
testCase.verifyLessThan(max(abs(act - 11)), 1e-6);

testCase.verifyTrue(exist(out_file, 'file') == 2);
written = spm_load(out_file);
exp = positions.name;
act = written.name;
testCase.verifyTrue(isequal(exp, act));
end


function test_spm_opm_simulate_array_4(testCase)
% A layout put back on the head it was made on gives the same positions,
% with the axes chosen

spm('defaults','eeg');

mesh = spm_eeg_inv_mesh([], 1);
S = [];
S.mesh = mesh;
S.n_sensors = 8;
S.plot = false;
[positions, lay] = spm_opm_simulate_array(S);

S = [];
S.mesh = mesh;
S.lay = lay;
S.axes = 2;
S.plot = false;
positions_dual = spm_opm_simulate_array(S);

exp = 2 * height(positions);
act = height(positions_dual);
testCase.verifyTrue(isequal(exp, act));

act = max(abs([positions.Px - positions_dual.Px(1:2:end); ...
    positions.Py - positions_dual.Py(1:2:end); ...
    positions.Pz - positions_dual.Pz(1:2:end)]));
testCase.verifyLessThan(act, 1e-6);
end
end % methods (Test)

end % classdef
