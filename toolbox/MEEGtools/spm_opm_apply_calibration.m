function [cal_D] = spm_opm_apply_calibration(S)
% [cal_D] = spm_opm_apply_calibration(S) applies position,
% orientation and gain calibration to an existing dataset. Gain applied is
% relative to the existing (assumed) calibration of the original data. If
% positions and orientations are already available for channels then the
% function will attempt to align the new calibration to the existing space.
%   Input:
%       S : structure with fields
%           - D:         MEEG D object to be calibrated (required)
%           - positions: positions file with gain column
%           - chunkSize: max memory usage (for large datasets; default: 512)
%           - prefix:    string, calibrated object prefix (default: 'cal_')
%           - balance:   boolean, update existing forward model/s (default:
%                            false)
%
% Author: Nicholas Alexander (n.alexander@ucl.ac.uk)
%     Adapted from spm_opm_hfc
% Copyright: Department of Imaging Neuroscience, UCL, 2026
%
% TODO: Add delay correction
%       Improve stability of alignment
%       Set bad channels if uncalibrated MEG

if ~isfield(S, 'D') || isempty(S.D)
    error('Specify D object')
end
if ~isfield(S, 'positions') || isempty(S.positions)
    error('Specify positions table')
end
if ~isfield(S, 'chunkSize') || isempty(S.chunkSize)
    S.chunkSize = 512;
end
if ~isfield(S, 'prefix') || isempty(S.prefix)
    S.prefix = 'cal_';
end
if ~isfield(S, 'balance') || isempty(S.balance)
    S.balance = false;
end

% Check if there is existing position information
grad = S.D.sensors('MEG');
positions_cal = S.positions;

% Attempt to align calibrated positions to existing space
if ~isempty(grad)
    % Positions format
    positions_original = table(grad.label(:), ...
        grad.chanpos(:,1), grad.chanpos(:,2), grad.chanpos(:,3), ...
        grad.chanori(:,1), grad.chanori(:,2), grad.chanori(:,3), ...
        'VariableNames', {'name','Px','Py','Pz','Ox','Oy','Oz'} ...
    );
    
    diff_tol = 25;
    positions_cal_aligned = align_positions(positions_cal, positions_original, diff_tol, true);

    if isempty(positions_cal_aligned)
        warning(['Failed to align to existing positions with tolerance of ', num2str(diff_tol)])
    else
        positions_cal = positions_cal_aligned;
    end
end

% Apply gains (adapted from spm_opm_hfc)
pos = [positions_cal.Px, positions_cal.Py, positions_cal.Pz];
ori = [positions_cal.Ox, positions_cal.Oy, positions_cal.Oz];
cl = positions_cal.name;
channels = S.D.chanlabels;
[~, sel2] = match_str(cl,channels);

fprintf('Creating output dataset\n'); 
outname = fullfile(path(S.D),[S.prefix fname(S.D)]);
cal_D = clone(S.D,outname);
cal_D.save();

grad = [];
grad.label = cl;
grad.coilpos = pos;
grad.coilori = ori;
grad.tra = eye(numel(grad.label));
grad.chanunit = repmat({'T'}, numel(grad.label), 1);
grad.chantype = repmat({'MEGMAG'}, numel(grad.label), 1);
grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');
cal_D = sensors(cal_D, 'MEG', grad);
save(cal_D);

gains = positions_cal.gain;
diag_gains = eye(length(gains)) .* gains;
Yinds = indchannel(S.D,channels(sel2));
if (size(Yinds,1)~=size(diag_gains,1))
    error('data size ~= number of sensors with gain information');
end

chunkSamples = round(S.chunkSize / (8 * size(S.D,1)) * 1e6);
begs = 1:chunkSamples:size(S.D, 2);
ends = (begs + chunkSamples - 1);
if(ends(end) > size(S.D, 2))
    ends(end) = size(S.D, 2);
end

fprintf('%-40s: %30s\n','Processing Data',spm('time'));
for j=1:size(S.D, 3)
    for i =1:length(begs)
        inds = begs(i):ends(i);
        out = S.D(:, inds, j);
        Y = out(Yinds, :);
        out(Yinds,:) = diag_gains * Y;
        cal_D(:, inds, j) = out;
    end
end

% Update forward modelling information
if (S.balance)
    fprintf('%-40s: %30s\n','Updating Sensor Information', spm('time'));
    grad = cal_D.sensors('MEG');
    tmpTra= eye(size(grad.coilori, 1));
    tmpTra(sinds, sinds) = M;
    grad.tra = tmpTra * grad.tra;
    grad.balance.previous = grad.balance.current;
    grad.balance.current = 'cal';
    cal_D = sensors(cal_D,'MEG',grad);
    % Check if any information in D.inv needs updating.
    % TODO: Update to support multiple invs/forwards/modalities
    if isfield(cal_D, 'inv')
        if isfield(cal_D.inv{1}, 'gainmat')
            fprintf(['Clearing current forward model, please recalculate '...
                'with spm_eeg_lgainmat\n']);
            cal_D.inv{1} = rmfield(cal_D.inv{1}, 'gainmat');
        end
        if isfield(cal_D.inv{1},'datareg')
            cal_D.inv{1}.datareg.sensors = grad;
        end
        if isfield(cal_D.inv{1},'forward')
            voltype = cal_D.inv{1}.forward.voltype;
            cal_D.inv{1}.forward = [];
            cal_D.inv{1}.forward.voltype = voltype;
            cal_D = spm_eeg_inv_forward(cal_D, 1);
        end
    end
    cal_D.save();
end
end

%==========================================================================
% - A L I G N   P O S I T I O N S
%==========================================================================
function [positions_aligned, positions_fixed_matched,positions_to_align] = align_positions(positions_to_align, positions_fixed, outlier_dist_threshold, varargin)
% Method to align position information from two different positions files
% by aligning matched pairs of channels. Positions which deviate from the
% general properties of either array are ignored, to improve overall fit.
% 
% Author: Nicholas Alexander (n.alexander@ucl.ac.uk)
% Copyright: Department of Imaging Neuroscience, UCL (2026)

if isempty(varargin)
    plot_output = false;
else
    plot_output = varargin{1};
end

positions_aligned = positions_to_align;
positions_fixed_matched = positions_fixed;

% Match index
[~,b] = ismember(positions_aligned.name,positions_fixed_matched.name);

% Depends on whether slot2sens has been recorded and applied
if ~isempty(b)
    positions_fixed_matched = positions_fixed_matched(b(b ~= 0),:);
    positions_aligned = positions_aligned(b ~= 0, :);
else 
    % This can work when there are minimal unmatched points
    % Get relative distances between channels
    pfm_pos = [positions_fixed_matched.Px, positions_fixed_matched.Py, positions_fixed_matched.Pz];
    pa_pos = [positions_aligned.Px, positions_aligned.Py, positions_aligned.Pz];
    pfm_dists = pdist(pfm_pos);
    pa_dists = pdist(pa_pos);
    
    % Sort the distances and find closest matches
    [~, idx_pfm] = sort(pfm_dists);
    [~, idx_pa] = sort(pa_dists);
    max_matches = min(length(idx_pfm), length(idx_pa));
    matched_idx_pfm = idx_pfm(1:max_matches);
    matched_idx_pta = idx_pa(1:max_matches);
    
    % Identify pairs
    [row_pfm, col_pfm] = ind2sub(size(pfm_pos, 1), matched_idx_pfm);
    [row_pa, col_pa] = ind2sub(size(pa_pos, 1), matched_idx_pta);
    points_pfm = unique([row_pfm; col_pfm]);
    points_pa = unique([row_pa; col_pa]);
    
    % Update table
    num_matches = min(length(points_pfm), length(points_pa));
    positions_fixed_matched = positions_fixed_matched(points_pfm(1:num_matches), :);
    positions_aligned = positions_aligned(points_pa(1:num_matches), :);
end

% Get relative distances between channels
pfm_pos = [positions_fixed_matched.Px, positions_fixed_matched.Py, positions_fixed_matched.Pz];
pa_pos = [positions_aligned.Px, positions_aligned.Py, positions_aligned.Pz];
pfm_dists = squareform(pdist(pfm_pos));
pa_dists = squareform(pdist(pa_pos));

% Identify shared and outlier positions
diff_dists = abs(pfm_dists - pa_dists);
shared_pos = diff_dists == 0;
num_shared_pos = sum(shared_pos,1);
summed_diff = sum(diff_dists,1) ./ (length(diff_dists) - num_shared_pos);
bad_idx_pa = summed_diff > outlier_dist_threshold;
positions_aligned = positions_aligned(~bad_idx_pa,:);

% Match index (again)
[~,b] = ismember(positions_aligned.name,positions_fixed_matched.name);

if ~isempty(b)
    positions_fixed_matched = positions_fixed_matched(b,:);
else
    % Get relative distances between channels
    pfm_pos = [positions_fixed_matched.Px, positions_fixed_matched.Py, positions_fixed_matched.Pz];
    pa_pos = [positions_aligned.Px, positions_aligned.Py, positions_aligned.Pz];
    pfm_dists = pdist(pfm_pos);
    pa_dists = pdist(pa_pos);
    
    % Sort the distances and find closest matches
    [~, idx_pfm] = sort(pfm_dists);
    [~, idx_pa] = sort(pa_dists);
    max_matches = min(length(idx_pfm), length(idx_pa));
    matched_idx_pfm = idx_pfm(1:max_matches);
    matched_idx_pta = idx_pa(1:max_matches);
    
    % Identify pairs
    [row_pfm, col_pfm] = ind2sub(size(pfm_pos, 1), matched_idx_pfm);
    [row_pa, col_pa] = ind2sub(size(pa_pos, 1), matched_idx_pta);
    points_pfm = unique([row_pfm; col_pfm]);
    points_pa = unique([row_pa; col_pa]);
    
    % Update table
    num_matches = min(length(points_pfm), length(points_pa));
    positions_fixed_matched = positions_fixed_matched(points_pfm(1:num_matches), :);
    positions_aligned = positions_aligned(points_pa(1:num_matches), :);
end

if height(positions_aligned) < 10
    return
end

% Subtract centroids
pta_cent(1) = mean(positions_aligned.Px);
pta_cent(2) = mean(positions_aligned.Py);
pta_cent(3) = mean(positions_aligned.Pz);
positions_aligned.Px = positions_aligned.Px - pta_cent(1);
positions_aligned.Py = positions_aligned.Py - pta_cent(2);
positions_aligned.Pz = positions_aligned.Pz - pta_cent(3);

pfm_cent(1) = mean(positions_fixed_matched.Px);
pfm_cent(2) = mean(positions_fixed_matched.Py);
pfm_cent(3) = mean(positions_fixed_matched.Pz);
positions_fixed_matched.Px = positions_fixed_matched.Px - pfm_cent(1);
positions_fixed_matched.Py = positions_fixed_matched.Py - pfm_cent(2);
positions_fixed_matched.Pz = positions_fixed_matched.Pz - pfm_cent(3);

% Also update original
positions_to_align.Px = positions_to_align.Px - pta_cent(1);
positions_to_align.Py = positions_to_align.Py - pta_cent(2);
positions_to_align.Pz = positions_to_align.Pz - pta_cent(3);
pta_pos = [positions_to_align.Px, positions_to_align.Py, positions_to_align.Pz];

% Kabsch
pfm_pos = [positions_fixed_matched.Px, positions_fixed_matched.Py, positions_fixed_matched.Pz];
pa_pos = [positions_aligned.Px, positions_aligned.Py, positions_aligned.Pz];

C = pfm_pos' * pa_pos;
[U, ~, V] = svd(C);
R = V * U';
if det(R) < 0
    V(:,3) = -V(:,3);
    R = V * U';
end

pa_pos = pa_pos * R;
pta_pos = pta_pos * R;

% Optimise transform to minimise error
initial_transform_params = [0, 0, 0, 0, 0, 0];
objective_fun = @(translation_params) error_from_transform(translation_params, pa_pos, pfm_pos);
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
optimal_transform_params = fminunc(objective_fun, initial_transform_params, options);

% Translate
x_translate = optimal_transform_params(1);
y_translate = optimal_transform_params(2);
z_translate = optimal_transform_params(3);

pa_pos(:,1) = pa_pos(:,1) + x_translate;
pa_pos(:,2) = pa_pos(:,2) + y_translate;
pa_pos(:,3) = pa_pos(:,3) + z_translate;

pta_pos(:,1) = pta_pos(:,1) + x_translate;
pta_pos(:,2) = pta_pos(:,2) + y_translate;
pta_pos(:,3) = pta_pos(:,3) + z_translate;

% Rotate
theta_x = optimal_transform_params(4);
theta_y = optimal_transform_params(5);
theta_z = optimal_transform_params(6);
R_x = [1, 0, 0; 0, cos(theta_x), -sin(theta_x); 0, sin(theta_x), cos(theta_x)];
R_y = [cos(theta_y), 0, sin(theta_y); 0, 1, 0; -sin(theta_y), 0, cos(theta_y)];
R_z = [cos(theta_z), -sin(theta_z), 0; sin(theta_z), cos(theta_z), 0; 0, 0, 1];
R = R_x * R_y * R_z;
pa_pos = pa_pos * R;
pta_pos = pta_pos * R;

% Report error
distances = sqrt(sum((pa_pos - pfm_pos).^2, 2));
dist_error = mean(distances);
disp(['Mean position error: ', num2str(dist_error)]);

% Update aligned positions
positions_aligned.Px = pa_pos(:,1) + pfm_cent(1);
positions_aligned.Py = pa_pos(:,2) + pfm_cent(2);
positions_aligned.Pz = pa_pos(:,3) + pfm_cent(3);
C = [positions_aligned.Ox, positions_aligned.Oy, positions_aligned.Oz];
C = C * R;
positions_aligned.Ox = C(:,1);
positions_aligned.Oy = C(:,2);
positions_aligned.Oz = C(:,3);

% And original
positions_to_align.Px = pta_pos(:,1) + pfm_cent(1);
positions_to_align.Py = pta_pos(:,2) + pfm_cent(2);
positions_to_align.Pz = pta_pos(:,3) + pfm_cent(3);
C = [positions_to_align.Ox, positions_to_align.Oy, positions_to_align.Oz];
C = C * R;
positions_to_align.Ox = C(:,1);
positions_to_align.Oy = C(:,2);
positions_to_align.Oz = C(:,3);

% Reapply offset to fixed positions
positions_fixed_matched.Px = pfm_pos(:,1) + pfm_cent(1);
positions_fixed_matched.Py = pfm_pos(:,2) + pfm_cent(2);
positions_fixed_matched.Pz = pfm_pos(:,3) + pfm_cent(3);

% Measure orientation error
U = [positions_fixed_matched.Ox, positions_fixed_matched.Oy, positions_fixed_matched.Oz];
V = [positions_aligned.Ox, positions_aligned.Oy, positions_aligned.Oz];

dot_products = sum(U .* V, 2); 
normU = sqrt(sum(U.^2, 2));
normV = sqrt(sum(V.^2, 2));

cos_theta = dot_products ./ (normU .* normV);
cos_theta = max(min(cos_theta, 1), -1);

angle_radians = acos(cos_theta);
angles = rad2deg(angle_radians);
idxX = startsWith(positions_aligned.name, 'X');
idxY = startsWith(positions_aligned.name, 'Y');
idxZ = startsWith(positions_aligned.name, 'Z');

angle_error(1) = mean(angles(idxX));
angle_error(2) = mean(angles(idxY));
angle_error(3) = mean(angles(idxZ));
disp(['Mean  orientation error: x, ', num2str(angle_error(1)),...
        '; y, ', num2str(angle_error(2)), '; z, ', num2str(angle_error(3))]);

% Visualise
if plot_output
    figure
    hold on
    scale = 10;
    quiver3(positions_fixed.Px, positions_fixed.Py, positions_fixed.Pz,...
        positions_fixed.Ox*scale, positions_fixed.Oy*scale , positions_fixed.Oz*scale,...
        'off','Color',[166,206,227]./255);
    quiver3(positions_fixed_matched.Px, positions_fixed_matched.Py, positions_fixed_matched.Pz,...
        positions_fixed_matched.Ox*scale, positions_fixed_matched.Oy*scale , positions_fixed_matched.Oz*scale,...
        'off','Color',[31,120,180]./255);
    quiver3(positions_to_align.Px, positions_to_align.Py, positions_to_align.Pz,...
        positions_to_align.Ox*scale, positions_to_align.Oy*scale, positions_to_align.Oz*scale,...
        'off','Color',[178,223,138]./255);
    quiver3(positions_aligned.Px, positions_aligned.Py, positions_aligned.Pz,...
        positions_aligned.Ox*scale, positions_aligned.Oy*scale, positions_aligned.Oz*scale,...
        'off','Color',[51,160,44]./255);
    legend({'Unmatched Fixed','Fixed','Bad aligned','Good aligned'})
    axis equal

    pf_pos = [positions_fixed.Px, positions_fixed.Py, positions_fixed.Pz];
    limits = (max(abs(pf_pos)).*1.1);
    xlim([-limits(1) limits(1)])
    ylim([-limits(2) limits(2)])
    zlim([-limits(3) limits(3)])
    hold off
    rotate3d
end
end

%==========================================================================
% - E R R O R   F R O M   T R A N S F O R M
%==========================================================================
function error = error_from_transform(transform_params, pos_1, pos_2)
    % Translate
    xTran = transform_params(1);
    yTran = transform_params(2);
    zTran = transform_params(3);

    pos_1(:,1) = pos_1(:,1) + xTran;
    pos_1(:,2) = pos_1(:,2) + yTran;
    pos_1(:,3) = pos_1(:,3) + zTran;

    % Rotate
    thetaX = transform_params(4);
    thetaY = transform_params(5);
    thetaZ = transform_params(6);
    Rx = [1, 0, 0; 0, cos(thetaX), -sin(thetaX); 0, sin(thetaX), cos(thetaX)];
    Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
    Rz = [cos(thetaZ), -sin(thetaZ), 0; sin(thetaZ), cos(thetaZ), 0; 0, 0, 1];
    Ropt = Rx * Ry * Rz;
    transformed_pos_1 = pos_1 * Ropt;

    % Compute error
    distances = sqrt(sum((transformed_pos_1 - pos_2).^2, 2));
    error = mean(distances);
end