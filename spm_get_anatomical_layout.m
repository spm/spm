function [lay] = spm_get_anatomical_layout(sensor_positions, sensor_labels, head_positions, fiducials, varargin)
% Produce anatomically valid 2D representation of 3D sensor positions
% FORMAT [lay] = spm_get_anatomical_layout(sensor_positions, sensor_labels, head_positions, fiducials, plot_output)
%
% Input Parameters:
%     sensor_positions: An nx3 matrix representing the 3D Cartesian
%                       coordinates (x, y, z) of n sensors.
%     sensor_labels:    An nx1 cell array containing labels for each sensor.
%     head_positions:   An mx3 matrix representing the coordinates of
%                       positions on the head.
%     fiducials:        Struct containing coordinates of fiducials.
%         fiducials.NAS = [1 x 3]
%         fiducials.LPA = [1 x 3]
%         fiducials.RPA = [1 x 3]
%         fiducials.INI = [1 x 3] - Optional
%                       Note, if fiducials.INI is not specified it will be
%                       estimated.
%     plot_output:      (optional) A logical value indicating whether to generate a
%                       plot of measurements taken in 3D, as well as the
%                       output layout.
%     unit:             (optional) A string, e.g. 'mm', will be estimated otherwise.
%     sensor_orientations:  (optional) An nx3 matrix of orientation vectors of n
%                           sensors
%
% Output:
%     lay: A structure representing the anatomical sensor layout in the
%         FieldTrip style. It contains the following fields:
%         lay.pos:      A nx2 matrix representing the 2D coordinates of the
%                       sensors in the layout. Each column contains the
%                       [x, y] coordinates of a sensor relative to the
%                       vertex (Cz).
%         lay.label:    A cell array containing labels for each sensor.
%         lay.outline:  A cell array of matrices representing the 2D
%                       coordinates of the head surface outline, ears and
%                       nose.
%         lay.mask:     A cell array containing a matrix of positions which
%                       draw a convex hull around lay.pos to mask grid
%                       positions which would otherwise be extrapolated to
%                       when plotting.
%                       To allow extrapolation to a full circle, try:
%                       lay.mask{1} = lay.outline{1}.
%         lay.ori:      If sensor_orientations are provided then equivalent
%                       2D orientations are calculated. 
%
%__________________________________________________________________________
%
% Further help:
%
% spm_get_anatomical_layout is a function that defines the scalp surface
% according to the approximate polar grid used to place electrodes in the
% 10-20 system. The method then measures the position of on-scalp sensors
% in relation to this polar grid (angle and eccentricity) and applies this
% to a standard 2D polar grid. For full details, or to cite this method see
% Alexander et al (2025) available here: https://doi.org/10.1111/ejn.70060.
%
% The function performs the following steps:
%
%   Get Anterior-Posterior and Left-Right Vectors:
%     Vectors defining the anterior-posterior and left-right directions on
%     the head surface are calculated. The vertex (Cz) position on the head
%     is then estimated, such that, for measurements taken across the scalp,
%     vertex->left, vertex->right are of equal length and vertex->anterior,
%     vertex-posterior are of equal lengths.
%
%   Create the approximate polar grid across the scalp.
%     Lines across the scalp are made based on the vectors described above.
%     These lines are then shortened according to the 10-20 method to
%     produce a grid with Fz, Oz, T3, T4 positions are the periphery. The
%     circumference around these points is then defined.
%
%   Define Sensor Position Relative to Cz:
%     Sensor positions are projected onto the scalp surface and their
%     position in relation to the vertex and circumference are defined.
%     Specifically, the distance from Cz to projected sensor position
%     relative to the distance from Cz to the circumference along the same
%     axis is taken. The relative distance from the point this axis crosses
%     the circumference to the nearest peripheral landmarks (e.g. Fz, T3)
%     determines the angle.
%
%   Define Sensor Position on a 2D Circle:
%     Using the eccentricity and angle measurements from the previous step,
%     the position of each sensor is reproduced on a 2D polar grid. This is
%     then formatted as a FieldTrip style layout structure with a nose,
%     ears, outline and mask.
%
%   (optional) Represent Sensor Orientations in 2D:
%     A tangential basis is defined in both 3D and 2D as an line (across the
%     scalp) from RPA to LPA, passing through the sensor position. The
%     local direction of the arc at the sensor position defined the primary
%     tangential axis. Radial is defined as a vector towards the origin.
%     Secondary primary axis as the cross product of primary and radial
%     axes. In 2D, magnitude is scaled by radial similarity. 
%_________________________________________________________________________

% Nicholas Alexander
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


%==========================================================================
% - S P M   G E T   A N A T O M I C A L   L A Y O U T
%==========================================================================
%-TIDY up the inputs
%==========================================================================
%-REMOVE positions that are unlikely to be on the head
%--------------------------------------------------------------------------
% For example, shoulders from a 3D scan.
% Origin is between the ears. This is fairly arbitrary, but in line with
% many coordinate systems.
origin = mean([fiducials.LPA; fiducials.RPA], 1);

% Get the distance from the origin to each sensor position
sensor_distances = sqrt(sum((sensor_positions - origin).^2, 2));

% Set a maximum distance proportional to the max distance to sensor
max_distance = max(sensor_distances) * 1.2;

% Calculate the Euclidean distance between the origin and each position
head_distances = sqrt(sum((head_positions - origin).^2, 2));

% Find the positions that are within the maximum distance
head_positions = head_positions(head_distances <= max_distance, :);

if (length(head_distances(:, 1)) - length(head_positions(:, 1))) > 0
    disp(['Removed ', num2str((length(head_distances(:, 1)) - length(head_positions(:, 1)))), ...
        ' head positions that were far from the centre of the head.'])
    disp('If you are having issues, try plotting the head positions along with the sensors.');
end

% Unpack varargin
if numel(varargin) >= 1
    plot_output = logical(varargin{1});
else
    plot_output = false;
end
if numel(varargin) >= 2
    unit = varargin{2};
else
    unit = [];
end
if numel(varargin) >= 3
    sensor_orientations = varargin{3};
else
    sensor_orientations = [];
end

% Estimate units
if isempty(unit)
    largest_dim = 2 * max(sensor_distances);
    if largest_dim > 200
        unit = 'mm';
        unit_dec_adjust = 3;
    elseif largest_dim > 20
        unit = 'cm';
        unit_dec_adjust = 2;
    elseif largest_dim > 0.2
        unit = 'm';
        unit_dec_adjust = 0;
    else
        unit = 'm';
        unit_dec_adjust = 0;
        warning('Unit estimation uncertain. Defaulting to metres.')
    end
else
    switch unit
        case 'mm'
            unit_dec_adjust = 3;
        case 'cm'
            unit_dec_adjust = 2;
        case 'm'
            unit_dec_adjust = 0;
        otherwise
            disp('Incorrect unit specified.')
    end
end
dec_point = 6 - unit_dec_adjust;

% No 2D position should be outside -1:1
max_image_size = 2;

%-ESTIMATE the position of the inion if not provided
%--------------------------------------------------------------------------
if ~isfield(fiducials, 'INI')
    % Start with the anterior vector (i.e. origin to NAS)
    a_vec = (fiducials.NAS - origin) / norm((fiducials.NAS - origin));

    % Use the right vector as an axis of rotation
    r_vec = (fiducials.RPA - origin) / norm(fiducials.RPA - origin);

    % 160 degrees seems about right, only based on 10 or so tests, so take
    % with a pinch of salt, or label the inion.
    theta = 160 * pi / 180;

    % Build the rotation matrix
    c = cos(theta);
    s = sin(theta);
    t = 1 - c;
    x = r_vec(1);
    y = r_vec(2);
    z = r_vec(3);

    rotation_matrix = [t * x * x + c, t * x * y - z * s, t * x * z + y * s;
                    t * x * y + z * s, t * y * y + c, t * y * z - x * s;
                    t * x * z - y * s, t * y * z + x * s, t * z * z + c];

    % Estimate the posterior vector
    p_vec = (rotation_matrix * a_vec')';

    % Make a model of the scalp from head positions and the fiducials.
    head_positions = [head_positions; fiducials.NAS; fiducials.LPA; fiducials.RPA];
    tri = convhulln([head_positions(:, 1), head_positions(:, 2), head_positions(:, 3)]);

    % Find where a_vec intersects with the mesh. This also adds in the INI to
    % the headPositions.
    [fiducials.INI, head_positions, ~] = find_line_triangle_intersection(p_vec, origin, head_positions, tri, dec_point);
else
    % Otherwise, just update the head_positions to include fiducials.
    head_positions = [head_positions; fiducials.NAS; fiducials.INI; fiducials.LPA; fiducials.RPA];
end

%-GET the set of vectors used to describe measurements across the scalp
%--------------------------------------------------------------------------
% left (origin to LPA)
l_vec = (fiducials.LPA - origin) / norm(fiducials.LPA - origin);

% right (origin to RPA)
r_vec = (fiducials.RPA - origin) / norm(fiducials.RPA - origin);

% anterior (origin to NAS)
a_vec = (fiducials.NAS - origin) / norm((fiducials.NAS - origin));

% posterior (origin to INI)
p_vec = (fiducials.INI - origin) / norm((fiducials.INI - origin));

%-CREATE a head model using a convex hull method
%--------------------------------------------------------------------------
tri = convhulln([head_positions(:, 1), head_positions(:, 2), head_positions(:, 3)]);

% Edges are much easier to work with
[edges] = tris_to_edges(tri);

% Check that the fiducials are on the boundary of the convhull. Adjust
% accordingly and report to user
used_head_position_idx = unique(tri);
unused_head_positions = head_positions;
unused_head_positions(used_head_position_idx, :) = [];
if ismember(fiducials.LPA, unused_head_positions, 'rows')
    tmpPos = fiducials.LPA;
    [fiducials.LPA, head_positions, tri] = find_line_triangle_intersection(l_vec, origin, head_positions, tri, dec_point);
    disp(['LPA not on boundary of convex head surface. Adjusted by ', num2str(pdist([tmpPos;fiducials.LPA])), unit])
end
if ismember(fiducials.RPA, unused_head_positions, 'rows')
    tmpPos = fiducials.RPA;
    [fiducials.RPA, head_positions, tri] = find_line_triangle_intersection(r_vec, origin, head_positions, tri, dec_point);
    disp(['RPA not on boundary of convex head surface. Adjusted by ', num2str(pdist([tmpPos;fiducials.RPA])), unit])
end
if ismember(fiducials.INI, unused_head_positions, 'rows')
    tmpPos = fiducials.INI;
    [fiducials.INI, head_positions, tri] = find_line_triangle_intersection(p_vec, origin, head_positions, tri, dec_point);
    disp(['INI not on boundary of convex head surface. Adjusted by ', num2str(pdist([tmpPos;fiducials.INI])), unit])
end
if ismember(fiducials.NAS, unused_head_positions, 'rows')
    tmpPos = fiducials.NAS;
    [fiducials.NAS, head_positions, tri] = find_line_triangle_intersection(a_vec, origin, head_positions, tri, dec_point);
    disp(['NAS not on boundary of convex head surface. Adjusted by ', num2str(pdist([tmpPos;fiducials.NAS])), unit])
end

% Update origin and vectors as they may have changed. Probably pointless?
origin = mean([fiducials.LPA; fiducials.RPA], 1);
l_vec = (fiducials.LPA - origin) / norm(fiducials.LPA - origin);
r_vec = (fiducials.RPA - origin) / norm(fiducials.RPA - origin);
a_vec = (fiducials.NAS - origin) / norm((fiducials.NAS - origin));
p_vec = (fiducials.INI - origin) / norm((fiducials.INI - origin));

%-LABEL the superior position
%--------------------------------------------------------------------------
% Estimate the superior vector as the mean of possible solutions
s_vec(1, :) = -cross(l_vec, a_vec);
s_vec(2, :) = cross(l_vec, p_vec);
s_vec(3, :) = cross(r_vec, a_vec);
s_vec(4, :) = -cross(r_vec, p_vec);
s_vec = mean(s_vec, 1) / norm(mean(s_vec, 1));

% Refine estimate of superior point by trying to make left and right sides
% equal length and posterior and anterior equal length
sup_error_threshold = 0.001 * 10^unit_dec_adjust;
rl_adjustment_mod = sup_error_threshold * 10;
pa_adjustment_mod = sup_error_threshold * 10;

% Set an upper limit on iterations
iteration = 0;
iteration_limit = 500;
while rl_adjustment_mod > sup_error_threshold || pa_adjustment_mod > sup_error_threshold && iteration < iteration_limit
    iteration = iteration + 1;

    % Add the current estimate of the superior position
    [fiducials.SUP, tmp_head_pos, tmp_tri] = find_line_triangle_intersection(s_vec, origin, head_positions, tri, dec_point);
    s_vec = (fiducials.SUP - origin) / norm(fiducials.SUP - origin);
    % Redefine tri as edges
    [tmp_edges] = tris_to_edges(tmp_tri);
    
    % Take measurements across the scalp from Cz, emulating 10-20 setup
    % left
    l_intersect_points = get_surface_points_about_plane(cross(l_vec, s_vec), origin, ...
                        cross(l_vec, p_vec), origin, tmp_head_pos, tmp_edges, fiducials.LPA, false, dec_point);
    tmp_start = find_fraction_along_line(l_intersect_points, fiducials.LPA, dec_point);
    tmp_end = find_fraction_along_line(l_intersect_points, fiducials.SUP, dec_point);
    l_intersect_points = cut_line_ends(l_intersect_points, [tmp_start, tmp_end], dec_point);

    % right
    r_intersect_points = get_surface_points_about_plane(cross(r_vec, s_vec), origin, ...
                        cross(r_vec, a_vec), origin, tmp_head_pos, tmp_edges, fiducials.RPA, false, dec_point);
    tmp_start = find_fraction_along_line(r_intersect_points, fiducials.RPA, dec_point);
    tmp_end = find_fraction_along_line(r_intersect_points, fiducials.SUP, dec_point);
    r_intersect_points = cut_line_ends(r_intersect_points, [tmp_start, tmp_end], dec_point);

    % anterior
    a_intersect_points = get_surface_points_about_plane(-cross(a_vec, s_vec), origin, ...
                        cross(a_vec, l_vec), origin, tmp_head_pos, tmp_edges, fiducials.NAS, false, dec_point);
    tmp_start = find_fraction_along_line(a_intersect_points, fiducials.NAS, dec_point);
    tmp_end = find_fraction_along_line(a_intersect_points, fiducials.SUP, dec_point);
    a_intersect_points = cut_line_ends(a_intersect_points, [tmp_start, tmp_end], dec_point);

    % posterior
    p_intersect_points = get_surface_points_about_plane(-cross(p_vec, s_vec), origin, ...
                        cross(p_vec, r_vec), origin, tmp_head_pos, tmp_edges, fiducials.INI, false, dec_point);
    tmp_start = find_fraction_along_line(p_intersect_points, fiducials.INI, dec_point);
    tmp_end = find_fraction_along_line(p_intersect_points, fiducials.SUP, dec_point);
    p_intersect_points = cut_line_ends(p_intersect_points, [tmp_start, tmp_end], dec_point);

    % Check lengths
    lLength = length_of_lines(l_intersect_points);
    rLength = length_of_lines(r_intersect_points);
    aLength = length_of_lines(a_intersect_points);
    pLength = length_of_lines(p_intersect_points);

    % Calculate error. Note, adjustment divided by 2 to account for curvature
    % of the scalp. Otherwise the adjustment will overshoot.
    rl_adjustment_mod = abs(rLength-lLength);
    pa_adjustment_mod = abs(pLength-aLength);

    if rl_adjustment_mod > sup_error_threshold
        % Left/right adjustment
        if lLength > rLength
            fiducials.SUP = fiducials.SUP + (l_vec*rl_adjustment_mod)/2;
        else
            fiducials.SUP = fiducials.SUP + (r_vec*rl_adjustment_mod)/2;
        end
    end

    if pa_adjustment_mod > sup_error_threshold
        % Anterior/posterior adjustment
        if aLength > pLength
            fiducials.SUP = fiducials.SUP + (a_vec*pa_adjustment_mod)/2;
        else
            fiducials.SUP = fiducials.SUP + (p_vec*pa_adjustment_mod)/2;
        end
    end

    % Update superior estimate and head positions.
    if rl_adjustment_mod <= sup_error_threshold || pa_adjustment_mod <= sup_error_threshold
        head_positions = tmp_head_pos;
        edges = tmp_edges;
    end

    s_vec = (fiducials.SUP - origin) / norm(fiducials.SUP - origin);
end

if iteration >= iteration_limit
    warning('Iteration limit for superior vector estimation reached.')
end

% Plot progress if requested.
if plot_output
    figure; hold on; axis equal;

    quiver3(origin(1), origin(2), origin(3), l_vec(1), l_vec(2), l_vec(3), 200);
    quiver3(origin(1), origin(2), origin(3), r_vec(1), r_vec(2), r_vec(3), 200);
    quiver3(origin(1), origin(2), origin(3), a_vec(1), a_vec(2), a_vec(3), 200);
    quiver3(origin(1), origin(2), origin(3), p_vec(1), p_vec(2), p_vec(3), 200);
    quiver3(origin(1), origin(2), origin(3), s_vec(1), s_vec(2), s_vec(3), 200);

    % Add fiducials to plot
    fid_labels = {'NAS', 'INI', 'LPA', 'RPA', 'SUP'};
    for i = 1:length(fid_labels)
        scatter3(fiducials.(fid_labels{i})(1), fiducials.(fid_labels{i})(2), fiducials.(fid_labels{i})(3), 'r');
    end
end

%-MAKE the polar grid on the head
%==========================================================================
%-GET 10-20 electrode positions for the circumference
%--------------------------------------------------------------------------
% Move up by 10% either side of the pa and lr intersects
l_intersect_points = cut_line_ends(l_intersect_points, [0.2, 1], dec_point);
r_intersect_points = cut_line_ends(r_intersect_points, [0.2, 1], dec_point);
a_intersect_points = cut_line_ends(a_intersect_points, [0.2, 1], dec_point);
p_intersect_points = cut_line_ends(p_intersect_points, [0.2, 1], dec_point);

% Now get the points around the circumference
elec.Oz = p_intersect_points(1, :);
elec.FPz = a_intersect_points(1, :);
elec.T3 = l_intersect_points(1, :);
elec.T4 = r_intersect_points(1, :);
elec.Cz = round(fiducials.SUP, dec_point, 'decimals');

% Add these to the plot, if there is one
if plot_output    
    plot_ordered_points(l_intersect_points);
    plot_ordered_points(r_intersect_points);
    plot_ordered_points(a_intersect_points);
    plot_ordered_points(p_intersect_points);

    elec_labels = fieldnames(elec);

    % Add fiducials
    for i = 1:length(elec_labels)
        scatter3(elec.(elec_labels{i})(1), elec.(elec_labels{i})(2), elec.(elec_labels{i})(3), 'k');
    end
end

% Move to an electrode structure (replacing fiducials as reference)
origin_elec = mean([elec.T3; elec.T4], 1);
s_vec = (elec.Cz - origin_elec) / norm(elec.Cz - origin_elec);
a_vec = (elec.FPz - origin_elec) / norm(elec.FPz - origin_elec);
p_vec = (elec.Oz - origin_elec) / norm(elec.Oz - origin_elec);
l_vec = (elec.T3 - origin_elec) / norm(elec.T3 - origin_elec);
r_vec = (elec.T4 - origin_elec) / norm(elec.T4 - origin_elec);

%-UPDATE the head definition to include these new points
%--------------------------------------------------------------------------
% left
% Create new edges for this intersection
startEdgeIdx = 1 + length(head_positions(:, 1));
endEdgeIdx = startEdgeIdx + length(l_intersect_points(:, 1)) - 1;
baseEdges = [(startEdgeIdx:1:(endEdgeIdx - 1))', ((startEdgeIdx + 1):1:endEdgeIdx)'];

% Add the points to the surface
head_positions = [head_positions; l_intersect_points];

% And concatenate edge lists
edges = [edges; baseEdges];

% right
% Create new edges for this intersection
startEdgeIdx = 1 + length(head_positions(:, 1));
endEdgeIdx = startEdgeIdx + length(r_intersect_points(:, 1)) - 1;
baseEdges = [(startEdgeIdx:1:(endEdgeIdx - 1))', ((startEdgeIdx + 1):1:endEdgeIdx)'];

% Add the points to the surface
head_positions = [head_positions; r_intersect_points];

% And concatenate edge lists
edges = [edges; baseEdges];

% anterior
% Create new edges for this intersection
startEdgeIdx = 1 + length(head_positions(:, 1));
endEdgeIdx = startEdgeIdx + length(a_intersect_points(:, 1)) - 1;
baseEdges = [(startEdgeIdx:1:(endEdgeIdx - 1))', ((startEdgeIdx + 1):1:endEdgeIdx)'];

% Add the points to the surface
head_positions = [head_positions; a_intersect_points];

% And concatenate edge lists
edges = [edges; baseEdges];

% posterior
% Create new edges for this intersection
startEdgeIdx = 1 + length(head_positions(:, 1));
endEdgeIdx = startEdgeIdx + length(p_intersect_points(:, 1)) - 1;
baseEdges = [(startEdgeIdx:1:(endEdgeIdx - 1))', ((startEdgeIdx + 1):1:endEdgeIdx)'];

% Add the points to the surface
head_positions = [head_positions; p_intersect_points];

% And concatenate edge lists
edges = [edges; baseEdges];

%-CREATE the circumference intersections
%--------------------------------------------------------------------------
% left to anterior
la_vec = mean([l_vec; a_vec]);
la_intersect_points = get_surface_points_about_plane(cross(a_vec, l_vec), ...
                    origin_elec, la_vec, origin_elec, head_positions, edges, elec.Oz, false, dec_point);
tmp_start = find_fraction_along_line(la_intersect_points, elec.T3, dec_point);
tmp_end = find_fraction_along_line(la_intersect_points, elec.FPz, dec_point);
la_intersect_points = cut_line_ends(la_intersect_points, [tmp_start, tmp_end], dec_point);

% anterior to right
ar_vec = mean([a_vec; r_vec]);
ar_intersect_points = get_surface_points_about_plane(cross(r_vec, a_vec), ...
                    origin_elec, ar_vec, origin_elec, head_positions, edges, elec.T3, false, dec_point);
tmp_start = find_fraction_along_line(ar_intersect_points, elec.FPz, dec_point);
tmp_end = find_fraction_along_line(ar_intersect_points, elec.T4, dec_point);
ar_intersect_points = cut_line_ends(ar_intersect_points, [tmp_start, tmp_end], dec_point);

% right to posterior
rp_vec = mean([r_vec; p_vec]);
rp_intersect_points = get_surface_points_about_plane(cross(p_vec, l_vec), ...
                    origin_elec, rp_vec, origin_elec, head_positions, edges, elec.FPz, false, dec_point);
tmp_start = find_fraction_along_line(rp_intersect_points, elec.T4, dec_point);
tmp_end = find_fraction_along_line(rp_intersect_points, elec.Oz, dec_point);
rp_intersect_points = cut_line_ends(rp_intersect_points, [tmp_start, tmp_end], dec_point);

% posterior to left
pl_vec = mean([p_vec; l_vec]);
pl_intersect_points = get_surface_points_about_plane(cross(l_vec, p_vec), ...
                    origin_elec, pl_vec, origin_elec, head_positions, edges, elec.T4, false, dec_point);
tmp_start = find_fraction_along_line(pl_intersect_points, elec.Oz, dec_point);
tmp_end = find_fraction_along_line(pl_intersect_points, elec.T3, dec_point);
pl_intersect_points = cut_line_ends(pl_intersect_points, [tmp_start, tmp_end], dec_point);

% Add the circumference to the plot
if plot_output
    plot_ordered_points(la_intersect_points);
    plot_ordered_points(ar_intersect_points);
    plot_ordered_points(rp_intersect_points);
    plot_ordered_points(pl_intersect_points);
end

% Add the points to the surface
head_positions = [head_positions; la_intersect_points; ar_intersect_points; ...
                rp_intersect_points; pl_intersect_points];


%-REPRESENT the sensor positions in 2D based on the polar grid
%==========================================================================

% Prepare output
lay = [];
num_sens = length(sensor_positions(:, 1));
basis_frame = nan(num_sens,9);
% Loop through each sensor position
for sens_idx = 1:num_sens
    % Make a vector from the origin to the sensor
    sens_vec = (sensor_positions(sens_idx, :) - origin_elec) / norm(sensor_positions(sens_idx, :) - origin_elec);

    % Get the point sens_vec crosses the scalp
    [sens_pos, tmp_head_pos, tmp_tri] = find_line_triangle_intersection(sens_vec, origin_elec, head_positions, tri, dec_point);

    % Redefine tri as edges
    [tmp_edges] = tris_to_edges(tmp_tri);

    % Add sensor positions to the plot, including projection
    if plot_output
        scatter3(sensor_positions(sens_idx, 1), sensor_positions(sens_idx, 2), sensor_positions(sens_idx, 3), 50, 'b');
        scatter3(sens_pos(1), sens_pos(2), sens_pos(3), 100);
        text(sensor_positions(sens_idx, 1), sensor_positions(sens_idx, 2), sensor_positions(sens_idx, 3), sensor_labels{sens_idx});
    end

    % Make a slice from Cz to the sensor
    sens_sup_vec = mean([sens_vec;s_vec]) / norm(mean([sens_vec;s_vec]));
    tmpVec = cross(sens_vec, s_vec) / norm(cross(sens_vec, s_vec));
    [sens_intersect_points] = get_surface_points_about_plane(tmpVec, origin_elec, sens_sup_vec, origin_elec, tmp_head_pos, tmp_edges, elec.Cz, false, dec_point);
    [tmp_start, ~] = find_fraction_along_line(sens_intersect_points, elec.Cz, dec_point);
    [tmp_end, ~] = find_fraction_along_line(sens_intersect_points, sens_pos, dec_point);
    sens_intersect_points_orig = sens_intersect_points;
    sens_intersect_points = cut_line_ends(sens_intersect_points, [tmp_start, tmp_end], dec_point);

    % Add the sensor intersection to the plot
    if plot_output
        plot_ordered_points(sens_intersect_points)
    end

    % Find where the sens intersect crosses the circumference
    % la
    tmp_length = length(la_intersect_points);
    la_cross_point = get_plane_surface_intersect([(1:tmp_length-1)', (2:tmp_length)'], la_intersect_points, cross(sens_vec, s_vec), origin_elec, dec_point);
    la_cross_point = remove_points_below_plane(la_cross_point, sens_sup_vec, origin_elec);

    % ar
    tmp_length = length(ar_intersect_points);
    ar_cross_point = get_plane_surface_intersect([(1:tmp_length-1)', (2:tmp_length)'], ar_intersect_points, cross(sens_vec, s_vec), origin_elec, dec_point);
    ar_cross_point = remove_points_below_plane(ar_cross_point, sens_sup_vec, origin_elec);

    % rp
    tmp_length = length(rp_intersect_points);
    rp_cross_point = get_plane_surface_intersect([(1:tmp_length-1)', (2:tmp_length)'], rp_intersect_points, cross(sens_vec, s_vec), origin_elec, dec_point);
    rp_cross_point = remove_points_below_plane(rp_cross_point, sens_sup_vec, origin_elec);

    % pl
    tmp_length = length(pl_intersect_points);
    pl_cross_point = get_plane_surface_intersect([(1:tmp_length-1)', (2:tmp_length)'], pl_intersect_points, cross(sens_vec, s_vec), origin_elec, dec_point);
    pl_cross_point = remove_points_below_plane(pl_cross_point , sens_sup_vec, origin_elec);
    

    if sum(~cellfun('isempty', {la_cross_point, ar_cross_point, rp_cross_point, pl_cross_point})) ~= 1
        % Decide which way to go
        dir_vec = (sens_pos - elec.Cz) / norm((sens_pos - elec.Cz));

        % List of fiducial names
        elec_names = {'T4', 'FPz', 'T4', 'Oz'};
        
        % Calculate cosine similarities for each direction
        cos_similarities = zeros(1, numel(elec_names));
        for i = 1:numel(elec_names)
            fiducial = elec.(elec_names{i});

            magnitude_fiducial = sqrt(sum(fiducial.^2, 2));
            dot_product = sum(fiducial .* dir_vec, 2);

            magnitude_dir_vec = sqrt(sum(dir_vec.^2, 2));
            cos_similarities(i) = dot_product ./ (magnitude_fiducial .* magnitude_dir_vec);
        end
        
        % Very rarely the plane can slice exactly at an intersection
        if sum(~cellfun('isempty', {la_cross_point, ar_cross_point, rp_cross_point, pl_cross_point})) == 0
            la_cross_point = elec.T4;
            ar_cross_point = elec.FPz;
            rp_cross_point = elec.T3;
            pl_cross_point = elec.Oz;
            
        end

        if cos_similarities(4) > cos_similarities(2) % posterior > anterior
            ar_cross_point = [];
            la_cross_point = [];
            if cos_similarities(1) > cos_similarities(3) ||...
                (~isempty(rp_cross_point) && isempty(pl_cross_point)) % left > right
                % rp
                pl_cross_point = [];
            else
                % pl
                rp_cross_point = [];
            end
        else
            pl_cross_point = [];
            rp_cross_point = [];
            if cos_similarities(1) > cos_similarities(3) ||...
                (~isempty(rp_cross_point) && isempty(la_cross_point))
                % ar
                la_cross_point = [];
            else
                % la
                rp_cross_point = [];
            end
        end
    end

    % Set a starting point about the circle
    c_intersect_points = [];
    cross_point = [];
    circ_mod = NaN;
    if ~isempty(la_cross_point)
        circ_mod = 0;
        cross_point = la_cross_point;
        c_intersect_points = la_intersect_points;
    elseif ~isempty(ar_cross_point)
        circ_mod = 0.25;
        cross_point = ar_cross_point;
        c_intersect_points = ar_intersect_points;
    elseif ~isempty(rp_cross_point)
        circ_mod = 0.5;
        cross_point = rp_cross_point;
        c_intersect_points = rp_intersect_points;
    elseif ~isempty(pl_cross_point)
        circ_mod = 0.75;
        cross_point = pl_cross_point;
        c_intersect_points = pl_intersect_points;
    end

    % Find the fraction along the circumference the sensor vector crossed
    circ_frac = find_fraction_along_line(c_intersect_points, cross_point, dec_point);
    circ_frac = (circ_frac / 4) + circ_mod;

    % And the distance form the centre as a fraction of distance to circ.
    [tmp_start, sens_intersect_points_orig] = find_fraction_along_line(sens_intersect_points_orig, elec.Cz, dec_point);
    [tmp_end, sens_intersect_points_orig] = find_fraction_along_line(sens_intersect_points_orig, cross_point, dec_point);
    c_intersect_points = cut_line_ends(sens_intersect_points_orig, [tmp_start, tmp_end], dec_point);
    
    c_cist_cumsum = cumsum(sqrt(sum(diff(c_intersect_points, [], 1).^2, 2)));
    c_dist = c_cist_cumsum(end);

    % and the sensor distance
    sens_dist_cumsum = cumsum(sqrt(sum(diff(sens_intersect_points, [], 1).^2, 2)));
    sens_dist = sens_dist_cumsum(end);

    cent_frac = sens_dist / c_dist;

    % Find point on circle
    theta = 2 * pi * circ_frac;
    x = -0.5 * cos(theta);
    y = 0.5 * sin(theta);

    % Plot point on line at distFrac from center
    x_dist = cent_frac * x / max_image_size;
    y_dist = cent_frac * y / max_image_size;

    lay.pos(sens_idx, 1:2) = [x_dist, y_dist];
    lay.label{sens_idx, 1} = sensor_labels{sens_idx};
    lay.width(sens_idx, 1) = 0.076/ max_image_size;
    lay.height(sens_idx, 1) = 0.064/ max_image_size;

    if ~isempty(sensor_orientations)
        % 2D orientations based on tangential basis set
        % Primary tan slice
        lr_vec = (fiducials.LPA - fiducials.RPA) / norm(fiducials.LPA - fiducials.RPA);
        sens_vec_lr = (sens_pos - origin) / norm(sens_pos - origin);
        sens_vec_lr_perp_to_lr = sens_vec_lr * dot(lr_vec, lr_vec) - lr_vec * dot(lr_vec, sens_vec_lr);
        l_intersect_points = get_surface_points_about_plane(cross(lr_vec, sens_vec_lr), origin, ...
                            sens_vec_lr_perp_to_lr, origin, tmp_head_pos, tmp_edges, fiducials.LPA, false, dec_point);
        [tmp_start, ~] = find_fraction_along_line(l_intersect_points, fiducials.LPA, dec_point);
        [tmp_end, ~] = find_fraction_along_line(l_intersect_points, fiducials.RPA, dec_point);
        if tmp_start > tmp_end
            l_intersect_points = flipud(l_intersect_points);
        end
    
        % Find tan axis as line segment
        [~,cur_sens_idx] = ismember(sens_pos, l_intersect_points, 'rows');
    
        tmp = (l_intersect_points(cur_sens_idx-1,:) - l_intersect_points(cur_sens_idx,:)) / norm((l_intersect_points(cur_sens_idx-1,:) - l_intersect_points(cur_sens_idx,:)));
        tmp2 = (l_intersect_points(cur_sens_idx,:) - l_intersect_points(cur_sens_idx+1,:)) / norm((l_intersect_points(cur_sens_idx,:) - l_intersect_points(cur_sens_idx+1,:)));
        
        initial_primary_tan_vec = mean([tmp;tmp2]) / norm(mean([tmp;tmp2]));
    
        % Radial vector
        sens_vec = (sens_pos - origin_elec) / norm(sens_pos - origin_elec);
        radial_vec = -sens_vec;
        radial_vec = radial_vec / norm(radial_vec);
    
        % Orthonormal frame
        primary_tan_vec = initial_primary_tan_vec - dot(initial_primary_tan_vec, radial_vec) * radial_vec;
        primary_tan_vec = primary_tan_vec / norm(primary_tan_vec);  % Now orthogonal and unit length
        secondary_tan_vec = cross(radial_vec, primary_tan_vec);
    
        basis_frame(sens_idx,1:9) = [radial_vec(:)' primary_tan_vec(:)' secondary_tan_vec(:)'];
    end
end

%-STRUCTURE the output layout
%==========================================================================
% Add an outline 10% out for 10:20 compatibility
samples = 128;
theta = linspace(0, 2 * pi, samples);
outer_r = 1.25;
rad2 = 0.5 * outer_r / max_image_size;
x2 =  rad2 * cos(theta);
y2 =  rad2 * sin(theta);
lay.outline{1, 1} = [x2', y2'];

% Add mask
mask_idx = convhull(lay.pos(:, 1), lay.pos(:, 2));
lay.mask{1} = [lay.pos(mask_idx, :);lay.pos(mask_idx(1), :)];

% Add a nose to the outline
nose_scale = rad2 * 1.1;
nose_half_width = 0.07 / max_image_size;
[~, index] = min(abs(theta - nose_half_width));
lay.outline{1, 2}(1, 1:2) = lay.outline{1}((samples/4) - (index - 1), :);
lay.outline{1, 2}(2, 1:2) = [0, nose_scale];
lay.outline{1, 2}(3, 1:2) = lay.outline{1, 2}(1, 1:2) .* [-1 1];

% Add left ear
ear_half_width = 0.16 / max_image_size;
[~, index] = min(abs(theta-ear_half_width));
lay.outline{1, 3}(1, 1:2) = lay.outline{1}((samples / 2) - (index - 1), :); % top start
lay.outline{1, 3}(2, 1:2) = [-rad2 * 1.037, rad2 * 0.235];
lay.outline{1, 3}(3, 1:2) = [-rad2 * 1.07, rad2 * 0.2];
lay.outline{1, 3}(4, 1:2) = [-rad2 * 1.08, rad2 * 0.165];

lay.outline{1, 3}(5, 1:2) = [-rad2 * 1.05, rad2 * 0.02]; % left middle
lay.outline{1, 3}(6, 1:2) = [-rad2 * 1.09, -rad2 * 0.15];

lay.outline{1, 3}(7, 1:2) = [-rad2 * 1.075, -rad2 * 0.23];
lay.outline{1, 3}(7, 1:2) = [-rad2 * 1.065, -rad2 * 0.21];
lay.outline{1, 3}(8, 1:2) = lay.outline{1, 3}(1, 1:2) .* [1 -1]; % bottom end

% Add right ear
lay.outline{1, 4} = lay.outline{1, 3} .* [-1 1];

% Build equivalent 2D tangential basis for each position
if ~isempty(sensor_orientations)
    LPA_2D = [0.5 * outer_r / max_image_size, 0];
    
    R_2D_basis = zeros(length(lay.pos(:,1)), 3, 2);
    for sens_idx = 1:length(lay.pos(:,1))
        sens_pos_2D = lay.pos(sens_idx,1:2);
        [R_2D_basis(sens_idx,2,1:2), R_2D_basis(sens_idx,3,1:2), ~] = get_scalp_vectors(sens_pos_2D, LPA_2D, -LPA_2D);
    end
	
	% Fit real orientation to basis
	R_2D = zeros(num_sens, 2);
	for sens_idx = 1:num_sens    
    	cur_basis_frame = reshape(basis_frame(sens_idx,:),[3,3]);
    	cur_primary_tangential_basis_2D = squeeze(R_2D_basis(sens_idx,2,:));
    	cur_secondary_tangential_basis_2D = squeeze(R_2D_basis(sens_idx,3,:));
    	
    	% real 3D orientation vector
    	sens_ori = sensor_orientations(sens_idx,1:3)';
	
    	similarity_coeffs = cur_basis_frame.' * sens_ori;
    	c_radial  = similarity_coeffs(1);
    	c_primary_tan = similarity_coeffs(2);
    	c_secondary_tan = similarity_coeffs(3);
	
    	sens_ori_2D = c_primary_tan * cur_primary_tangential_basis_2D + c_secondary_tan * cur_secondary_tangential_basis_2D;
    	radial_weight = 1 - abs(c_radial) / norm(sens_ori);
    	radial_weight = max(radial_weight, 0);
    	
    	R_2D(sens_idx,1:2) = radial_weight * sens_ori_2D;
	end
	lay.ori = R_2D;
end
    

end


%==========================================================================
% - G E T   S U R F A C E   P O I N T S   A B O U T   P L A N E
%==========================================================================
function [intersect_points] = get_surface_points_about_plane(slice_vector, slice_origin, base_vector, base_origin, surface_points, edges, start_pos, plot, dec_point)
% Extract intersection points of two planes with a 3D surface and orders them
% to form lines across the surface.
%
% Syntax:
%   [intersect_points] = get_surface_points_about_plane(slice_vector, slice_origin, 
%                                                       base_vector, base_origin, 
%                                                       surface_points, edges, 
%                                                       start_pos, plot)
%
% Input:
%   slice_vector     - The normal vector of the slice plane.
%   slice_origin     - A point on the slice plane.
%   base_vector      - The normal vector of the base plane.
%   base_origin      - A point on the base plane.
%   surface_points  - An n x 3 matrix of points on a 3D surface.
%   edges            - An array defining edges between points on the surface.
%   start_pos        - The index of the starting point on the surface.
%   plot             - true/false whether to plot the lines formed by the 
%                      intersection points.
%
% Output:
%   intersect_points - An ordered set of intersection points between the two planes
%                     and the 3D surface, forming lines across the surface.
%
% Description:
%   This function calculates the intersection points between two planes and a 3D
%   surface, such as a slice plane and a base plane intersecting the surface.
%   It orders these intersection points to define lines across the surface.

% Normalise vectors
slice_vector = slice_vector / norm(slice_vector);
base_vector = base_vector / norm(base_vector);

% Add points where the base plane intersect with the surface to it.
% Find the intersect points
base_intersect_points = get_plane_surface_intersect(edges, surface_points, base_vector, base_origin, dec_point);

% Order the intersect points
base_intersect_points = order_points_on_plane(base_intersect_points, dec_point);

% Create new edges for the base intersect
start_edge_idx = 1 + length(surface_points(:, 1));
end_edge_idx = start_edge_idx + length(base_intersect_points(:, 1)) - 1;
base_edges = [(start_edge_idx:1:end_edge_idx)', [(start_edge_idx + 1):1:end_edge_idx, start_edge_idx]'];

% Add the points to the surface
surface_points = [surface_points; base_intersect_points];

% And concatenate edge lists
edges = [edges; base_edges];

% Get the slice interceptPoints, above the base
% Find the intersect points
slice_intersect_points = get_plane_surface_intersect(edges, surface_points, slice_vector, slice_origin, dec_point);

% Order the slice intersect
slice_intersect_points = order_points_on_plane(slice_intersect_points, dec_point);

% Get the two points that crossed the base plane
slice_edges = [(1:length(slice_intersect_points) - 1)', (2:length(slice_intersect_points))'; length(slice_intersect_points), 1];
slice_end_points = get_plane_surface_intersect(slice_edges, slice_intersect_points, base_vector, base_origin, dec_point);

% Sometimes there are two points which are very close. 
if length(slice_end_points(:, 1)) > 2
    while length(slice_end_points(:, 1)) > 2
        dec_point = dec_point - 1;
        slice_end_points = unique(round(slice_end_points, dec_point, 'decimals'), 'rows', 'stable');
    end
end

if (length(slice_end_points(:, 1)) < 2)
    error('Less than two points crossing. Something is wrong if this is happening.')
end

% Remove points below the base plane
[slice_intersect_points, ~] = remove_points_below_plane(slice_intersect_points, base_vector, base_origin);

% Find which endpoint is closest to the start point
if norm(start_pos - slice_end_points(1, :)) > norm(start_pos - slice_end_points(2, :))
    % Check the endpoints were not removed
    slice_intersect_points(ismember(slice_intersect_points, slice_end_points(1, :), 'rows'), :) = [];
    slice_intersect_points = [slice_end_points(1, :); slice_intersect_points];
    slice_intersect_points(ismember(slice_intersect_points, slice_end_points(2, :), 'rows'), :) = [];
    slice_intersect_points = [slice_end_points(2, :); slice_intersect_points];
else
    % Check the endpoints were not removed
    slice_intersect_points(ismember(slice_intersect_points, slice_end_points(2, :), 'rows'), :) = [];
    slice_intersect_points = [slice_end_points(2, :); slice_intersect_points];
    slice_intersect_points(ismember(slice_intersect_points, slice_end_points(1, :), 'rows'), :) = [];
    slice_intersect_points = [slice_end_points(1, :); slice_intersect_points];
end

% Order the slice points to make lines
slice_intersect_points = order_points_on_plane(slice_intersect_points, dec_point);
[slice_intersect_points] = reindex_line_ends(slice_intersect_points, slice_end_points);

intersect_points = slice_intersect_points;

% Plot if specified
if plot
    plot_ordered_points(intersect_points);
end

end


%==========================================================================
% - G E T   P L A N E   S U R F A C E   I N T E R S E C T 
%==========================================================================
function [intersect_pos] = get_plane_surface_intersect(edges, pos, plane_vec, plane_pos, dec_point)
% Calculate the intersection points between a plane and a set of line segments.
%
% Syntax:
%   [intersect_pos] = get_plane_surface_intersect(edges, pos, plane_vec, plane_pos)
%
% Input:
%   edges      - An array defining line segments by their vertex indices.
%   pos        - An array of vertex positions.
%   plane_vec  - The normal vector of the intersecting plane.
%   plane_pos  - A point on the intersecting plane.
%
% Output:
%   intersect_pos - An array containing unique intersection points between the
%                   plane and line segments.
%
% Description:
%   This function calculates the intersection points between a plane defined by
%   its normal vector and a set of line segments defined by their vertex indices.
%   It checks each line segment to see if it intersects with the plane and returns
%   the unique intersection points.
%
%   Note: The function assumes that the input line segments are defined as pairs
%   of vertex indices, and the plane is defined by a normal vector and a point
%   on the plane.

% Can't really preallocate
intersect_pos = zeros(0, 3);
count = 0;

% Need to go both directions to mirror edges
edges = [edges; flipud(edges')'];

% Go through each triangle
for i = 1:length(edges(:, 1))
    % Get the associated points
    p1 = pos(edges(i, 1), 1:3);
    p2 = pos(edges(i, 2), 1:3);

    % Define a line between them
    v = p2 - p1;

    % Calculate the dot product of the plane and vector
    dot_prod = dot(v, plane_vec);

    % If the plane lies between the points check where
    if dot_prod > 0 || abs(dot_prod) < 0.001
        % Calculate the parameter t for the line containing the segment
        t = dot(plane_vec, plane_pos - p1) / dot(plane_vec, v);

        % Calculate the intersection point
        p = p1 + t * v;

        % Check if the intersection point lies on the line segment
        t = round(dot(p - p1, v) / dot(v, v), dec_point, 'decimal');

        if t >= 0 && t <= 1
            % Store that point
            count = count + 1;
            intersect_pos(count, 1:3) = p;
        end
    end

end

intersect_pos = unique(round(intersect_pos, dec_point, 'decimals'), 'rows', 'stable');
end


%==========================================================================
% - F I N D   L I N E   T R I A N G L E   I N T E R S E C T I O N
%==========================================================================
function [intersection_point, vertices, tris] = find_line_triangle_intersection(direction, origin, vertices, tris, dec_point)
% Finds the intersection point between a line and a set of triangles in 3D space.
%
% Syntax:
%   [intersection_point, vertices, tris] = find_line_triangle_intersection(direction, origin, vertices, tris)
%
% Input:
%   direction - The direction vector of the line.
%   origin    - The origin point of the line.
%   vertices  - An array of vertices representing the triangle mesh.
%   tris      - An array defining the triangles using vertex indices.
%
% Output:
%   intersection_point - The intersection point between the line and a triangle (if any).
%   vertices           - Updated vertices array with the intersection point (if added).
%   tris               - Updated triangles array (if a triangle was intersected and split).
%
% Description:
%   This function calculates the intersection point between a line defined by its
%   direction and origin and a set of triangles defined by their vertices. It checks
%   if the line intersects with any of the triangles and, if so, returns the
%   intersection point and updates the triangle mesh if necessary.

num_triangles = size(tris, 1);
intersection_point = [];
intersect_index = [];

for i = 1:num_triangles
    % Get the vertices of the current triangle
    v1 = vertices(tris(i, 1), :);
    v2 = vertices(tris(i, 2), :);
    v3 = vertices(tris(i, 3), :);

    % Calculate the normal vector of the current triangle
    normal = cross(v2 - v1, v3 - v1);

   % Calculate the point of intersection with the infinite plane
    P0 = origin + dot(v1 - origin, normal) / dot(direction, normal) * direction;

    % Check if the intersection point lies within the triangle
    condition1 = dot(cross(v2 - v1, P0 - v1), normal) >= 0;
    condition2 = dot(cross(v3 - v2, P0 - v2), normal) >= 0;
    condition3 = dot(cross(v1 - v3, P0 - v3), normal) >= 0;
    condition4 = dot(P0 - origin, direction) > 0;

    if condition1 && condition2 && condition3 && condition4
        % Check there are not multiple
        if ~isempty(intersection_point)
            warning('Multiple triangles intersect with vector.');
        end

        % Calculate the intersection point within the plane
        intersection_point = P0;
        intersect_index = i;

    end
end

intersection_point = round(intersection_point, dec_point, 'decimals');

% If there is an intersection point, add it into the convhull
if ~isempty(intersection_point)
    % Remove the intersecting triangle
    originalTris = tris(intersect_index, :);
    tris(intersect_index, :) = [];

    % Make 3 new triangles to include the new point
    vertices(end + 1, :) = intersection_point;

    tris(end + 1, :) = [originalTris(1), originalTris(2), length(vertices(:, 1))];
    tris(end + 1, :) = [originalTris(2), originalTris(3), length(vertices(:, 1))];
    tris(end + 1, :) = [originalTris(1), originalTris(3), length(vertices(:, 1))];
end
end


%==========================================================================
% - O R D E R   P O I N T S   O N   A   P L A N E
%==========================================================================
function ordered_points = order_points_on_plane(points, dec_point)
    % Inputs:
    %   - points: nx3 matrix representing the coordinates of the points
    % Output:
    %   - ordered_points: nx3 matrix with the points sorted in a circular
    %        order around the origin, corrected for plane orientation
    points = unique(round(points, dec_point, 'decimals'), 'rows', 'stable');
    
    % If origin and normal are not provided, calculate them from the points
    origin = mean(points, 1);
    [~, ~, V] = svd(points - origin);
    normal = V(:, end)';

    % Calculate the vector from the origin to each point
    vectors = points - origin;

    % Calculate the angle between the vectors and a reference vector (e.g., [1, 0, 0])
    reference_vector = [1, 0, 0];
    cos_angles = dot(vectors, repmat(reference_vector, size(points, 1), 1), 2) ./ (vecnorm(vectors, 2, 2) * norm(reference_vector));
    angles = acosd(cos_angles);

    % Calculate the cross product between the vectors and the plane's normal vector
    cross_products = cross(vectors, repmat(normal, size(points, 1), 1));

    % Determine the correct sign of the angles using the cross product
    angle_signs = sign(dot(cross_products, repmat(reference_vector, size(points, 1), 1), 2));

    % Adjust the angles based on the sign of the angles
    angles = angles .* angle_signs;

    % Convert negative angles to positive angles in the range [0, 360)
    angles(angles < 0) = angles(angles < 0) + 360;

    % Sort the points based on their angles
    [~, sorted_indices] = sort(angles);
    ordered_points = points(sorted_indices, :);

    % Reindex to the original starting position
    first_point_index = find(ismember(ordered_points, points(1, :), 'rows'));
    if first_point_index > 1
        ordered_points = [ordered_points(first_point_index:end, :); ordered_points(1:first_point_index - 1, :)];
    end
end


%==========================================================================
% - R E M O V E   P O I N T S   B E L O W   P L A N E
%==========================================================================
function [pos, removed_pos] = remove_points_below_plane(pos, plane_vec, plane_origin)
% Remove any points below the base
rm_idx = zeros(length(pos(:, 1)), 1);
for i = 1:length(pos(:, 1))
    distance = dot(pos(i, 1:3) - plane_origin, plane_vec);
    
    if distance >= 0
        rm_idx(i, 1) = 1;
    end
end
removed_pos = pos(rm_idx == false, :);
pos(rm_idx == false, :) = [];

end

%==========================================================================
% - C U T   L I N E   E N D S 
%==========================================================================
function [shorter_lines] = cut_line_ends(lines, range, dec_point)
% Syntax:
%   [shorter_lines] = cut_line_ends(lines, range)
%
% Input:
%   lines      - An n x 3 matrix of positions in 3D space, forming a chain of lines.
%   range      - A 2-element array [start_frac, end_frac] specifying the
%                fractional positions along the total length of the line chain
%                to define the new endpoints for the trimmed line.
%
% Output:
%   shorter_lines - A line chain with its ends trimmed according to the specified
%                   range, resulting in a shorter line segment.
%
% Description:
%   This function trims the ends of a chain of 3D lines based on a specified
%   range given as fractional positions along the total length of the line chain.
%   It calculates the new endpoints for the trimmed line and returns the shorter
%   line segment.
%
%   Note: The function assumes that 'range' is specified in increasing order,
%   and 'lines' forms a continuous chain of lines in 3D space.

flipback = false;
if range(1) > range(2)
    flipback = true;
    range = flip(range);
end

% Get the total length of the line
total_length = 0;
for i = 1:size(lines, 1)-1
    total_length = total_length + norm(lines(i+1, :) - lines(i, :));
end

% Compute the Euclidean distance between each pair of adjacent points
d = sqrt(sum(diff(lines, [], 1).^2, 2));

% Compute the cumulative sum of the distances
cumd = [0; cumsum(d)];

% Normalise the cumulative distance
fracd = cumd / total_length;
fracd = round(fracd, dec_point, 'decimals');
% Find the new endpoints
idx = zeros(1, length(range));
end_points = zeros(length(range), 3);
for i = 1:length(range)
    % Find the index of the nearest point
    idx(i) = find(fracd >= range(i), 1, 'first');

    if idx(i) ~= 1
        % Compute the coordinates of the nearest point
        p1 = lines(idx(i)-1, :);
        p2 = lines(idx(i), :);
        d1 = fracd(idx(i)-1);
        d2 = fracd(idx(i));
        frac = (range(i) - d1) / (d2 - d1);
        end_points(i, :) = p1 + frac * (p2 - p1);
    else
        end_points(i, :) = lines(idx(i), :);
    end
end

end_points = round(end_points, dec_point, 'decimals');

% Trim
lines(idx(2):end, :) = [];
lines(1:idx(1)-1, :) = [];
lines = [end_points(1, :); lines; end_points(2, :)];
shorter_lines = lines;

if flipback
    shorter_lines = flip(shorter_lines);
end
end


%==========================================================================
% - F I N D   F R A C T I O N   A L O N G   L I N E 
%==========================================================================

function [frac, lines] = find_fraction_along_line(lines, point, dec_point)
% Calculates the fraction of the total length of a 3D line chain where a
% specified point lies.
%
% Syntax:
%   [frac, lines] = find_fraction_along_line(lines, point)
%
% Input:
%   lines  - An n x 3 matrix of positions in 3D space representing a chain of lines.
%   point  - A point in 3D space for which the fraction along the line chain
%            is to be determined.
%
% Output:
%   frac   - The fraction of the total length of the line chain where the specified
%            point lies. The fraction is a value between 0 and 1.
%   lines  - An updated line chain, which may include the specified point if it
%            was not already part of the chain.
%
%   Note: The function assumes that 'lines' forms a continuous chain of lines in
%   3D space and that 'point' is a 3D point in the same space.

lines = round(lines, dec_point, 'decimals');
end_points = [lines(1, :); lines(end, :)];
point = round(point, dec_point, 'decimals');

% Reorder the line, including this new point
if ~any(ismember(lines, point, 'rows'))
    lines(end + 1, :) = point;
    lines = unique(lines, 'rows', 'stable');
    lines = order_points_on_plane(lines, dec_point);
    lines = reindex_line_ends(lines, end_points);
end

% Compute the Euclidean distance between each pair of adjacent points
d = sqrt(sum(diff(lines, [], 1).^2, 2));

% Compute the cumulative sum of the distances
cumd = [0; cumsum(d)];

% Compute the total length of the line
total_length = cumd(end);

% Find the point in the new lines array
point_idx = find(ismember(lines, point, 'rows'), 1, 'first');

% If the point is no longer there, a rounding issue has occured.
if isempty(point_idx)
    % Find the closest point (it should be clear).
    [~, closestIdx] = min(vecnorm(lines - point,2,2));
    point = lines(closestIdx,:);

    % Find this replacement point in the new lines array
    point_idx = find(ismember(lines, point, 'rows'), 1, 'first');
end

dist = cumd(point_idx);

frac = dist / total_length;
end


%==========================================================================
% - P L O T   O R D E R E D   P O I N T S
%==========================================================================
function [] = plot_ordered_points(ordered_points)
% Takes some of the pain out of using plot3 and shows order as colour
% change. 

% Breakout
x = ordered_points(:, 1);
y = ordered_points(:, 2);
z = ordered_points(:, 3);

% Set a gradual change in colour to show directionality
num_points = size(ordered_points, 1);
cmap = jet(num_points);

% Plot
for i = 1:num_points - 1
    color = cmap(i, :); 
    plot3([x(i);x(i+1)], [y(i); y(i+1)], [z(i); z(i+1)], 'LineWidth', 3.5, 'Color', color);
    hold on;
end
end


%==========================================================================
% - R E - I N D E X   L I N E   E N D S
%==========================================================================
function [reordered_points] = reindex_line_ends(points, end_points)
% Reorders a set of circular points to start and end at specified points.
%
% Syntax:
%   [reordered_points] = reindex_line_ends(points, end_points)
%
% Input:
%   points         - An n x 3 matrix of 3D points forming a circular sequence.
%   end_points     - A pair of 3D points that represent the desired start and
%                    end points for reordering the circular sequence.
%
% Output:
%   reordered_points - The circular sequence of points reordered to start and
%                      end at the specified points.
%
% Description:
%   Note: The function assumes that the input points form a continuous circular
%   sequence in 3D space and that 'end_points' consists of two distinct points
%   within the set of 'points'.

% Find the points of interest
start_point_idx = find(ismember(points, end_points(1, :), 'rows'), inf, 'last');
end_point_idx = find(ismember(points, end_points(2, :), 'rows'), inf, 'first');

% Check if it needs flipping
flipped = false;
if start_point_idx > end_point_idx
    tmp = start_point_idx;
    start_point_idx = end_point_idx;
    end_point_idx = tmp;
    flipped = true;
end

% See which route between circular points is shortest
num_points = size(points, 1);
forward_distance = mod(end_point_idx - start_point_idx, num_points);
backward_distance = mod(start_point_idx - end_point_idx, num_points);

% Make correction as required
if forward_distance > backward_distance
    reordered_points = points(start_point_idx:end_point_idx, :);
else
    reordered_points = flip([points(end_point_idx:end, :); points(1:start_point_idx, :)]);
end

% Correct flipping
if flipped
    reordered_points = flip(reordered_points);
end
end


%==========================================================================
% - T R I S   T O   E D G E S
%==========================================================================
function [edges] = tris_to_edges(tris)
% Converts a set of triangles into a list of unique edges.
%
% Syntax:
%   [edges] = tris_to_edges(tris)
%
% Input:
%   tris  - An array representing triangles using vertex indices.
%
% Output:
%   edges - A matrix containing unique edges extracted from the input triangles.

% Prep combinations of tris
tri_comb = nchoosek([1, 2, 3], 2); 
tri_comb = [tri_comb; flipud(tri_comb')'];

edges = zeros(length(tris(:, 1)*length(tri_comb(:, 1))), 2);
count = 0;
for i = 1:length(tris(:, 1))
    % Then each side of it
    for j = 1:length(tri_comb(:, 1))
        count = count + 1;
        % Get the associated points
        edges(count, 1) = tris(i, tri_comb(j, 1));
        edges(count, 2) = tris(i, tri_comb(j, 2));
    end
end
end


%==========================================================================
% - L E N G T H   O F   L I N E S
%==========================================================================
function [total_length, d] = length_of_lines(lines)
% Useful method to get the cumulative length of a line defined by many
% points.
d = sqrt(sum(diff(lines, [], 1).^2, 2));

% Compute the cumulative sum of the distances
cumd = [0; cumsum(d)];

% Compute the total length of the line
total_length = cumd(end);

end

%==========================================================================
% - G E T   S C A L P   V E C T O R S
%==========================================================================
function [primary_tangent, secondary_tangent, arc_points] = get_scalp_vectors(sensor_position_2D, LPA, RPA)
% Makes orthogonal scalp vectors using tangential basis
A = LPA;
B = sensor_position_2D;
C = RPA;

% Check if colinear
area2 = abs(det([B - A; C - A]));
scale = norm(B - A) * norm(C - A);
tol   = 1e-6;

if area2 / scale < tol
    % Correct for left/right/middle
    base_dir = C - A;
    base_dir = base_dir / norm(base_dir);

    if dot(B - A, base_dir) < 0
        primary_tangent = -base_dir;
        secondary_tangent = [primary_tangent(2), -primary_tangent(1)];
    elseif dot(B - C, base_dir) > 0
        primary_tangent = -base_dir;
        secondary_tangent = [primary_tangent(2), -primary_tangent(1)];
    else
        primary_tangent = base_dir;
        secondary_tangent = [primary_tangent(2), -primary_tangent(1)];
    end

    arc_points = [linspace(A(1), C(1), 2)', ...
                  linspace(A(2), C(2), 2)'];
    return
end

% Circle through points
mid_AB = (A + B) / 2;
mid_BC = (B + C) / 2;

dir_AB = B - A;
dir_BC = C - B;

perp_AB = [-dir_AB(2), dir_AB(1)];
perp_BC = [-dir_BC(2), dir_BC(1)];

M = [perp_AB(:), -perp_BC(:)];
rhs = (mid_BC - mid_AB)';
params = M \ rhs;

centre = mid_AB + params(1) * perp_AB;
r = norm(A - centre);

% Find it on the arc
theta_A = atan2(A(2) - centre(2), A(1) - centre(1));
theta_B = atan2(B(2) - centre(2), B(1) - centre(1));
theta_C = atan2(C(2) - centre(2), C(1) - centre(1));

theta_A = mod(theta_A, 2*pi);
theta_B = mod(theta_B, 2*pi);
theta_C = mod(theta_C, 2*pi);

% Arc direction
ccw_AC = mod(theta_C - theta_A, 2*pi);
ccw_AB = mod(theta_B - theta_A, 2*pi);

if ccw_AB <= ccw_AC
    theta = linspace(theta_A, theta_A + ccw_AC, 200);
else
    theta = linspace(theta_A, theta_A - (2*pi - ccw_AC), 200);
end
  
arc_points = [centre(1) + r * cos(theta(:)), centre(2) + r * sin(theta(:))];

% Tangents
[~, idx] = min(vecnorm(arc_points - sensor_position_2D, 2, 2));

if idx > 1 && idx < size(arc_points, 1)
    primary_tangent = arc_points(idx + 1, :) - arc_points(idx - 1, :);
else
    primary_tangent = arc_points(min(idx + 1, end), :) ...
            - arc_points(max(idx - 1, 1), :);
end

primary_tangent = primary_tangent / norm(primary_tangent);
secondary_tangent = [primary_tangent(2), -primary_tangent(1)];
end
