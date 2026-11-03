function [positions, precalculated] = spm_inv_anatomical_layout(sensor_2d_positions, head_positions, fiducials, varargin)
% Inverse of spm_get_anatomical_layout
% FORMAT [sensor_positions] = spm_inv_anatomical_layout(sensor_positions, sensor_labels, head_positions, fiducials, plot_output)
%
% Input Parameters:
%     lay: A structure representing the anatomical sensor layout in the
%         FieldTrip style. It contains the following fields:
%         sensor_2d_pos: A nx2 matrix representing the 2D coordinates of the
%                        sensors in the layout. Each column contains the
%                        [x, y] coordinates of a sensor relative to the
%                        vertex (Cz), on the scale of
%                        spm_get_anatomical_layout (the 10-20
%                        circumference at radius 0.25).
%     sensor_labels:     An nx1 cell array containing labels for each sensor.
%     head_positions:    An mx3 matrix representing the coordinates of
%                        positions on the head.
%     fiducials:         Struct containing coordinates of fiducials.
%         fiducials.NAS = [1 x 3]
%         fiducials.LPA = [1 x 3]
%         fiducials.RPA = [1 x 3]
%         fiducials.INI = [1 x 3] - Optional
%                       Note, if fiducials.INI is not specified it will be
%                       estimated.
%     plot_output:      A logical value indicating whether to generate a
%                       plot of measurements taken in 3D, as well as the
%                       output layout.
%     unit:             A string, e.g. 'mm', will be estimated otherwise.
%     precalculated:    Struct with precalculated values for optimisation.
%
% Output:
%     sensor_positions: An nx3 matrix representing the 3D Cartesian
%                       coordinates (x, y, z) of n sensors.
%__________________________________________________________________________
%
% Further help:
%
% spm_inv_anatomical_layout is a function that defines the scalp surface
% according to the approximate polar grid used to place electrodes in the
% 10-20 system. The method also receives a 2D layout - i.e. a set of
% positions in 2D projected space. It applies the inverse of
% spm_get_anatomical_layout to find the 3D positions on the scalp for each
% 2D position.

% The method then measures the position of on-scalp sensors
% in relation to this polar grid (angle and eccentricity) and applies this
% to a standard 2D polar grid. For full details see Alexander et al (2025)
% https://doi.org/10.1111/ejn.70060.
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
%   Primary tangential set comes from making a slice from opposing polar grid
%     positions (e.g. T3-T4; FPz-Oz), normal to the origin and passing through
%     the sensor position. The primary tangential axis simply follows along
%     that line while orthogonal to the radial axis. The secondary tangential
%     axis is the cross product of the primary tangential and radial axes. 
%_________________________________________________________________________

% Nicholas Alexander
% Copyright (C) 2026 Wellcome Centre for Human Neuroimaging


%==========================================================================
% - S P M   G E T   L A Y O U T   3 D   P O S I T I O N S
%==========================================================================
%-TIDY up the inputs
%==========================================================================
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
    precalculated = varargin{3};
else
    precalculated = [];
end

% The head positions as given, before landmarks and lines are added, for
% the scalp normals
head_points = head_positions;

if isempty(precalculated)
	origin = mean([fiducials.LPA; fiducials.RPA], 1);
	
	
	% Estimate units
	if isempty(unit)
    	largest_dim = max(max(head_positions) - min(head_positions));
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


	% Radius of the scalp patch for the radial normal
	normal_radius = 0.015 * 10^unit_dec_adjust;

	precalculated = [];
	precalculated.normal_radius         = normal_radius;
	precalculated.origin_elec           = origin_elec;
	precalculated.origin = origin;
	precalculated.s_vec                 = s_vec;
	precalculated.dec_point             = dec_point;
	
	precalculated.head_positions        = head_positions;
	precalculated.tri                   = tri;
	
	precalculated.elec                  = elec;
	
	precalculated.ar_intersect_points   = ar_intersect_points;
	precalculated.rp_intersect_points   = rp_intersect_points;
	precalculated.pl_intersect_points   = pl_intersect_points;
	precalculated.la_intersect_points   = la_intersect_points;

end

% Use precalculated values if available
if ~isempty(precalculated)
	normal_radius         = precalculated.normal_radius;
	origin_elec           = precalculated.origin_elec;
	s_vec                 = precalculated.s_vec;
	dec_point             = precalculated.dec_point;
	origin = precalculated.origin;
	head_positions        = precalculated.head_positions;
	tri                   = precalculated.tri;
	
	elec                  = precalculated.elec;
	
	ar_intersect_points   = precalculated.ar_intersect_points;
	rp_intersect_points   = precalculated.rp_intersect_points;
	pl_intersect_points   = precalculated.pl_intersect_points;
	la_intersect_points   = precalculated.la_intersect_points;

end

%-REPRESENT the 2D layout positions in 3D based on the polar grid
%==========================================================================


% Prepare output
num_sens = length(sensor_2d_positions(:,1));
sensor_positions = nan(num_sens,3);
R = zeros(num_sens,9);
head_positions = round(head_positions,dec_point,'decimals');
% Loop through each sensor positions
for sens_idx = 1:num_sens
    % extract position
    x_dist = sensor_2d_positions(sens_idx,1);
    y_dist = sensor_2d_positions(sens_idx,2);
    
	if (x_dist == 0) && (y_dist == 0)
		sensor_positions(sens_idx,1:3) = elec.Cz;
		tmp_head_pos = head_positions;
	else
    	% Distance from centre as a fraction of the distance to the 10-20
    	% circumference, at radius 0.25 as in spm_get_anatomical_layout
    	r_dist   = norm([x_dist, y_dist]) / 0.25;
    	
    	% Point on circumference
    	x = x_dist / r_dist;
    	y = y_dist / r_dist;
	
    	% Which quadrant
    	theta = atan2(y, x);
    	theta = mod(pi/2 - theta, 2*pi);
    	
    	quadrant = floor(theta / (pi/2)) + 1;
    	
    	if quadrant == 1
        	c_intersect_points = ar_intersect_points;
    	elseif quadrant == 2
        	c_intersect_points = rp_intersect_points;
    	elseif quadrant == 3
        	c_intersect_points = pl_intersect_points;
    	elseif quadrant == 4
        	c_intersect_points = la_intersect_points;
    	end
	
    	% Find the point along the circumference we need
    	fraction = mod(theta, pi/2) / (pi/2);
	
    	circ_pos = find_point_along_line(c_intersect_points, fraction);
	
    	% Make a slice which passes through the position on the circumference
    	circ_vec = (circ_pos - origin_elec) / norm(circ_pos - origin_elec);
	
    	% Get the point sens_vec crosses the scalp
    	[circ_pos, tmp_head_pos, tmp_tri] = find_line_triangle_intersection(circ_vec, origin_elec, head_positions, tri, dec_point);
	
    	% Redefine tri as edges
    	[tmp_edges] = tris_to_edges(tmp_tri);
	
    	% Make a slice from Cz to the circumference position
    	% circ_sup_vec = mean([circ_vec;s_vec]) / norm(mean([circ_vec;s_vec]));
		if r_dist < 1
    		circ_sup_vec = mean([circ_vec;s_vec]) / norm(mean([circ_vec;s_vec]));
		else
    		circ_sup_vec = (circ_vec + 0.3 * s_vec) / norm(circ_vec + 0.3 * s_vec);
		end
	
    	tmpVec = cross(circ_vec, s_vec) / norm(cross(circ_vec, s_vec));
    	[circ_intersect_points] = get_surface_points_about_plane(tmpVec, origin_elec, circ_sup_vec, origin_elec, tmp_head_pos, tmp_edges, elec.Cz, false, dec_point);
    	[tmp_start, ~] = find_fraction_along_line(circ_intersect_points, elec.Cz, dec_point);
    	[tmp_end, ~] = find_fraction_along_line(circ_intersect_points, circ_pos, dec_point);
    	circ_intersect_points_orig = circ_intersect_points;
    	circ_intersect_points = cut_line_ends(circ_intersect_points, [tmp_start, tmp_end], dec_point);
	
		% Add the sensor intersection to the plot
    	% Find the point along that line which is equal to our distance from the
    	% centre
    	if r_dist < 1
        	cur_sens_pos = find_point_along_line(circ_intersect_points, r_dist);
    	else
        	[tmp_start, ~] = find_fraction_along_line(circ_intersect_points_orig, elec.Cz, dec_point);
        	if tmp_start > tmp_end
				circ_intersect_points_orig = cut_line_ends(circ_intersect_points_orig, [tmp_start, 0], dec_point);
			else
				circ_intersect_points_orig = cut_line_ends(circ_intersect_points_orig, [tmp_start, 1], dec_point);
        	end
			to_circ_length = length_of_lines(circ_intersect_points);
			sens_length = r_dist * to_circ_length;
	
			tot_length = length_of_lines(circ_intersect_points_orig);
			cur_sens_pos = find_point_along_line(circ_intersect_points_orig, sens_length / tot_length);
        	
    	end
    	% [tmp_start, ~] = find_fraction_along_line(circ_intersect_points_orig, elec.Cz, dec_point);
    	% [tmp_end, ~] = find_fraction_along_line(circ_intersect_points_orig, cur_sens_pos, dec_point);
    	% sens_intersect_points = cut_line_ends(circ_intersect_points_orig, [tmp_start, tmp_end], dec_point);
    	
    	if plot_output
	%         plot_ordered_points(sens_intersect_points)
        	scatter3(cur_sens_pos(1), cur_sens_pos(2), cur_sens_pos(3), 50, 'b','filled');
    	end
	
		cur_sens_pos =	round(cur_sens_pos,dec_point, 'decimals');
    	sensor_positions(sens_idx,1:3) = cur_sens_pos;
	end

	% Add orientations based on tangential basis set
	% Add it into the convhull
	[sens_pos, ~, ~, tri_idx] = closest_point_on_mesh(sensor_positions(sens_idx, :), head_positions, tri);
	sens_pos = round(sens_pos, dec_point, 'decimals');

	% Remove the intersecting triangle
	tmp_tri = tri;
	originalTris = tmp_tri(tri_idx, :);
	tmp_tri(tri_idx, :) = [];

	% Make 3 new triangles to include the new point
	tmp_head_pos(end, :) = sens_pos;
	tmp_tri(end + 1, :) = [originalTris(1), originalTris(2), length(tmp_head_pos(:, 1))]; %#ok<AGROW> 
	tmp_tri(end + 1, :) = [originalTris(2), originalTris(3), length(tmp_head_pos(:, 1))]; %#ok<AGROW> 
	tmp_tri(end + 1, :) = [originalTris(1), originalTris(3), length(tmp_head_pos(:, 1))]; %#ok<AGROW> 

    % Redefine tri as edges
    [tmp_edges] = tris_to_edges(tmp_tri);

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

	if plot_output
		plot_ordered_points(l_intersect_points)
	end

	% Find tan axis as line segment
	[~,cur_sens_idx] = ismember(sens_pos, l_intersect_points, 'rows');

	tmp = (l_intersect_points(cur_sens_idx-1,:) - l_intersect_points(cur_sens_idx,:)) / norm((l_intersect_points(cur_sens_idx-1,:) - l_intersect_points(cur_sens_idx,:)));
	tmp2 = (l_intersect_points(cur_sens_idx,:) - l_intersect_points(cur_sens_idx+1,:)) / norm((l_intersect_points(cur_sens_idx,:) - l_intersect_points(cur_sens_idx+1,:)));
	
	initial_primary_tan_vec = mean([tmp;tmp2]) / norm(mean([tmp;tmp2]));

	% Scalp normal at the sensor
	radial_vec = local_normal(sens_pos, head_points, origin_elec, normal_radius);

	% Orthonormal frame
	primary_tan_vec = initial_primary_tan_vec - dot(initial_primary_tan_vec, radial_vec) * radial_vec;
	primary_tan_vec = primary_tan_vec / norm(primary_tan_vec);
	secondary_tan_vec = cross(radial_vec, primary_tan_vec);

	R(sens_idx,1:9) = [radial_vec(:)', primary_tan_vec(:)', secondary_tan_vec(:)'];
	if plot_output
		quiver3(cur_sens_pos(1), cur_sens_pos(2), cur_sens_pos(3), radial_vec(1), radial_vec(2), radial_vec(3), 10, 'r', 'LineWidth', 1.5);
		quiver3(cur_sens_pos(1), cur_sens_pos(2), cur_sens_pos(3), primary_tan_vec(1), primary_tan_vec(2), primary_tan_vec(3), 10, 'g', 'LineWidth', 1.5);
		quiver3(cur_sens_pos(1), cur_sens_pos(2), cur_sens_pos(3), secondary_tan_vec(1), secondary_tan_vec(2), secondary_tan_vec(3), 10, 'b', 'LineWidth', 1.5);
	end
end

% Make positions table
n = size(sensor_positions, 1);
Px = repelem(sensor_positions(:,1), 3);
Py = repelem(sensor_positions(:,2), 3);
Pz = repelem(sensor_positions(:,3), 3);
Ox = reshape(R(:, [1 4 7]).', [], 1);
Oy = reshape(R(:, [2 5 8]).', [], 1);
Oz = reshape(R(:, [3 6 9]).', [], 1);
idx = repelem((1:n).', 3);
axis_labels = repmat(["X"; "Y"; "Z"], n, 1);
name = strcat(string(idx), "_", axis_labels);
positions = table(name, Px, Py, Pz, Ox, Oy, Oz);
end


%==========================================================================
% - L O C A L   N O R M A L
%==========================================================================
function n = local_normal(p, points, inside, radius)
% Unit normal of the surface through points near p
near = points(vecnorm(points - p, 2, 2) < radius, :);
if size(near, 1) < 6
	[~, order] = sort(vecnorm(points - p, 2, 2));
	near = points(order(1:min(12, end)), :);
end
[~, ~, V] = svd(near - mean(near, 1), 'econ');
n = V(:, 3)';
if dot(n, inside - p) < 0
	n = -n;
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
	tmp_dec_point = dec_point;
    while length(slice_end_points(:, 1)) > 2
        tmp_dec_point = tmp_dec_point - 1;
        slice_end_points = unique(round(slice_end_points, tmp_dec_point, 'decimals'), 'rows', 'stable');
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
% - F I N D   P O I N T   A L O N G   L I N E 
%==========================================================================

function point = find_point_along_line(lines, frac)
% Calculates the 3D point located at a given fraction along a line chain.
%
% Syntax:
%   point = find_point_along_line(lines, frac)
%
% Input:
%   lines - An n x 3 matrix of positions in 3D space representing a chain of lines
%   frac  - Fraction along the total length of the line chain (0–1)
%
% Output:
%   point - 1 x 3 vector giving the interpolated point on the line

% Clamp fraction safely
frac = max(0, min(1, frac));

% Segment lengths
d = sqrt(sum(diff(lines, 1, 1).^2, 2));

% Cumulative distances
cumd = [0; cumsum(d)];
total_length = cumd(end);

% Target distance along the line
target_dist = frac * total_length;

% Find the segment containing the target distance
seg_idx = find(cumd <= target_dist, 1, 'last');

% Handle exact end-point
if seg_idx == size(lines, 1)
    point = lines(end, :);
    return
end

% Fraction within this segment
seg_start_dist = cumd(seg_idx);
seg_length = d(seg_idx);
local_frac = (target_dist - seg_start_dist) / seg_length;

% Linear interpolation
p0 = lines(seg_idx, :);
p1 = lines(seg_idx + 1, :);

point = p0 + local_frac * (p1 - p0);

end

%==========================================================================
% - C L O S E S T   P O I N T   O N   M E S H 
%==========================================================================

function [closest_point, closest_normal, min_dist, tri_idx] = closest_point_on_mesh(point, verts, tris)
    % Initialise minimum distance
    min_dist = inf;
    closest_point = [NaN, NaN, NaN];
    closest_normal = [NaN, NaN, NaN];
    tri_idx = NaN;

    % Loop through each triangle
    for i = 1:size(tris, 1)
        tri = verts(tris(i, :), :); 

        % Get closest point on this triangle
        [p_closest, ~] = point_triangle_distance(tri, point);
        dist2 = sum((p_closest - point).^2);
        if dist2 < min_dist
            min_dist = dist2;
            closest_point = p_closest;
            tri_idx = i;
            v1 = tri(2, :) - tri(1, :);
            v2 = tri(3, :) - tri(1, :);
            n = cross(v1, v2);
            closest_normal = n / norm(n);
        end
    end

    % Return actual distance
    min_dist = sqrt(min_dist);
end


%==========================================================================
% - P O I N T   T R I A N G L E   D I S T A N C E 
%==========================================================================
function [closest, is_inside] = point_triangle_distance(tri, p)
    % vertices
    a = tri(1, :);
    b = tri(2, :);
    c = tri(3, :);

    % EdgeS
    ab = b - a;
    ac = c - a;
    ap = p - a;

    % Check if point in vertex region outside A
    d1 = dot(ab, ap);
    d2 = dot(ac, ap);
    if d1 <= 0 && d2 <= 0
        closest = a;
        is_inside = false;
        return;
    end

    % Outside B
    bp = p - b;
    d3 = dot(ab, bp);
    d4 = dot(ac, bp);
    if d3 >= 0 && d4 <= d3
        closest = b;
        is_inside = false;
        return;
    end

    % AB
    vc = d1*d4 - d3*d2;
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1 / (d1 - d3);
        closest = a + v * ab;
        is_inside = false;
        return;
    end

    % C
    cp = p - c;
    d5 = dot(ab, cp);
    d6 = dot(ac, cp);
    if d6 >= 0 && d5 <= d6
        closest = c;
        is_inside = false;
        return;
    end

    % AC
    vb = d5*d2 - d1*d6;
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2 / (d2 - d6);
        closest = a + w * ac;
        is_inside = false;
        return;
    end

    % BC
    va = d3*d6 - d5*d4;
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        closest = b + w * (c - b);
        is_inside = false;
        return;
    end

    % Point inside face region
    denom = 1 / (va + vb + vc);
    v = vb * denom;
    w = vc * denom;
    closest = a + ab * v + ac * w;
    is_inside = true;
end