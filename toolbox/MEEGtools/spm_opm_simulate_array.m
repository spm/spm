function [positions, lay, info] = spm_opm_simulate_array(S)
% Make an OPM sensor array and put it on a head, ready for spm_opm_sim
% FORMAT [positions, lay, info] = spm_opm_simulate_array(S)
%
% Arrays are laid out on the 2D anatomical layout of 
% spm_get_anatomical_layout and put on  a head with spm_inv_anatomical_layout.
% See Alexander et al (2025) https://doi.org/10.1111/ejn.70060 for details.
% An array made once, for example on the SPM canonical head, keeps its 
% structure, matched to the 10-20 landmarks, on every head it is put on.
%
% Without S.lay a layout is made. With a head (S.mesh) the layout is put on
% it, giving sensor positions and orientations.
%
%   S               - input structure
% Optional fields of S:
% HEAD
%   S.mesh          - Scalp: a gifti or its file name, a patch structure
%                     (.faces, .vertices), a triangulation, or an SPM mesh
%                     from spm_eeg_inv_mesh (its scalp and fiducials)
%                     Default: [] (layout only)
%   S.fid           - Fiducials on the scalp, .NAS, .LPA, .RPA and
%                     optionally .INI, in its coordinates
%                     Default: [] (from an SPM mesh, otherwise picked)
%   S.pick_inion    - Also pick the inion - Default: 0
%   S.projection_vertices - Scalp downsampled to this many vertices for
%                     the projection, 0 for all - Default: 3000
% LAYOUT
%   S.lay           - Layout to put on the head: a structure, or a file
%                     read by ft_prepare_layout (.mat holding layout)
%                     Default: [] (make one)
%   S.method        - 'equidistant': hexagonal lattice filling the mask,
%                     in columns parallel to the midline
%                     'picked': clicked on the scalp
%                     Default: 'equidistant'
%   S.n_sensors     - Number of sensors, equidistant - Default: 32
%   S.mask          - Mx2 outline the array sits inside (layout units), or
%                     'draw' to draw it on the scalp
%                     Default: [] (drawn when picking, otherwise the circle
%                     through the fiducials)
%   S.symmetry      - 'left' or 'right': symmetric about the midline (when
%                     picking, this side is picked and mirrored); 'none'
%                     Default: 'left'
%   S.mask_symmetry - Drawn mask: mirror this side - Default: S.symmetry
%   S.midline       - Allow sensors on the midline - Default: 1
%   S.midline_snap  - Picks this near the midline go on it (mm) - Default: 5
%   S.midline_margin - Without midline sensors, keep the others this far
%                     from it (mm), e.g. room for a helmet's hinge
%                     Default: 0 (half a spacing)
%   S.edge_margin   - Equidistant: keep sensors this far inside the mask
%                     (mm) - Default: 5
%   S.relax         - Equidistant: even the lattice out by centroidal
%                     Voronoi relaxation, on the scalp if given, otherwise
%                     in 2D. 'area': each sensor covers an equal share of
%                     the mask; 'edge': the outer row pushed out to the edge
%                     margin; 'none' - Default: 'none'
%   S.relax_iterations - Maximum iterations - Default: 200
%   S.relax_grid    - Spacing of the grid mapped onto the scalp once for
%                     the relaxation and the picking checks (layout units)
%                     Default: 0.05
% CHECKS (on the scalp when it was mapped, otherwise 2D at S.mm_per_unit)
%   S.mm_per_unit   - mm per layout unit on an adult head - Default: 550
%   S.min_spacing   - Warn if neighbours are closer (mm) - Default: 13
%   S.check_margin  - Helmet: warn if sensors are nearer the mask edge or
%                     midline (mm) - Default: 10
%   S.helmet        - Check for a split helmet - Default: 0
%   S.footprint_radius - Picking: circle drawn round each pick (mm)
%                     Default: S.min_spacing/2
%   S.footprint_lift - Picking: circles this far off the scalp (mm)
%                     Default: 0
% SENSORS
%   S.axes          - 'radial' (1), 'dual' (2) or 'triaxial' (3)
%                     Default: 'radial'
%   S.offset        - Scalp to sensitive volume (mm) - Default: 6.5
%   S.output_file   - Write the positions here (.tsv) - Default: []
%   S.name          - For titles - Default: 'array'
%   S.plot          - Plot the layout and the array - Default: 1
%
% Output:
%   positions       - Table, one row per channel: name (label_X radial,
%                     _Y and _Z tangential), Px, Py, Pz, Ox, Oy, Oz, with
%                     orientations into the head. For spm_opm_sim, write it
%                     with S.output_file and pass the file as S.positions.
%                     [] without a head.
%   lay             - Layout, FieldTrip style
%   info            - .fid, .label, .pos (sensors), .scalp_pos, .radial,
%                     .tan1, .tan2, .nearest_mm (on this head), .mask_3d,
%                     .split_3d (on the scalp)
%__________________________________________________________________________
%
% Example:
% Build template layout and positions on canonical head
% S = [];
% S.mesh = spm_eeg_inv_mesh([], 2);
% S.n_sensors = 64;
% S.relax = 'area';
% S.axes = 3;
% S.mask = 'draw';
% S.offset = 11;
% S.output_file = [output_dir,'/positions_canonical.tsv'];
% [positions_canonical, lay] = spm_opm_simulate_array(S);
% 
% % Apply layout to individual head geometry
% S = [];
% S.mesh = spm_eeg_inv_mesh(individual_mri, 2);
% S.lay = lay;
% S.offset = 11;
% S.output_file = [output_dir,'/positions_individual.tsv'];
% [positions_individual] = spm_opm_simulate_array(S);
% 
% % Run simulations
% S = [];
% S.positions = [output_dir,'/positions_canonical.tsv'];
% S.fname = [output_dir,'sim_opm_canonical'];
% S.fs = 1000;
% S.nSamples = 1000;
% S.lead = 1;
% [D_canonical, L_canonical] = spm_opm_sim(S);
% 
% S.positions = [output_dir,'/positions_individual.tsv'];
% S.fname = [output_dir,'sim_opm_individual'];
% S.sMRI = individual_mri;
% [D_individual, L_individual] = spm_opm_sim(S);
%__________________________________________________________________________

% Nicholas Alexander
% Copyright (C) 2026 Department of Imaging Neuroscience


spm('FnBanner', mfilename);

% Set default values
%--------------------------------------------------------------------------
Sd = [];
Sd.mesh = [];
Sd.fid = [];
Sd.pick_inion = false;
Sd.projection_vertices = 3000;
Sd.lay = [];
Sd.method = 'equidistant';
Sd.n_sensors = 32;
Sd.mask = [];
Sd.symmetry = 'left';
Sd.mask_symmetry = [];
Sd.midline = true;
Sd.midline_snap = 5;
Sd.midline_margin = 0;
Sd.edge_margin = 5;
Sd.relax = 'none';
Sd.relax_iterations = 200;
Sd.relax_grid = 0.05;
Sd.mm_per_unit = 550;
Sd.min_spacing = 13;
Sd.check_margin = 10;
Sd.helmet = false;
Sd.footprint_radius = [];
Sd.footprint_lift = 0;
Sd.axes = 'radial';
Sd.offset = 10;
Sd.output_file = [];
Sd.name = 'array';
Sd.plot = true;

if ~nargin || isempty(S)
    S = struct();
end
unknown = setdiff(fieldnames(S), fieldnames(Sd));
if ~isempty(unknown)
    warning('Unknown settings ignored: %s', strjoin(unknown, ', '));
end
for f = fieldnames(Sd)'
    if ~isfield(S, f{1}), S.(f{1}) = Sd.(f{1}); end
end
if isempty(S.mask_symmetry)
	S.mask_symmetry = S.symmetry;
end
if isempty(S.footprint_radius), S.footprint_radius = S.min_spacing / 2; end

switch S.axes
    case {'radial', 1}
        n_axes = 1;
    case {'dual', 2}
        n_axes = 2;
    case {'triaxial', 3}
        n_axes = 3;
    otherwise
        error('Unknown axes %s, use radial, dual or triaxial.', num2str(S.axes));
end
if ~any(strcmp(S.symmetry, {'left', 'right', 'none'}))
    error('Unknown symmetry %s, use left, right or none.', S.symmetry);
end
symmetric = ~strcmp(S.symmetry, 'none');
picking = isempty(S.lay) && strcmp(S.method, 'picked');
drawing = isempty(S.lay) && (strcmp(S.mask, 'draw') || (picking && isempty(S.mask)));
have_head = ~isempty(S.mesh);
if (picking || drawing) && ~have_head
    error('Picking sensors or drawing the mask needs S.mesh.');
end

% Head
%--------------------------------------------------------------------------
% The projection slows with the vertex count, so it runs on a downsampled
% scalp and points are snapped back onto the full one
if have_head
    [M, fid] = read_mesh(S.mesh);
    if ~isempty(S.fid), fid = S.fid; end
    tri = triangulation(M.faces, M.vertices);
    M_proj = M;
    if S.projection_vertices > 0 && size(M.vertices,1) > S.projection_vertices
        M_proj = spm_mesh_reduce(M, 2 * S.projection_vertices);
    end
    head_points = double(M_proj.vertices);
    if isempty(fid)
        % Click a fiducial's button, then its place on the scalp; close
        % the window when all are placed
        fid_names = {'LPA', 'RPA', 'NAS'};
        if S.pick_inion, fid_names{end+1} = 'INI'; end
        fprintf('Pick the fiducials for %s: %s\n', S.name, strjoin(fid_names, ', '));
        P = spm_mesh_select(M, fid_names);
        if any(isnan(P(:)))
            error('Fiducial picking was closed before all fiducials were picked.');
        end
        fid = cell2struct(num2cell(P', 2), fid_names, 1);
    end
    missing = setdiff({'NAS', 'LPA', 'RPA'}, fieldnames(fid));
    if ~isempty(missing)
        error('Fiducials missing: %s.', strjoin(missing, ', '));
    end
end

% Layout to scalp map, made when picking or relaxing on the head
to_head = [];
cache   = [];


% Layout
%--------------------------------------------------------------------------
if ~isempty(S.lay)
    lay = read_layout(S.lay);
    mask = [];
    if isfield(lay, 'mask') && ~isempty(lay.mask)
        mask = lay.mask;
        if iscell(mask), mask = mask{1}; end
        mask = open_polygon(mask);
    end
    % The split must run past the mask at both ends
    if ~isempty(mask)
        if ~isfield(lay, 'split') || isempty(lay.split)
            lay.split = midline_split(mask);
        elseif any(inpolygon(lay.split([1 end],1), lay.split([1 end],2), mask(:,1), mask(:,2)))
            warning('The layout split does not run past the mask, using a straight one back to front.');
            lay.split = midline_split(mask);
        end
    end
    if ~isfield(lay, 'neighbours') || isempty(lay.neighbours)
        lay.neighbours = layout_neighbours(lay.pos);
    end
else
    lay = struct();

    % Mask
    %----------------------------------------------------------------------
    mask = S.mask;
    if iscell(mask)
        mask = mask{1};
    end
    if drawing
        % Drawn on the scalp and projected into 2D
        outline_picked = pick_points(tri, fid, 'Draw the edge of the array around the ears, face and neck');
        outline_lay = spm_get_anatomical_layout(outline_picked, number_labels(size(outline_picked,1)), ...
            head_points, fid, false, 'mm', [], inf);
        mask = symmetric_outline(outline_lay.pos, S.mask_symmetry);
        lay.outline_picked = outline_picked;
        lay.outline_picked_2d = outline_lay.pos;
    elseif isempty(mask)
        % Circle through the fiducials: the 10-20 circumference (radius
        % 0.25) is 40% of the way from Cz to them, so they are at 50%
        theta = linspace(0, 2*pi, 101)';
        mask = 0.25 * 0.5/0.4 * [sin(theta(1:end-1)), cos(theta(1:end-1))];
    end
    mask = open_polygon(mask);
    split_2d = midline_split(mask);

    % Sensors
    %----------------------------------------------------------------------
    switch S.method
        case 'picked'
            % Mask and midline shown on the scalp, a footprint circle round
            % each pick, and picks too near the edge, the midline or
            % another sensor (or its mirror) flagged in red
            [to_head, cache] = head_map(mask, head_points, fid, cache, M, S.relax_grid);
            H = layout_to_head([mask; split_2d], head_points, fid, cache, M);
            fine_2d = layout_grid(mask, 0.005);
            chk = [];
            chk.fine_2d = fine_2d;
            chk.fine_3d = to_head(fine_2d);
            chk.mask = mask;
            chk.edge_3d = to_head(resample_polygon(mask, 1000));
            chk.midline_3d = to_head(midline_points(mask, 500));
            chk.to_head = to_head;
            chk.edge_margin = S.check_margin;
            chk.midline_margin = max(S.check_margin, S.midline_margin);
            chk.min_spacing = S.min_spacing;
            chk.side = S.symmetry;
            chk.allow_midline = S.midline;
            chk.midline_snap = S.midline_snap;

            prompt = 'Pick the sensor positions';
            if symmetric
                prompt = sprintf('Pick the sensor positions on the %s side (mirrored to the other)', S.symmetry);
            end
            if S.midline
                prompt = sprintf('%s, within %g mm of the midline goes on it', prompt, S.midline_snap);
            end
            P = [];
            P.outline = H.pos(1:size(mask,1),:);
            P.split = H.pos(size(mask,1)+1:end,:);
            P.radius = S.footprint_radius;
            P.lift = S.footprint_lift;
            P.check = @(p) check_picks(p, chk);
            picked = pick_points(tri, fid, prompt, false, P);
            if isempty(picked)
                error('No sensors picked.');
            end
            picked_lay = spm_get_anatomical_layout(picked, number_labels(size(picked,1)), ...
                head_points, fid, false, 'mm', [], inf);
            pos = picked_lay.pos;
            if S.midline
                pos(abs(pos(:,1)) < S.midline_snap / S.mm_per_unit, 1) = 0;
            end
            [lay.pos, lay.label] = picked_labels(pos, S.symmetry, S.midline);
            lay.sensors_picked = picked;
            lay.sensors_picked_2d = picked_lay.pos;

        case 'equidistant'
            [lay.pos, lay.label, lay.spacing] = equidistant_points(mask, S.n_sensors, symmetric, S.midline, ...
                S.edge_margin / S.mm_per_unit, S.midline_margin / S.mm_per_unit);

            if ~strcmp(S.relax, 'none')
                % Measured on the head if there is one
                if have_head
                    [to_head, cache] = head_map(mask, head_points, fid, cache, M, S.relax_grid);
                    relax_map = to_head;
                else
                    relax_map = @(p) [p, zeros(size(p,1),1)] * S.mm_per_unit;
                end
                R = [];
                R.mode = S.relax;
                R.edge_margin = S.edge_margin;
                R.midline_margin = S.midline_margin;
                R.iterations = S.relax_iterations;
                R.pairs = [find(startsWith(lay.label, 'L')), find(startsWith(lay.label, 'R'))];
                R.mid = find(startsWith(lay.label, 'M'));
                [lay.pos, lay.relax] = relax_layout(lay.pos, mask, relax_map, R);
                lay.relax.mode = S.relax;
            end

        otherwise
            error('Unknown layout method %s, use equidistant or picked.', S.method);
    end

    lay.neighbours = layout_neighbours(lay.pos);
    
    spacing = min(nearest_other(lay.pos));
    lay.width    = ones(size(lay.pos,1), 1) * spacing * 0.8;
    lay.height   = ones(size(lay.pos,1), 1) * spacing * 0.6;
    lay.mask     = {mask};
    lay.outline  = {mask};
    lay.split    = split_2d;
    lay.method   = S.method;
    lay.symmetry = S.symmetry;
    lay.midline  = S.midline;
    lay.name     = S.name;
    lay.created  = char(datetime('now'));

    % Check
    %----------------------------------------------------------------------
    % On the head where it was mapped, otherwise a 2D heuristic. 
    if isempty(to_head)
        m = measure(lay.pos, mask, @(p) [p, zeros(size(p,1),1)] * S.mm_per_unit);
    else
        m = measure(lay.pos, mask, to_head);
    end
    lay.metrics = m;
    if any(m.edge_mm < 0)
        warning('%d sensors are outside the mask: %s', sum(m.edge_mm < 0), strjoin(lay.label(m.edge_mm < 0), ', '));
    end
    close_ones = m.nearest_mm < S.min_spacing;
    if any(close_ones)
        warning('%d sensors are closer than %g mm to a neighbour: %s', sum(close_ones), S.min_spacing, strjoin(lay.label(close_ones), ', '));
    end
    if S.helmet
        near_edge = m.edge_mm >= 0 & m.edge_mm < S.check_margin;
        if any(near_edge)
            warning('%d sensors are within %g mm of the edge: %s', sum(near_edge), S.check_margin, strjoin(lay.label(near_edge), ', '));
        end
    end
    midline_limit = S.midline_margin;
    if S.helmet
        midline_limit = max(midline_limit, S.check_margin);
    end
    near_midline = m.midline_mm < midline_limit;
    if any(near_midline)
        warning('%d sensors are within %g mm of the midline, where the split and hinge run: %s', ...
            sum(near_midline), midline_limit, strjoin(lay.label(near_midline), ', '));
    end

    %-Plot
    %----------------------------------------------------------------------
    if S.plot
        figure;
        hold on
        axis equal
        grid on
        plot_layout_reference();
        if isfield(lay, 'outline_picked_2d')
            plot(lay.outline_picked_2d([1:end 1],1), lay.outline_picked_2d([1:end 1],2), 'r:', 'linewidth', 1);
        end
        plot(mask([1:end 1],1), mask([1:end 1],2), 'r-', 'linewidth', 1.5);
        plot(split_2d(:,1), split_2d(:,2), 'g:', 'linewidth', 1);
        plot_neighbours(lay.pos, lay.neighbours);
        flagged = close_ones | m.edge_mm < 0;
        scatter(lay.pos(~flagged,1), lay.pos(~flagged,2), 30, 'b', 'filled');
        scatter(lay.pos(flagged,1), lay.pos(flagged,2), 30, 'm', 'filled');
        text(lay.pos(:,1), lay.pos(:,2), lay.label, 'fontsize', 7, 'verticalalignment', 'bottom');
        title(sprintf('%s: %d sensors, %s', S.name, size(lay.pos,1), S.method), ...
            'interpreter', 'none');
    end
end
lay.label = cellstr(lay.label(:));


% Put the array on the head
%--------------------------------------------------------------------------
positions = [];
info = [];
info.label = lay.label;
if ~have_head
    if ~isempty(S.output_file)
        warning('No head (S.mesh) to put the array on, nothing written.');
    end
    return
end

% The projection runs on the downsampled scalp, so positions are snapped
% back onto the full scalp, then moved out along the radial axis (the
% scalp normal) by the offset
[H, cache] = layout_to_head(lay.pos, head_points, fid, cache, M);
pos = H.pos - S.offset * H.radial;
n = size(pos,1);
axis_names = ["X", "Y", "Z"];
axis_ori = cat(3, H.radial, H.tan1, H.tan2);
axis_ori = axis_ori(:,:,1:n_axes);

% One row per channel, sensor by sensor
positions = table();
positions.name = cellstr(reshape(string(lay.label)' + "_" + axis_names(1:n_axes)', [], 1));
positions.Px = repelem(pos(:,1), n_axes);
positions.Py = repelem(pos(:,2), n_axes);
positions.Pz = repelem(pos(:,3), n_axes);
ori = reshape(permute(axis_ori, [3 1 2]), [], 3);
positions.Ox = ori(:,1);
positions.Oy = ori(:,2);
positions.Oz = ori(:,3);

nearest_mm = nearest_other(pos);
fprintf('%s on the head: %d sensors, %d channels, spacing about %.1f mm (%.1f - %.1f)\n', ...
    S.name, n, height(positions), median(nearest_mm), min(nearest_mm), max(nearest_mm));
% Far outside the 10-20 circle the inverse layout runs out of scalp and
% several 2D points land on the same place
stacked = nearest_mm < 1;
if any(stacked)
    warning('%d sensors landed on top of another (beyond the reach of the projection): %s', ...
        sum(stacked), strjoin(lay.label(stacked), ', '));
end
too_close = nearest_mm < S.min_spacing & ~stacked;
if any(too_close)
    warning('%d sensors are closer than %g mm on this head: %s', sum(too_close), S.min_spacing, ...
        strjoin(lay.label(too_close), ', '));
end

info.fid = fid;
info.pos = pos;
info.scalp_pos = H.pos;
info.radial = H.radial;
info.tan1 = H.tan1;
info.tan2 = H.tan2;
info.nearest_mm = nearest_mm;
info.mask_3d = [];
info.split_3d = [];
if ~isempty(mask)
    H = layout_to_head([mask; lay.split], head_points, fid, cache, M);
    info.mask_3d = H.pos(1:size(mask,1),:);
    info.split_3d = H.pos(size(mask,1)+1:end,:);
end

if ~isempty(S.output_file)
    writetable(positions, S.output_file, 'FileType', 'text', 'delimiter', '\t');
    fprintf('Written to %s\n', S.output_file);
end

if S.plot
    figure;
    hold on
    axis equal
    trisurf(tri, 'EdgeAlpha', 0, 'FaceColor', [0.9 0.75 0.65]);
    if ~isempty(info.mask_3d)
        plot3(info.mask_3d([1:end 1],1), info.mask_3d([1:end 1],2), info.mask_3d([1:end 1],3), 'r-', 'linewidth', 1.5);
    end
    scatter3(pos(:,1), pos(:,2), pos(:,3), 20, 'k', 'filled');
    axis_colours = 'rgb';
    for a = 1:n_axes
        quiver3(pos(:,1), pos(:,2), pos(:,3), axis_ori(:,1,a)*10, axis_ori(:,2,a)*10, axis_ori(:,3,a)*10, 0, axis_colours(a));
    end
    text(pos(:,1), pos(:,2), pos(:,3), lay.label, 'fontsize', 7);
    view(3);
    follow_light(gca);
    title(sprintf('%s: %d sensors, %s (red radial, green and blue tangential)', S.name, n, num2str(S.axes)), ...
        'interpreter', 'none');
end
end


%==========================================================================
% Head and layout input
%==========================================================================
function [M, fid] = read_mesh(mesh)
% Scalp as a patch structure, and the fiducials of an SPM mesh
fid = [];
if ischar(mesh), mesh = gifti(mesh); end
if isstruct(mesh) && isfield(mesh, 'tess_scalp')
    labels = upper(mesh.fid.fid.label);
    for k = 1:numel(labels)
        if any(strcmp(labels{k}, {'NAS', 'LPA', 'RPA'}))
            fid.(labels{k}) = mesh.fid.fid.pnt(k,:);
        end
    end
    mesh = gifti(mesh.tess_scalp);
end
if isa(mesh, 'triangulation')
    M = struct('faces', mesh.ConnectivityList, 'vertices', mesh.Points);
elseif isa(mesh, 'gifti') || isfield(mesh, 'faces')
    M = struct('faces', double(mesh.faces), 'vertices', double(mesh.vertices));
elseif isfield(mesh, 'tri')
    M = struct('faces', double(mesh.tri), 'vertices', double(mesh.pos));
else
    error('S.mesh should be a gifti, a patch structure, a triangulation or an SPM mesh.');
end

end

function lay = read_layout(given)
if ischar(given) && endsWith(given, '.mat')
    tmp = load(given, 'layout');
    given = tmp.layout;
end
cfg = [];
cfg.layout    = given;
cfg.skipscale = 'yes';
cfg.skipcomnt = 'yes';
cfg.outline   = 'no';
cfg.mask      = 'no';
cfg.feedback  = 'no';
lay = ft_prepare_layout(cfg);
if isstruct(given)
    for f = {'mask', 'outline', 'split', 'neighbours', 'name'}
        if isfield(given, f{1}) && (~isfield(lay, f{1}) || isempty(lay.(f{1})))
            lay.(f{1}) = given.(f{1});
        end
    end
end
end


%==========================================================================
% Projection
%==========================================================================
function [H, cache] = layout_to_head(pos2d, head_points, fid, cache, M)
% 2D layout positions on the head with spm_inv_anatomical_layout. The
% projection works on the convex hull of the head points, so in hollows
% (eye sockets, under the nose). To address this points are moved back onto
% the original scalp mesh.
if isempty(cache)
    [T, cache] = spm_inv_anatomical_layout(pos2d, head_points, fid, false, 'mm');
else
    T = spm_inv_anatomical_layout(pos2d, head_points, fid, false, 'mm', cache);
end
H.pos    = [T.Px(1:3:end), T.Py(1:3:end), T.Pz(1:3:end)];
H.radial = [T.Ox(1:3:end), T.Oy(1:3:end), T.Oz(1:3:end)];
H.tan1   = [T.Ox(2:3:end), T.Oy(2:3:end), T.Oz(2:3:end)];
H.tan2   = [T.Ox(3:3:end), T.Oy(3:3:end), T.Oz(3:3:end)];
H.pos    = snap_to_mesh(H.pos, cache.origin_elec - H.pos, M);
end

function pts = snap_to_mesh(pts, dirs, M)
% Move each point along its direction, into the head, to where the line
% first crosses the mesh. Points start outside or on the scalp, bar a
% little slack for the downsampled scalp they came from. Points whose line
% misses within max_shift mm are left as they are.
slack = 1;
max_shift = 30;
V = [M.vertices(M.faces(:,1),:)'; M.vertices(M.faces(:,2),:)'; M.vertices(M.faces(:,3),:)'];
missed = false(size(pts,1),1);
for k = 1:size(pts,1)
    d = dirs(k,:) / norm(dirs(k,:));
    [~, ~, t] = spm_mesh_ray_intersect(V, struct('orig', (pts(k,:) - slack * d)', 'vec', d'));
    t = t(t <= max_shift + slack);
    if isempty(t)
        missed(k) = true;
    else
        pts(k,:) = pts(k,:) + (min(t) - slack) * d;
    end
end
if any(missed)
    warning('%d points did not meet the scalp within %g mm and were not snapped.', sum(missed), max_shift);
end
end

function [to_head, cache] = head_map(mask, head_points, fid, cache, M, step)
% Layout to scalp map: a grid over the mask is put on the head once and
% interpolated, so later steps (checks while picking, the relaxation) need
% no more projections
grid_2d = layout_grid(mask, step);
fprintf('Mapping %d grid points onto the head\n', size(grid_2d,1));
[H, cache] = layout_to_head(grid_2d, head_points, fid, cache, M);
F = cell(1,3);
for c = 1:3
    F{c} = scatteredInterpolant(grid_2d(:,1), grid_2d(:,2), H.pos(:,c), 'natural', 'nearest');
end
to_head = @(p) [F{1}(p(:,1), p(:,2)), F{2}(p(:,1), p(:,2)), F{3}(p(:,1), p(:,2))];
end

function g = layout_grid(mask, step)
% Square grid over the mask and up to 0.1 beyond its edge
r = max(vecnorm(mask, 2, 2)) + 0.1;
[gx, gy] = meshgrid(-r:step:r);
g = [gx(:), gy(:)];
[~, d] = nearest_point(resample_polygon(mask, 1000), g);
g = g(inpolygon(g(:,1), g(:,2), mask(:,1), mask(:,2)) | d < 0.1, :);
end


%==========================================================================
% Layouts
%==========================================================================
function [pos, label, spacing] = equidistant_points(mask, n, symmetric, midline, edge_gap, midline_gap)
% Hexagonal lattice in columns parallel to the midline, at the largest
% spacing that fits at least n points edge_gap inside the mask, then
% trimmed to n from nearest the mask edge. Symmetric layouts are built on
% the left and mirrored, so they are exactly symmetric whatever the mask.
% With midline points the lattice has a column on the midline. Without,
% the nearest columns are half a spacing either side, so neighbours across
% the midline are one spacing apart, or midline_gap from it if that is
% more.
if symmetric && ~midline && mod(n, 2) == 1
    error('A symmetric layout without midline sensors needs an even number of sensors.');
end
r_max = max(vecnorm(mask, 2, 2));
boundary = resample_polygon(mask, 2000);
tol = 1e-9;

pos = [];
for spacing = 0.3 * 0.998.^(0:3000)
    h = spacing * sqrt(3) / 2;
    n_col = ceil(r_max / h) + 2;
    n_row = ceil(r_max / spacing) + 2;
    rows = (-n_row:n_row)' * spacing;

    if symmetric && ~midline
        col_x = -(max(spacing/2, midline_gap) + (0:n_col) * h);
        col_shift = mod(0:n_col, 2) * spacing / 2;
    elseif symmetric
        col_x = -(0:n_col) * h;
        col_shift = mod(0:n_col, 2) * spacing / 2;
    else
        col_x = (-n_col:n_col) * h;
        col_shift = mod(-n_col:n_col, 2) * spacing / 2;
    end
    cand = [repelem(col_x(:), numel(rows)), reshape(rows + col_shift, [], 1)];
    cand = cand(inpolygon(cand(:,1), cand(:,2), mask(:,1), mask(:,2)),:);
    [~, d] = nearest_point(boundary, cand);
    cand = cand(d >= edge_gap,:);
    d = d(d >= edge_gap);

    if symmetric
        on_mid = abs(cand(:,1)) < tol;
        left = cand(~on_mid,:);
        d_left = d(~on_mid);
        mid = [zeros(sum(on_mid),1), cand(on_mid,2)];
        d_mid = d(on_mid);
        extra = 2 * size(left,1) + size(mid,1) - n;
        if extra < 0 || (mod(extra, 2) == 1 && isempty(mid))
            continue
        end
        % Trim from nearest the edge: a midline sensor if the excess is
        % odd, then mirrored pairs
        if mod(extra, 2) == 1
            [~, worst] = min(d_mid);
            mid(worst,:) = [];
            extra = extra - 1;
        end
        [~, order] = sort(d_left);
        left(order(1:extra/2),:) = [];
        % Front to back, then outwards
        left = sortrows([-left(:,2), -left(:,1)]);
        left = [-left(:,2), -left(:,1)];
        mid = sortrows(mid, -2);
        n_half = size(left,1);
        pos = [left; -left(:,1), left(:,2); mid];
        label = [compose('L%d', (1:n_half)'); compose('R%d', (1:n_half)'); compose('M%d', (1:size(mid,1))')];
    else
        extra = size(cand,1) - n;
        if extra < 0
            continue
        end
        [~, order] = sort(d);
        cand(order(1:extra),:) = [];
        pos = sortrows([-cand(:,2), cand(:,1)]);
        pos = [pos(:,2), -pos(:,1)];
        label = compose('%d', (1:n)');
    end
    break
end
if isempty(pos)
    error('Could not fit %d sensors inside the mask.', n);
end

end

function [pos, info] = relax_layout(pos, mask, to_head, R)
% Spread layout positions evenly over the mask by centroidal Voronoi
% relaxation (Lloyd's algorithm), measured on the head through to_head.
% Fine samples cover the mask, and each goes to its nearest sensor on the
% head, weighted by the scalp area it stands for, and each sensor moves to
% the weighted centre of its samples until they settle. Every sensor then
% covers an equal share, so the array is both even and fills the mask.
%   R.mode           'area': equal shares of the mask, the outer row about
%                    half a spacing in; 'edge': the outer row pushed out to
%                    the margins, still evenly spread
%   R.edge_margin    Sensors kept this far inside the mask edge (mm)
%   R.midline_margin Sensors off the midline kept this far from it (mm)
%   R.pairs          Kx2 mirrored left and right sensors, kept mirrored
%   R.mid            Sensors kept on the midline
%   R.iterations     Maximum iterations
%   info.iterations, info.spread (coefficient of variation of the nearest
%   neighbour distance on the head, at each iteration)
step = 0.003;
tol = 1e-4;
if ~isempty(R.mid) && R.midline_margin > 0
    warning('Sensors on the midline: the midline margin is ignored.');
    R.midline_margin = 0;
end

% Samples on a grid symmetric about the midline, so each has its mirror
m = ceil((max(vecnorm(mask, 2, 2)) + 0.1) / step);
[X, Y] = meshgrid((-m:m) * step);
mirror_full = reshape(fliplr(reshape(1:numel(X), size(X))), [], 1);

% Where each sample is on the head and how much scalp area it stands for
P = to_head([X(:), Y(:)]);
[dx1, dy1] = gradient(reshape(P(:,1), size(X)), step);
[dx2, dy2] = gradient(reshape(P(:,2), size(X)), step);
[dx3, dy3] = gradient(reshape(P(:,3), size(X)), step);
weight = vecnorm(cross([dx1(:), dx2(:), dx3(:)], [dy1(:), dy2(:), dy3(:)], 2), 2, 2);

% Distance on the head from the mask edge (positive inside) and from the
% midline, and how much room each sample has beyond the margins (slack,
% positive where a sensor may sit)
[~, edge_mm] = nearest_point(to_head(resample_polygon(mask, 2000)), P);
inside = inpolygon(X(:), Y(:), mask(:,1), mask(:,2));
edge_mm(~inside) = -edge_mm(~inside);
slack = edge_mm - R.edge_margin;
if R.midline_margin > 0
    [~, midline_mm] = nearest_point(to_head(midline_points(mask, 1000)), P);
    slack = min(slack, midline_mm - R.midline_margin);
end

keep = find(edge_mm > -50);
kept = zeros(numel(X),1);
kept(keep) = 1:numel(keep);
P = P(keep,:);
weight = weight(keep);
slack = slack(keep);
sample_2d = [X(keep), Y(keep)];

% Mirrored sensors must have room on both sides, so they use the room at
% each sample and its mirror image
mirror = kept(mirror_full(keep));
has_mirror = mirror > 0;
room = slack;
if ~isempty(R.pairs) || ~isempty(R.mid)
    room = -inf(size(slack));
    room(has_mirror) = min(slack(has_mirror), slack(mirror(has_mirror)));
end
allowed = room >= 0;
allowed_left = find(allowed & sample_2d(:,1) < 0);
allowed_right = find(allowed & sample_2d(:,1) > 0);
allowed_mid = find(allowed & sample_2d(:,1) == 0);

n = size(pos,1);
info.spread = nan(R.iterations,1);
for it = 1:R.iterations
    sensor_head = to_head(pos);
    nearest_mm = nearest_other(sensor_head);
    info.spread(it) = std(nearest_mm) / mean(nearest_mm);

    % The region shared out: up to the mask edge for 'area', or for 'edge'
    % half a spacing beyond the margins, so the outer sensors settle on them
    grow = R.edge_margin;
    if strcmp(R.mode, 'edge')
        grow = max(median(nearest_mm)/2, R.edge_margin);
    end
    in_region = slack > -grow;

    owner = nearest_point(sensor_head, P(in_region,:));
    w = weight(in_region);
    total = accumarray(owner, w, [n 1]);
    centre = [accumarray(owner, w .* sample_2d(in_region,1), [n 1]), ...
        accumarray(owner, w .* sample_2d(in_region,2), [n 1])] ./ total;
    empty = total == 0;
    centre(empty,:) = pos(empty,:);

    % Keep the symmetry: mirrored pairs averaged on the left, midline
    % sensors on it
    if ~isempty(R.pairs)
        l = R.pairs(:,1);
        r = R.pairs(:,2);
        left = (centre(l,:) + [-centre(r,1), centre(r,2)]) / 2;
        left(:,1) = min(left(:,1), -1e-3);
        centre(l,:) = left;
    end
    centre(R.mid,1) = 0;

    % Sensors without room are moved to the nearest place with room on
    % their own side, then the right of each pair mirrors its left
    c = min(max(round(centre / step) + m + 1, 1), 2*m + 1);
    near = kept(sub2ind(size(X), c(:,2), c(:,1)));
    short = true(n,1);
    short(near > 0) = room(near(near > 0)) < 0;
    short(R.pairs(:,2)) = false;
    for k = find(short)'
        if ismember(k, R.mid)
            candidates = allowed_mid;
        elseif ismember(k, R.pairs(:,1)) || centre(k,1) < 0
            candidates = allowed_left;
        else
            candidates = allowed_right;
        end
        if isempty(candidates)
            continue;
        end
        centre(k,:) = sample_2d(candidates(nearest_point(sample_2d(candidates,:), centre(k,:))),:);
    end
    centre(R.mid,1) = 0;
    if ~isempty(R.pairs)
        centre(R.pairs(:,2),:) = [-centre(R.pairs(:,1),1), centre(R.pairs(:,1),2)];
    end

    moved = max(vecnorm(centre - pos, 2, 2));
    pos = centre;
    if moved < tol, break; end
end
info.iterations = it;
info.spread = info.spread(1:it);

end

function [pos, label] = picked_labels(pos2d, symmetry, allow_midline)
% Positions picked on a head, projected into 2D. symmetry 'left' or
% 'right' keeps the sensors picked on that side and mirrors them, 'none'
% keeps them all. Sensors on the midline (x = 0) are kept once as M1.. if
% allow_midline, otherwise dropped.
tol = 1e-3;
if strcmp(symmetry, 'none')
    pos = pos2d;
    label = number_labels(size(pos,1));
    return
end
side = -1;
if strcmp(symmetry, 'right'), side = 1; end
on_side = side * pos2d(:,1) > tol;
on_mid = abs(pos2d(:,1)) <= tol;
if any(on_mid) && ~allow_midline
    warning('%d picked sensors are on the midline and were dropped.', sum(on_mid));
end
if any(side * pos2d(:,1) < -tol)
    fprintf('Mirroring the %s side: ignoring %d sensors picked on the other side\n', symmetry, sum(side * pos2d(:,1) < -tol));
end
half = pos2d(on_side,:);
half = sortrows([-half(:,2), abs(half(:,1))]);
left = [-half(:,2), -half(:,1)];
n_half = size(left,1);
pos = [left; -left(:,1), left(:,2)];
label = [compose('L%d', (1:n_half)'); compose('R%d', (1:n_half)')];
if allow_midline && any(on_mid)
    mid = sortrows([zeros(sum(on_mid),1), pos2d(on_mid,2)], -2);
    pos = [pos; mid];
    label = [label; compose('M%d', (1:size(mid,1))')];
end
end

function neighbours = layout_neighbours(pos)
% Neighbours from a Delaunay triangulation of the positions, leaving out
% edges much longer than the typical nearest neighbour spacing
neighbours = cell(size(pos,1),1);
if size(pos,1) < 3, return; end
tri = delaunay(pos(:,1), pos(:,2));
edges = unique(sort([tri(:,[1 2]); tri(:,[2 3]); tri(:,[3 1])], 2), 'rows');
len = vecnorm(pos(edges(:,1),:) - pos(edges(:,2),:), 2, 2);
edges = edges(len < 1.6 * median(nearest_other(pos)),:);
for e = 1:size(edges,1)
    neighbours{edges(e,1)}(end+1) = edges(e,2);
    neighbours{edges(e,2)}(end+1) = edges(e,1);
end
end

function m = measure(pos, mask, to_head)
% Spacing and margins of a layout measured through to_head (on the head, or
% 2D scaled to mm). Per sensor: nearest_mm to the nearest other, edge_mm
% inside the mask edge (negative outside), midline_mm from the midline.
% Summaries of the spacing: median, min and max, and spacing_cv, its spread
% relative to the median (0 for a perfectly even array).
h = to_head(pos);
m.nearest_mm = nearest_other(h);
[~, m.edge_mm] = nearest_point(to_head(resample_polygon(mask, 2000)), h);
outside = ~inpolygon(pos(:,1), pos(:,2), mask(:,1), mask(:,2));
m.edge_mm(outside) = -m.edge_mm(outside);
[~, m.midline_mm] = nearest_point(to_head(midline_points(mask, 1000)), h);
m.midline_mm(abs(pos(:,1)) < 1e-9) = 0;
m.median_spacing = median(m.nearest_mm);
m.min_spacing = min(m.nearest_mm);
m.max_spacing = max(m.nearest_mm);
m.spacing_cv = std(m.nearest_mm) / m.median_spacing;
end

function outline = symmetric_outline(outline, side)
% Make an outline symmetric by mirroring one side across the midline
if strcmp(side, 'none'), return; end
big = 10;
shape = polyshape(outline);
switch side
    case 'left'
        half = intersect(shape, polyshape([-big -big; 0 -big; 0 big; -big big]));
    case 'right'
        half = intersect(shape, polyshape([0 -big; big -big; big big; 0 big]));
    otherwise
        error('Unknown mask symmetry %s.', side);
end
half = largest_region(half);
% Snap the cut edge exactly onto the midline so the halves share it
vertices = half.Vertices;
vertices(abs(vertices(:,1)) < 1e-9, 1) = 0;
half = polyshape(vertices(:,1), vertices(:,2));
mirrored = polyshape(-vertices(:,1), vertices(:,2));
joined = rmholes(largest_region(union(half, mirrored)));
if area(joined) < 1.99 * area(half)
    warning('Mirrored mask lost %.1f%% of its area joining the halves.', 100 * (1 - area(joined) / (2 * area(half))));
end
outline = joined.Vertices;

end


function shape = largest_region(shape)
% The biggest piece of a polyshape
pieces = regions(shape);
[~, biggest] = max(arrayfun(@area, pieces));
shape = pieces(biggest);

end


%==========================================================================
% Geometry
%==========================================================================
function [idx, d] = nearest_point(X, Q)
% Nearest row of X to each row of Q, and the distance, in blocks
idx = zeros(size(Q,1),1);
d = zeros(size(Q,1),1);
block = max(1, floor(5e6 / size(X,1)));
for i = 1:block:size(Q,1)
    j = i:min(i + block - 1, size(Q,1));
    [d(j), idx(j)] = min(distances(Q(j,:), X), [], 2);
end

end

function d = nearest_other(X)
% Distance from each row of X to the nearest other row
D = distances(X, X);
D(1:size(X,1)+1:end) = inf;
d = min(D, [], 2);

end

function D = distances(A, B)
% Distances between the rows of A and the rows of B
D = sqrt(max(sum(A.^2, 2) + sum(B.^2, 2)' - 2 * (A * B'), 0));
end

function pts = resample_polygon(poly, n)
% n points evenly spaced along the perimeter of a closed polygon
poly = [poly; poly(1,:)];
s = [0; cumsum(vecnorm(diff(poly), 2, 2))];
[s, keep] = unique(s);
pts = interp1(s, poly(keep,:), linspace(0, s(end), n+1)');
pts = pts(1:end-1,:);
end

function poly = open_polygon(poly)
% A closed outline without its first point repeated at the end
if size(poly,1) > 1 && isequal(poly(1,:), poly(end,:))
    poly = poly(1:end-1,:);
end
end

function split = midline_split(mask)
% Line back to front through Cz, past the mask at both ends: the midline,
% and the split of a scannercast helmet
y = linspace(min(mask(:,2)) - 0.05, max(mask(:,2)) + 0.05, 21)';
split = [zeros(size(y)), y];
end

function p = midline_points(mask, n)
% n points along the midline, past the mask and the grids round it
r = max(vecnorm(mask, 2, 2)) + 0.1;
p = [zeros(n,1), linspace(-r, r, n)'];
end

function label = number_labels(n)
label = compose('%d', (1:n)');
end

%==========================================================================
% Picking
%==========================================================================
function pts = pick_points(mesh, fid, prompt, closed, opts)
% Click points on a mesh in order. Backspace removes the last point, enter
% finishes. Turn the rotate tool off before clicking. closed (default true)
% only changes how the picked line is drawn.
%   opts.outline  Nx3 closed line drawn in red (the mask)
%   opts.split    Mx3 line drawn in green (the midline)
%   opts.radius   circle of this radius drawn round each point, square to
%                 the surface (a sensor footprint); the points are not joined
%   opts.lift     circles drawn this far off the surface
%   opts.check    function of the picked points returning a message per
%                 point, '' if it is fine. Points with a message are drawn
%                 red with the message beside them.
if nargin < 4, closed = true; end
if nargin < 5, opts = struct(); end
footprint = isfield(opts, 'radius');

fig = figure;
hold on
axis equal
h = trisurf(mesh, 'EdgeAlpha', 0, 'FaceColor', [0.9 0.75 0.65]);
fid_names = fieldnames(fid);
for k = 1:numel(fid_names)
    p = fid.(fid_names{k});
    scatter3(p(1), p(2), p(3), 60, 'k', 'filled');
    text(p(1), p(2), p(3), fid_names{k}, 'fontweight', 'bold', 'verticalalignment', 'bottom');
end
if isfield(opts, 'outline')
    plot3(opts.outline([1:end 1],1), opts.outline([1:end 1],2), opts.outline([1:end 1],3), 'r-', 'linewidth', 2, 'hittest', 'off');
end
if isfield(opts, 'split')
    plot3(opts.split(:,1), opts.split(:,2), opts.split(:,3), 'g-', 'linewidth', 2, 'hittest', 'off');
end
line_style = 'r.-';
if footprint, line_style = 'k.'; end
line_h = plot3(nan, nan, nan, line_style, 'markersize', 20, 'linewidth', 1.5, 'hittest', 'off');
ok_h = plot3(nan, nan, nan, 'g-', 'linewidth', 1.5, 'hittest', 'off');
bad_h = plot3(nan, nan, nan, 'r-', 'linewidth', 1.5, 'hittest', 'off');
label_h = gobjects(0);
view(3);
follow_light(gca);
title({prompt, 'Click to add, backspace to undo, enter to finish'});

centre = mean(mesh.Points);
pts = zeros(0,3);
h.ButtonDownFcn = @add_point;
fig.KeyPressFcn = @key_press;
uiwait(fig);
if isvalid(fig), close(fig); end

    function add_point(~, event)
        pts(end+1,:) = event.IntersectionPoint;
        redraw();
    end

    function key_press(~, event)
        switch event.Key
            case 'backspace'
                pts = pts(1:end-1,:);
                redraw();
            case 'return'
                uiresume(fig);
        end
    end

    function redraw()
        shown = pts;
        if closed && ~footprint
            shown = pts([1:end min(1,end)],:);
        end
        set(line_h, 'XData', shown(:,1), 'YData', shown(:,2), 'ZData', shown(:,3));

        messages = repmat({''}, size(pts,1), 1);
        if isfield(opts, 'check') && ~isempty(pts)
            messages = opts.check(pts);
        end
        bad = ~cellfun(@isempty, messages);
        if footprint
            set(ok_h, {'XData', 'YData', 'ZData'}, circles(pts(~bad,:)));
            set(bad_h, {'XData', 'YData', 'ZData'}, circles(pts(bad,:)));
        end

        delete(label_h);
        label_h = gobjects(size(pts,1),1);
        for j = 1:size(pts,1)
            label = num2str(j);
            colour = 'k';
            if bad(j)
                label = [label, ': ', messages{j}]; %#ok<AGROW>
                colour = 'r';
            end
            label_h(j) = text(pts(j,1), pts(j,2), pts(j,3), ['  ', label], 'color', colour, ...
                'fontweight', 'bold', 'fontsize', 8, 'hittest', 'off');
        end

        status = 'Click to add, backspace to undo, enter to finish';
        if any(bad)
            status = sprintf('%d point(s) flagged in red. %s', sum(bad), status);
        end
        title({prompt, status});
    end

    function xyz = circles(c)
        % Circle round each point square to the surface normal, lifted off
        % the surface, separated by NaNs
        theta = [linspace(0, 2*pi, 41), nan]';
        xyz = zeros(0,3);
        for j = 1:size(c,1)
            nrm = vertexNormal(mesh, nearestNeighbor(mesh, c(j,:)));
            if dot(nrm, c(j,:) - centre) < 0, nrm = -nrm; end
            a = null(nrm)';
            ring = c(j,:) + opts.lift * nrm + opts.radius * (cos(theta) * a(1,:) + sin(theta) * a(2,:));
            xyz = [xyz; ring]; %#ok<AGROW>
        end
        xyz = {xyz(:,1), xyz(:,2), xyz(:,3)};
    end
end

function messages = check_picks(pts, chk)
% Rough checks on sensors as they are picked, one message per point ('' if
% it looks fine). Each pick's layout position is that of the nearest point
% of a fine grid mapped onto the head; distances are on the head.
q = chk.fine_2d(nearest_point(chk.fine_3d, pts),:);
[~, edge_mm] = nearest_point(chk.edge_3d, pts);
inside = inpolygon(q(:,1), q(:,2), chk.mask(:,1), chk.mask(:,2));
[~, midline_mm] = nearest_point(chk.midline_3d, pts);
on_mid = chk.allow_midline & midline_mm < chk.midline_snap;

spacing = distances(pts, pts);
spacing(1:size(pts,1)+1:end) = inf;
if ~strcmp(chk.side, 'none')
    % Mirrors of the picks, including each one's own, which a midline
    % sensor does not have
    mirror = distances(pts, chk.to_head([-q(:,1), q(:,2)]));
    mirror(sub2ind(size(mirror), find(on_mid), find(on_mid))) = inf;
    spacing = min(spacing, mirror);
end

messages = repmat({''}, size(pts,1), 1);
for k = 1:size(pts,1)
    issues = {};
    if ~inside(k)
        issues{end+1} = 'outside the edge'; %#ok<AGROW>
    elseif edge_mm(k) < chk.edge_margin
        issues{end+1} = sprintf('%.0f mm from the edge', edge_mm(k)); %#ok<AGROW>
    end
    if on_mid(k)
        % On the midline
    elseif (strcmp(chk.side, 'left') && q(k,1) > 0) || (strcmp(chk.side, 'right') && q(k,1) < 0)
        issues{end+1} = 'wrong side'; %#ok<AGROW>
    elseif midline_mm(k) < chk.midline_margin
        issues{end+1} = sprintf('%.0f mm from the midline', midline_mm(k)); %#ok<AGROW>
    end
    if min(spacing(k,:)) < chk.min_spacing
        issues{end+1} = sprintf('%.0f mm from a neighbour', min(spacing(k,:))); %#ok<AGROW>
    end
    messages{k} = strjoin(issues, ', ');
end
end

%==========================================================================
% Plotting
%==========================================================================
function follow_light(ax)
% Light that stays fixed relative to the camera, so the surface always
% shows some shading across the view as it is rotated
delete(findobj(ax, 'type', 'light'));
h = camlight(ax, 40, 30);
material(ax, 'dull');
lighting(ax, 'gouraud');
addlistener(ax, 'View', 'PostSet', @(~,~) move_light(h));
end

function move_light(h)
if isvalid(h)
    camlight(h, 40, 30);
end
end

function plot_layout_reference()
% Head outline of the anatomical layout as spm_get_anatomical_layout draws
% it.
samples = 128;
theta = linspace(0, 2*pi, samples)';
circ_r = 0.25;
rad2 = 0.3125;
head = rad2 * [cos(theta), sin(theta)];
[~, index] = min(abs(theta - 0.035));
nose = [head((samples/4) - (index - 1), :); 0, rad2 * 1.1];
nose(3,:) = nose(1,:) .* [-1 1];
[~, index] = min(abs(theta - 0.08));
ear = [head((samples/2) - (index - 1), :); ...
        -rad2 * 1.037, rad2 * 0.235; ...
        -rad2 * 1.07,  rad2 * 0.2; ...
        -rad2 * 1.08,  rad2 * 0.165; ...
        -rad2 * 1.05,  rad2 * 0.02; ...
        -rad2 * 1.09, -rad2 * 0.15; ...
        -rad2 * 1.065, -rad2 * 0.21];
ear(8,:) = ear(1,:) .* [1 -1];
plot(circ_r * cos(theta), circ_r * sin(theta), 'k--');
plot(head(:,1), head(:,2), 'k');
plot(nose(:,1), nose(:,2), 'k');
plot(ear(:,1), ear(:,2), 'k');
plot(-ear(:,1), ear(:,2), 'k');
text([0 0 -circ_r circ_r 0], [0 circ_r 0 0 -circ_r], {'Cz', 'FPz', 'T3', 'T4', 'Oz'}, ...
    'fontsize', 8, 'color', [0.4 0.4 0.4], 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
end

function plot_neighbours(pos, neighbours)
% Dashed lines between neighbours
for k = 1:numel(neighbours)
    for j = neighbours{k}(:)'
        if j > k
            plot(pos([k j],1), pos([k j],2), 'k--', 'linewidth', 0.5);
        end
    end
end
end
