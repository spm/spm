function [shape] = ft_read_headshape(filename, varargin)

% FT_READ_HEADSHAPE reads the fiducials and/or the measured headshape and/or meshes
% that describe headsurfaces, headmodels or cortical sheets from a variety of files
% (like CTF and Polhemus). The headshape and fiducials can for example be used for
% coregistration.
%
% Use as
%   [shape] = ft_read_headshape(filename, ...)
% or
%   [shape] = ft_read_headshape({filename1, filename2}, ...)
%
% If you specify the filename as a cell-array, the following situations are supported:
%  - a two-element cell-array with the file names for the left and right hemisphere,
%    e.g., FreeSurfer's {'lh.orig' 'rh.orig'}, or Caret's {'X.L.Y.Z.surf.gii' 'X.R.Y.Z.surf.gii'}
%  - a two-element cell-array points to files that represent the coordinates and topology
%    in separate files, e.g., Caret's {'X.L.Y.Z.coord.gii' 'A.L.B.C.topo.gii'}
% By default all information from the two files will be assumed to correspond to the
% left and right hemispeheres and concatenated. The option 'concatenate' can be set
% to 'no' to prevent them from being concatenated in a single structure.
%
% Additional options should be specified in key-value pairs and can include
%   'format'      = string, see below
%   'coordsys'    = string, e.g. 'head' or 'dewar' (only supported for CTF)
%   'unit'        = string, e.g. 'mm' (default is the native units of the file)
%   'concatenate' = 'yes' or 'no' (default = 'yes')
%   'image'       = path to corresponding image/jpg file
%   'surface'     = specific surface to be read (only for Caret spec files)
%   'refine'      = number, used for refining Structure Sensor meshes (default = 1)
%   'jmeshopt'    = cell-array with {'name', 'value'} pairs, options for reading JSON/JMesh files
%   'meshtype'    = string, which type of mesh to read in case the file contains multiple types, can be 'tri', 'tet' or 'hex'
%
% Supported input file formats include
%   'gifti'           see https://www.nitrc.org/projects/gifti/
%   'gmsh_ascii'      see https://gmsh.info
%   'gmsh_binary'     see https://gmsh.info
%   'matlab'          containing FieldTrip or BrainStorm headshapes or cortical meshes
%   'mne_tri'         MNE surface description in ASCII format
%   'mne_pos'         MNE source grid in ascii format, described as 3D points
%   'neurojson_jmesh' NeuroJSON ascii JSON-based format
%   'neurojson_bmesh' NeuroJSON binary JSON-based format
%   'obj'             Wavefront .obj file obtained with the Structure Sensor
%   'off'             see http://www.geomview.org/docs/html/OFF.html
%   'ply'             Stanford Polygon file format, for use with Paraview or Meshlab
%   'stl'             STereoLithography file format, for use with CAD and/or generic 3D mesh editing programs
%   'tck'             Mrtrix track file
%   'trk'             Trackvis trk file
%   'vista'           see http://www.cs.ubc.ca/nest/lci/vista/vista.html
%   'vtk'             Visualization ToolKit file format, for use with Paraview
%   'vtk_xml'         Visualization ToolKit file format
%   'itab_asc'
%   'ctf_*'
%   '4d_*'
%   'neuromag_*'
%   'yokogawa_*'
%   'yorkinstruments_hdf5'
%   'polhemus_*'
%   'freesurfer_*'
%   'mne_source'
%   'spmeeg_mat'
%   'netmeg'
%   'tet'
%   'tetgen_ele'
%   'caret_surf'
%   'caret_coord'
%   'caret_topo'
%   'caret_spec'
%   'brainvisa_mesh'
%   'brainsuite_dfs'
%
% See also FT_READ_HEADMODEL, FT_READ_SENS, FT_READ_ATLAS, FT_WRITE_HEADSHAPE

% Copyright (C) 2008-2025, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% get the options
annotationfile = ft_getopt(varargin, 'annotationfile');
useimage       = ft_getopt(varargin, 'useimage', true);      % use image if hasimage
concatenate    = ft_getopt(varargin, 'concatenate', 'yes');
coordsys       = ft_getopt(varargin, 'coordsys', 'head');    % for ctf or neuromag_mne coil positions, the alternative is dewar
fileformat     = ft_getopt(varargin, 'format');
unit           = ft_getopt(varargin, 'unit');
image          = ft_getopt(varargin, 'image');               % path to .jpg file
surface        = ft_getopt(varargin, 'surface');
refine_        = ft_getopt(varargin, 'refine', 1);           % do not confuse with the private/refine function
meshtype       = ft_getopt(varargin, 'meshtype');            % tri, tet or hex

% Check the input, if filename is a cell-array, call FT_READ_HEADSHAPE recursively and combine the outputs.
% This is used to read the left and right hemisphere of a Freesurfer cortical segmentation.
if iscell(filename)

  for i = 1:numel(filename)
    tmp       = ft_read_headshape(filename{i}, varargin{:});
    haspos(i) = isfield(tmp, 'pos') && ~isempty(tmp.pos);
    hastri(i) = isfield(tmp, 'tri') && ~isempty(tmp.tri);
    if ~haspos(i), tmp.pos = []; end
    if ~hastri(i), tmp.tri = []; end
    if ~isfield(tmp, 'unit'), tmp.unit = 'unknown'; end
    bnd(i) = tmp;
  end

  % Concatenate the meshes (only if 'concatenate' = 'yes' ) and if all
  % structures have non-empty vertices and triangles. If not, the input filenames
  % may have been caret-style coord and topo, which needs combination of
  % the pos and tri.

  if  numel(filename)>1 && all(haspos==1) && strcmp(concatenate, 'yes')
    if length(bnd)>2
      ft_error('Cannot concatenate more than two files') % no more than two files are taken for cancatenation
    else
      fprintf('Concatenating the meshes in %s and %s\n', filename{1}, filename{2});

      shape     = [];
      shape.pos = cat(1, bnd.pos);
      npos      = size(bnd(1).pos,1);

      if isfield(bnd(1), 'tri')  && isfield(bnd(2), 'tri')
        shape.tri = cat(1, bnd(1).tri, bnd(2).tri + npos);
      elseif ~isfield(bnd(1), 'tri') && ~isfield(bnd(2), 'tri')
        % this is ok
      else
        ft_error('not all input files seem to contain a triangulation');
      end

      % concatenate any other fields
      fnames = {'sulc' 'curv' 'area' 'thickness' 'atlasroi'};
      for k = 1:numel(fnames)
        if isfield(bnd(1), fnames{k}) && isfield(bnd(2), fnames{k})
          shape.(fnames{k}) = cat(1, bnd.(fnames{k}));
        elseif ~isfield(bnd(1), fnames{k}) && ~isfield(bnd(2), fnames{k})
          % this is ok
        else
          ft_error('not all input files seem to contain a "%s"', fnames{k});
        end
      end

      shape.brainstructure = []; % keeps track of the order of files in concatenation
      for h = 1:length(bnd)
        shape.brainstructure  = [shape.brainstructure; h*ones(length(bnd(h).pos), 1)];
        [p,f,e]               = fileparts(filename{h});

        % do an educated guess, otherwise default to the filename
        iscortexright = contains(f,'rh');
        iscortexright = iscortexright || contains(f,'.R.');
        iscortexright = iscortexright || contains(f,'Right');
        iscortexright = iscortexright || contains(f,'RIGHT');

        iscortexleft = contains(f,'lh');
        iscortexleft = iscortexleft || contains(f,'.L.');
        iscortexleft = iscortexleft || contains(f,'Left');
        iscortexleft = iscortexleft || contains(f,'LEFT');

        if iscortexright && iscortexleft
          % something strange is going on, default to the filename and let the user take care of this
          shape.brainstructurelabel{h,1} = f;
        elseif iscortexleft
          shape.brainstructurelabel{h,1} = 'CORTEX_LEFT';
        elseif iscortexright
          shape.brainstructurelabel{h,1} = 'CORTEX_RIGHT';
        else
          % nothing to be guessed
          shape.brainstructurelabel{h,1} = f;
        end
      end

    end
  elseif numel(filename)>1 && ~all(haspos==1)
    if numel(bnd)>2
      ft_error('Cannot combine more than two files') % no more than two files are taken for cancatenation
    else
      shape = [];
      if sum(haspos==1)==1
        fprintf('Using the vertex positions from %s\n', filename{haspos==1});
        shape.pos  = bnd(haspos==1).pos;
        shape.unit = bnd(haspos==1).unit;
      else
        ft_error('Don''t know what to do');
      end
      if sum(hastri==1)==1
        fprintf('Using the faces definition from %s\n', filename{hastri==1});
        shape.tri = bnd(hastri==1).tri;
      end
      if max(shape.tri(:))~=size(shape.pos,1)
        ft_error('mismatch in number of points in pos and tri');
      end
    end

  else
    % in case numel(filename)==1, or strcmp(concatenate, 'no')
    shape = bnd;
  end

  return
end % if iscell

% some of the code does not work well on matlab strings, (i.e. "" vs ''),
% specifically ["a" "b"] yields something different than ['a' 'b'].
if isstring(filename)
  filename = char(filename);
end

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

if isempty(fileformat)
  % only do the autodetection if the format was not specified
  fileformat = ft_filetype(filename);
end

% checks if there exists a .jpg file of 'filename'
[pathstr,name]  = fileparts(filename);
if useimage
  if exist(fullfile(pathstr,[name,'.jpg']), 'file')
    image    = fullfile(pathstr,[name,'.jpg']);
    hasimage = true;
  elseif exist(fullfile(pathstr,[name,'.JPG']), 'file')
    image    = fullfile(pathstr,[name,'.JPG']);
    hasimage = true;
  elseif exist(fullfile(pathstr,[name,'.png']), 'file')
    image    = fullfile(pathstr,[name,'.png']);
    hasimage = true;
  elseif exist(fullfile(pathstr,[name,'.PNG']), 'file')
    image    = fullfile(pathstr,[name,'.PNG']);
    hasimage = true;
  else
    hasimage = false;
  end
else
  hasimage = false;
end

if ~isempty(annotationfile) && ~strcmp(fileformat, 'mne_source')
  ft_error('extracting annotation information only works for ''mne_source'' files');
end

% start with an empty structure
shape     = [];
shape.pos = [];

switch fileformat
  case {'ctf_ds', 'ctf_hc', 'ctf_meg4', 'ctf_res4', 'ctf_old'}
    [p, f, x] = fileparts(filename);

    if strcmp(fileformat, 'ctf_old')
      fileformat = ft_filetype(filename);
    end

    if strcmp(fileformat, 'ctf_ds')
      filename = fullfile(p, [f x], [f '.hc']);
    elseif strcmp(fileformat, 'ctf_meg4')
      filename = fullfile(p, [f '.hc']);
    elseif strcmp(fileformat, 'ctf_res4')
      filename = fullfile(p, [f '.hc']);
    end

    orig = read_ctf_hc(filename);
    switch coordsys
      case 'head'
        shape.fid.pos = cell2mat(struct2cell(orig.head));
        shape.coordsys = 'ctf';
      case 'dewar'
        shape.fid.pos = cell2mat(struct2cell(orig.dewar));
        shape.coordsys = 'dewar';
      otherwise
        ft_error('incorrect coordsys specified');
    end
    shape.fid.label = fieldnames(orig.head);

  case 'ctf_shape'
    orig = read_ctf_shape(filename);
    shape.pos = orig.pos;

    % The file also contains fiducial information, but those are in MRI voxels and
    % inconsistent with the headshape itself.
    %
    % shape.fid.label = {'NASION', 'LEFT_EAR', 'RIGHT_EAR'};
    % shape.fid.pos = zeros(0,3); % start with an empty array
    % for i = 1:numel(shape.fid.label)
    %   shape.fid.pos = cat(1, shape.fid.pos, getfield(orig.MRI_Info, shape.fid.label{i}));
    % end

  case {'4d_xyz', '4d_m4d', '4d_hs', '4d', '4d_pdf'}
    [p, f, x] = fileparts(filename);
    if ~strcmp(fileformat, '4d_hs')
      filename = fullfile(p, 'hs_file');
    end
    [shape.pos, fid] = read_bti_hs(filename);

    % I'm making some assumptions here
    % which I'm not sure will work on all 4D systems

    % fid = fid(1:3, :);

    [junk, NZ] = max(fid(1:3,1));
    [junk, L]  = max(fid(1:3,2));
    [junk, R]  = min(fid(1:3,2));
    rest       = setdiff(1:size(fid,1),[NZ L R]);

    shape.fid.pos = fid([NZ L R rest], :);
    shape.fid.label = {'NZ', 'L', 'R'};
    if ~isempty(rest)
      for i = 4:size(fid,1)
        shape.fid.label{i} = ['fiducial' num2str(i)];
        % in a 5 coil configuration this corresponds with Cz and Inion, or
        % wherever the researcher has placed the coils
      end
    end
    shape.coordsys = '4d';

  case 'itab_asc'
    shape = read_itab_asc(filename);

  case 'gifti'
    ft_hastoolbox('gifti', 1);
    g = gifti(filename);
    if ~isfield(g, 'vertices')
      ft_error('%s does not contain a tesselated surface', filename);
    end
    shape.pos = ft_warp_apply(g.mat, g.vertices);
    shape.tri = g.faces;
    shape.unit = 'mm';  % defined in the GIFTI standard to be milimeter
    if isfield(g, 'cdata')
      shape.mom = g.cdata;
    end

  case {'caret_surf' 'caret_topo' 'caret_coord'}
    ft_hastoolbox('gifti', 1);
    g = gifti(filename);
    if ~isfield(g, 'vertices') && strcmp(fileformat, 'caret_topo')
      try
        % do a clever guess by replacing topo with coord
        g2 = gifti(strrep(filename, '.topo.', '.coord.'));
        vertices  = ft_warp_apply(g2.mat, g2.vertices);
      catch
        vertices  = [];
      end
    else
      vertices  = ft_warp_apply(g.mat, g.vertices);
    end
    if ~isfield(g, 'faces') && strcmp(fileformat, 'caret_coord')
      try
        % do a clever guess by replacing topo with coord
        g2 = gifti(strrep(filename, '.coord.', '.topo.'));
        faces = g2.faces;
      catch
        faces = [];
      end
    else
      faces = g.faces;
    end

    shape.pos = vertices;
    shape.tri = faces;
    if isfield(g, 'cdata')
      shape.mom = g.cdata;
    end

    % check whether there is curvature info etc
    filename    = strrep(filename, '.surf.', '.shape.');
    filename    = strrep(filename, '.topo.', '.shape.');
    filename    = strrep(filename, '.coord.', '.shape.');

    [p,f,e]     = fileparts(filename);
    tok         = tokenize(f, '.');
    if length(tok)>2
      tmpfilename = strrep(filename, tok{3}, 'sulc');
      if exist(tmpfilename, 'file'), g = gifti(tmpfilename); shape.sulc = g.cdata; end
      if exist(strrep(tmpfilename, 'sulc', 'curvature'), 'file'),  g = gifti(strrep(tmpfilename, 'sulc', 'curvature')); shape.curv = g.cdata; end
      if exist(strrep(tmpfilename, 'sulc', 'thickness'), 'file'),  g = gifti(strrep(tmpfilename, 'sulc', 'thickness')); shape.thickness = g.cdata; end
      if exist(strrep(tmpfilename, 'sulc', 'atlasroi'),  'file'),  g = gifti(strrep(tmpfilename, 'sulc', 'atlasroi'));  shape.atlasroi  = g.cdata; end
    end

  case 'caret_spec'
    [p, f, e] = fileparts(filename);
    [spec, headerinfo] = read_caret_spec(filename);
    fn = fieldnames(spec);

    % concatenate the filenames that contain coordinates
    % concatenate the filenames that contain topologies
    coordfiles = {};
    topofiles  = {};
    for k = 1:numel(fn)
      if ~isempty(strfind(fn{k}, 'topo'))
        topofiles = cat(1,topofiles, fullfile(p,spec.(fn{k})));
      end
      if ~isempty(strfind(fn{k}, 'coord'))
        coordfiles = cat(1,coordfiles, fullfile(p,spec.(fn{k})));
      end
    end
    if isempty(surface)
      [selcoord, ok] = listdlg('ListString',coordfiles,'SelectionMode','single','PromptString','Select a file describing the coordinates');
    else
      selcoord = find(contains(coordfiles, surface));
    end
    if numel(topofiles)>1
      [seltopo, ok]  = listdlg('ListString',topofiles,'SelectionMode','single','PromptString','Select a file describing the topology');
    else
      seltopo = 1;
    end

    % recursively call ft_read_headshape
    tmp1 = ft_read_headshape(coordfiles{selcoord});
    tmp2 = ft_read_headshape(topofiles{seltopo});

    % quick and dirty sanity check to see whether the indexing of the
    % points in the topology matches the number of points
    if max(tmp2.tri(:))~=size(tmp1.pos,1)
      ft_error('there''s a mismatch between the number of points used in the topology, and described by the coordinates');
    end

    shape     = tmp1;
    shape.tri = tmp2.tri;

  case 'neuromag_mex'
    [co,ki,nu] = hpipoints(filename);
    fid = co(:,ki==1)';

    [junk, NZ] = max(fid(:,2));
    [junk, L]  = min(fid(:,1));
    [junk, R]  = max(fid(:,1));

    shape.fid.pos = fid([NZ L R], :);
    shape.fid.label = {'NZ', 'L', 'R'};

  case 'mne_source'
    % read the source space from an MNE file
    ft_hastoolbox('mne', 1);

    try
      % a fif-file can also contain a source space that is volumetric, in which case the below function call will fail (due to the add_geom
      % being specified as true, but the file does not contain triangulation information. strictly speaking, the fif-file then
      % does not represent a headshape, but as a service to the casual user, let's support the reading of this type of file
      src = mne_read_source_spaces(filename, 1);
    catch
      src = mne_read_source_spaces(filename);
    end

    if numel(src)==2
      if ~isempty(annotationfile)
        ft_hastoolbox('freesurfer', 1);
        if numel(annotationfile)~=2
          ft_error('two annotationfiles expected, one for each hemisphere');
        end
        for k = 1:numel(annotationfile)
          [v{k}, label{k}, c(k)] = read_annotation(annotationfile{k}, 1);
        end

        % match the annotations with the src structures
        if src(1).np == numel(label{1}) && src(2).np == numel(label{2})
          src(1).labelindx = label{1};
          src(2).labelindx = label{2};
        elseif src(1).np == numel(label{2}) && src(1).np == numel(label{1})
          src(1).labelindx = label{2};
          src(2).labelindx = label{1};
        else
          ft_warning('incompatible annotation with triangulations, not using annotation information');
        end
        if ~isequal(c(1),c(2))
          ft_error('the annotation tables differ, expecting equal tables for the hemispheres');
        end
        c = c(1);
      end

      shape = [];
      % only keep the points that are in use
      inuse1 = src(1).inuse==1;
      inuse2 = src(2).inuse==1;
      shape.pos=[src(1).rr(inuse1,:); src(2).rr(inuse2,:)];

      % only keep the triangles that are in use; these have to be renumbered
      newtri1 = src(1).use_tris;
      newtri2 = src(2).use_tris;
      for i=1:numel(src(1).vertno)
        newtri1(newtri1==src(1).vertno(i)) = i;
      end
      for i=1:numel(src(2).vertno)
        newtri2(newtri2==src(2).vertno(i)) = i;
      end
      shape.tri  = [newtri1; newtri2 + numel(src(1).vertno)];
      if isfield(src(1), 'use_tri_area')
        shape.area = [src(1).use_tri_area(:); src(2).use_tri_area(:)];
      end
      if isfield(src(1), 'use_tri_nn')
        shape.nn = [src(1).use_tri_nn; src(2).use_tri_nn];
      end
      shape.orig.pos = [src(1).rr; src(2).rr];
      shape.orig.tri = [src(1).tris; src(2).tris + src(1).np];
      shape.orig.inuse = [src(1).inuse src(2).inuse]';
      shape.orig.nn    = [src(1).nn; src(2).nn];
      if isfield(src(1), 'labelindx')
        shape.orig.labelindx = [src(1).labelindx;src(2).labelindx];
        shape.labelindx      = [src(1).labelindx(inuse1); src(2).labelindx(inuse2)];
        shape.label          = c.struct_names;
        shape.annotation     = c.orig_tab; % to be able to recover which one
        shape.ctable         = c.table;
      end
    else
      ft_warning('the fif-file did not seem to contain triangulation information, probably you try to read a volumetric source space');

      shape        = [];
      shape.pos    = src.rr;
      shape.inside = src.inuse(:)>0;
      shape.dim    = pos2dim3d(src.rr);
    end

  case {'neuromag_fif' 'neuromag_mne'}

    orig = read_neuromag_hc(filename);
    switch coordsys
      case 'head'
        fidN=1;
        posN=1;
        for i=1:size(orig.head.pos,1)
          if strcmp(orig.head.label{i}, 'LPA') || strcmp(orig.head.label{i}, 'Nasion') || strcmp(orig.head.label{i}, 'RPA')
            shape.fid.pos(fidN,1:3) = orig.head.pos(i,:);
            shape.fid.label{fidN} = orig.head.label{i};
            fidN = fidN + 1;
          else
            shape.pos(posN,1:3) = orig.head.pos(i,:);
            shape.label{posN} = orig.head.label{i};
            posN = posN + 1;
          end
        end
        shape.coordsys = orig.head.coordsys;
      case 'dewar'
        fidN=1;
        posN=1;
        for i=1:size(orig.dewar.pos,1)
          if strcmp(orig.dewar.label{i}, 'LPA') || strcmp(orig.dewar.label{i}, 'Nasion') || strcmp(orig.dewar.label{i}, 'RPA')
            shape.fid.pos(fidN,1:3) = orig.dewar.pos(i,:);
            shape.fid.label{fidN} = orig.dewar.label{i};
            fidN = fidN + 1;
          else
            shape.pos(posN,1:3) = orig.dewar.pos(i,:);
            shape.label{posN} = orig.dewar.label{i};
            posN = posN + 1;
          end
        end
        shape.coordsys = orig.dewar.coordsys;
      otherwise
        ft_error('incorrect coordinates specified');
    end

  case {'ricoh_mrk', 'ricoh_ave', 'ricoh_con'}
    hdr = read_ricoh_header(filename);

    %% An exported file or an original one
    isexported = hdr.orig.digitize.info.done;

    %% Marker-coil positions
    mrk_pnt = hdr.orig.coregist.hpi;
    if any([mrk_pnt(:).meg_pos])
      mrk_pos = cat(1, mrk_pnt(1:end).meg_pos);
      mrk_label = transpose({mrk_pnt(1:end).label});
      sw_ind = [3 1 2];
      mrk_pos(1:3,:)= mrk_pos(sw_ind, :);
      mrk_pos = mrk_pos * 100; % unit: cm
      mrk_label(1:3,:)= mrk_label(sw_ind, :);
    else
      ft_error('No coil information found in the file');
    end

    %% Digitized points
    if ~isexported
      ft_info('The input file is an original one: only marker-coil positions are loaded\n');
      % The fiducial points are represented by the marker-coil positions.
      if any([mrk_pnt(:).meg_pos])
        shape.fid.pos = mrk_pos;   % unit: cm
        shape.fid.label = {'nas'; 'lpa'; 'rpa'; 'Marker4'; 'Marker5'};
      end
    else
      ft_info('The input file is a third-party-exported one including the digitized points\n');
      % All digitized points
      dig_pnt = hdr.orig.digitize.point;
      digitizer2meg = hdr.orig.digitize.info.digitizer2meg;
      R = digitizer2meg(1:3,1:3);
      T = digitizer2meg(1:3,4);
      % Transform to MEG coordinate:
      shape.pos = transpose( R * [[dig_pnt.x]; [dig_pnt.y]; [dig_pnt.z]] + repmat(T, 1, numel(dig_pnt))).*100; % unit: cm
      shape.label = transpose( deblank({dig_pnt(1:end).name}));
      % Fiducial points
      nas = find(strcmpi(shape.label, 'fidnz'));
      lpa = find(strcmpi(shape.label, 'fidt9'));
      rpa = find(strcmpi(shape.label, 'fidt10'));
      if ~isempty(nas) && ~isempty(lpa) && ~isempty(rpa)
        anatfid_pos = [ shape.pos(nas,:) ;  shape.pos(lpa,:) ; shape.pos(rpa,:) ];
        anatfid_label = {'nas'; 'lpa'; 'rpa'};
      end
      % HPIs
      HPI_1 = find(strcmpi(shape.label, 'HPI_1'));
      HPI_2 = find(strcmpi(shape.label, 'HPI_2'));
      HPI_3 = find(strcmpi(shape.label, 'HPI_3'));
      HPI_4 = find(strcmpi(shape.label, 'HPI_4'));
      HPI_5 = find(strcmpi(shape.label, 'HPI_5'));
      if ~isempty(HPI_1) && ~isempty(HPI_2) && ~isempty(HPI_3) && ~isempty(HPI_4) && ~isempty(HPI_5)
        HPI_pos = [ shape.pos(HPI_3,:) ;...
          shape.pos(HPI_1,:) ;...
          shape.pos(HPI_2,:) ;...
          shape.pos(HPI_4,:) ;...
          shape.pos(HPI_5,:) ];
        HPI_label = { 'HPI_3'; 'HPI_1'; 'HPI_2'; 'HPI_4'; 'HPI_5' };
      end
      shape.fid.pos = [ anatfid_pos; HPI_pos; mrk_pos ];
      shape.fid.label = [ anatfid_label; HPI_label; mrk_label ];
    end
    % 'cm' as a unit for 'pos':
    shape.unit = 'cm';

  case {'yokogawa_mrk', 'yokogawa_ave', 'yokogawa_con'}
    if ft_hastoolbox('yokogawa_meg_reader')
      hdr = read_yokogawa_header_new(filename);
      %% Marker-coil positions
      mrk_pnt = hdr.orig.coregist.hpi;

      % markers 1-3 identical to zero: try *.mrk file
      if ~any([mrk_pnt(:).meg_pos])
        ft_info('Reading marker-coil positions from a .mrk file\n');
        [p, f, x] = fileparts(filename);
        filename = fullfile(p, [f '.mrk']);
        if exist(filename, 'file')
          hdr_tmp = read_yokogawa_header_new(filename);
          mrk_pnt = hdr_tmp.orig.coregist.hpi;
        end
      end

      if any([mrk_pnt(:).meg_pos])
        mrk_pos = cat(1, mrk_pnt(1:end).meg_pos);
        mrk_label = transpose({mrk_pnt(1:end).label});
        sw_ind = [3 1 2];
        mrk_pos(1:3,:)= mrk_pos(sw_ind, :);
        mrk_pos = mrk_pos * 100; % unit: cm
        mrk_label(1:3,:)= mrk_label(sw_ind, :);
      else
        ft_error('No coil information found in the file');
      end

      %% An exported file or an original one
      isexported = hdr.orig.digitize.info.done;

      %% Digitized points
      if ~isexported
        ft_info('The input file is an original one: only marker-coil positions are loaded\n');
        % The fiducial points are represented by the marker-coil positions.
        if any([mrk_pnt(:).meg_pos])
          shape.fid.pos = mrk_pos;   % unit: cm
          shape.fid.label = {'nas'; 'lpa'; 'rpa'; 'Marker4'; 'Marker5'};
        end
      else
        ft_info('The input file is a third-party-exported one including the digitized points\n');
        % All digitized points
        dig_pnt = hdr.orig.digitize.point;
        digitizer2meg = hdr.orig.digitize.info.digitizer2meg;
        R = digitizer2meg(1:3,1:3);
        T = digitizer2meg(1:3,4);
        % Transform to MEG coordinate:
        shape.pos = transpose( R * [[dig_pnt.x]; [dig_pnt.y]; [dig_pnt.z]] + repmat(T, 1, numel(dig_pnt))).*100; % unit: cm
        shape.label = transpose( deblank({dig_pnt(1:end).name}));
        % Fiducial points
        nas = find(strcmpi(shape.label, 'fidnz'));
        lpa = find(strcmpi(shape.label, 'fidt9'));
        rpa = find(strcmpi(shape.label, 'fidt10'));
        if ~isempty(nas) && ~isempty(lpa) && ~isempty(rpa)
          anatfid_pos = [ shape.pos(nas,:) ;  shape.pos(lpa,:) ; shape.pos(rpa,:) ];
          anatfid_label = {'nas'; 'lpa'; 'rpa'};
        end
        % HPIs
        HPI_1 = find(strcmpi(shape.label, 'HPI_1'));
        HPI_2 = find(strcmpi(shape.label, 'HPI_2'));
        HPI_3 = find(strcmpi(shape.label, 'HPI_3'));
        HPI_4 = find(strcmpi(shape.label, 'HPI_4'));
        HPI_5 = find(strcmpi(shape.label, 'HPI_5'));
        if ~isempty(HPI_1) && ~isempty(HPI_2) && ~isempty(HPI_3) && ~isempty(HPI_4) && ~isempty(HPI_5)
          HPI_pos = [ shape.pos(HPI_3,:) ;...
            shape.pos(HPI_1,:) ;...
            shape.pos(HPI_2,:) ;...
            shape.pos(HPI_4,:) ;...
            shape.pos(HPI_5,:) ];
          HPI_label = { 'HPI_3'; 'HPI_1'; 'HPI_2'; 'HPI_4'; 'HPI_5' };
        end
        shape.fid.pos = [ anatfid_pos; HPI_pos; mrk_pos ];
        shape.fid.label = [ anatfid_label; HPI_label; mrk_label ];
      end
      % 'cm' as a unit for 'pos':
      shape.unit = 'cm';

    else  % the case that "yokogawa_meg_reader" is not available
      hdr = read_yokogawa_header(filename);
      marker = hdr.orig.matching_info.marker;
      % markers 1-3 identical to zero: try *.mrk file
      if ~any([marker(:).meg_pos])
        [p, f, x] = fileparts(filename);
        filename = fullfile(p, [f '.mrk']);
        if exist(filename, 'file')
          hdr = read_yokogawa_header(filename);
          marker = hdr.orig.matching_info.marker;
        end
      end

      % non zero markers 1-3
      if any([marker(:).meg_pos])
        shape.fid.pos = cat(1, marker(1:5).meg_pos);
        sw_ind = [3 1 2];
        shape.fid.pos(1:3,:)= shape.fid.pos(sw_ind, :);
        shape.fid.label = {'nas'; 'lpa'; 'rpa'; 'Marker4'; 'Marker5'};
      else
        ft_error('no coil information found in Yokogawa file');
      end

      % convert to the units of the grad, the desired default for yokogawa is centimeter.
      shape = ft_convert_units(shape, 'cm');
    end

  case 'yokogawa_raw'
    if ft_hastoolbox('yokogawa_meg_reader')
      hdr = read_yokogawa_header_new(filename);
      marker = hdr.orig.coregist.hpi;
    else
      hdr = read_yokogawa_header(filename);
      marker = hdr.orig.matching_info.marker;
    end

    % markers 1-3 identical to zero: try *.mrk file
    if ~any([marker(:).meg_pos])
      [p, f, x] = fileparts(filename);
      filename = fullfile(p, [f '.mrk']);
      if exist(filename, 'file')
        if ft_hastoolbox('yokogawa_meg_reader')
          hdr = read_yokogawa_header_new(filename);
          marker = hdr.orig.coregist.hpi;
        else
          hdr = read_yokogawa_header(filename);
          marker = hdr.orig.matching_info.marker;
        end
      end
    end

    % non zero markers 1-3
    if any([marker(:).meg_pos])
      shape.fid.pos = cat(1, marker(1:5).meg_pos);
      sw_ind = [3 1 2];
      shape.fid.pos(1:3,:)= shape.fid.pos(sw_ind, :);
      shape.fid.label = {'nas'; 'lpa'; 'rpa'; 'Marker4'; 'Marker5'};
    else
      ft_error('no coil information found in Yokogawa file');
    end

    % convert to the units of the grad, the desired default for yokogawa is centimeter.
    shape = ft_convert_units(shape, 'cm');

    %  case {'yokogawa_mrk', 'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw' }
    %    if ft_hastoolbox('yokogawa_meg_reader')
    %      hdr = read_yokogawa_header_new(filename);
    %      marker = hdr.orig.coregist.hpi;
    %    else
    %      hdr = read_yokogawa_header(filename);
    %      marker = hdr.orig.matching_info.marker;
    %    end
    %
    %    % markers 1-3 identical to zero: try *.mrk file
    %    if ~any([marker(:).meg_pos])
    %      [p, f, x] = fileparts(filename);
    %      filename = fullfile(p, [f '.mrk']);
    %      if exist(filename, 'file')
    %        if ft_hastoolbox('yokogawa_meg_reader')
    %          hdr = read_yokogawa_header_new(filename);
    %          marker = hdr.orig.coregist.hpi;
    %        else
    %          hdr = read_yokogawa_header(filename);
    %          marker = hdr.orig.matching_info.marker;
    %        end
    %      end
    %    end
    %
    %    % non zero markers 1-3
    %    if any([marker(:).meg_pos])
    %      shape.fid.pos = cat(1, marker(1:5).meg_pos);
    %      sw_ind = [3 1 2];
    %      shape.fid.pos(1:3,:)= shape.fid.pos(sw_ind, :);
    %      shape.fid.label = {'nas'; 'lpa'; 'rpa'; 'Marker4'; 'Marker5'};
    %    else
    %      ft_error('no coil information found in Yokogawa file');
    %    end
    %
    %    % convert to the units of the grad, the desired default for yokogawa is centimeter.
    %    shape = ft_convert_units(shape, 'cm');

  case 'yokogawa_coregis'
    in_str = textread(filename, '%s');
    nr_items = size(in_str,1);
    ind = 1;
    coil_ind = 1;
    shape.fid.pos = [];
    shape.fid.label = {};
    while ind < nr_items
      if strcmp(in_str{ind},'MEG:x=')
        shape.fid.pos = [shape.fid.pos; str2num(strtok(in_str{ind+1},[',','['])) ...
          str2num(strtok(in_str{ind+3},[',','['])) str2num(strtok(in_str{ind+5},[',','[']))];
        shape.fid.label = [shape.fid.label ; ['Marker',num2str(coil_ind)]];
        coil_ind = coil_ind + 1;
        ind = ind + 6;
      else
        ind = ind +1;
      end
    end
    if size(shape.fid.label,1) ~= 5
      ft_error('Wrong number of coils');
    end

    sw_ind = [3 1 2];

    shape.fid.pos(1:3,:)= shape.fid.pos(sw_ind, :);
    shape.fid.label(1:3)= {'nas', 'lpa', 'rpa'};

  case 'yokogawa_hsp'
    fid = fopen_or_error(filename, 'rt');

    fidstart = false;
    hspstart = false;

    % try to locate the fiducial positions
    while ~fidstart && ~feof(fid)
      line = fgetl(fid);
      if ~isempty(strmatch('//Position of fiducials', line))
        fidstart = true;
      end
    end
    if fidstart
      line_xpos = fgetl(fid);
      line_ypos = fgetl(fid);
      line_yneg = fgetl(fid);
      xpos = sscanf(line_xpos(3:end), '%f');
      ypos = sscanf(line_ypos(3:end), '%f');
      yneg = sscanf(line_yneg(3:end), '%f');
      shape.fid.pos = [
        xpos(:)'
        ypos(:)'
        yneg(:)'
        ];
      shape.fid.label = {
        'X+'
        'Y+'
        'Y-'
        };
    end

    % try to locate the fiducial positions
    while ~hspstart && ~feof(fid)
      line = fgetl(fid);
      if ~isempty(strmatch('//No of rows', line))
        hspstart = true;
      end
    end
    if hspstart
      line = fgetl(fid);
      siz = sscanf(line, '%f');
      shape.pos = zeros(siz(:)');
      for i=1:siz(1)
        line = fgetl(fid);
        shape.pos(i,:) = sscanf(line, '%f');
      end
    end

    fclose(fid);

  case 'yorkinstruments_hdf5'
    acquisition='default';
    try
      shape.pos=transpose(h5read(filename,  '/geometry/head_shape/head_shape'));
    catch
      error('Headshape data not found.');
    end
    shape.unit='mm';
    temp=h5info(filename,  '/geometry/fiducials/');
    Nfids=length(temp.Groups);
    for i=1:Nfids
      [null, shape.fid.label{i}, null]= fileparts(temp.Groups(i).Name);
      shape.fid.pos(i,1:3)=h5read(filename, strcat('/geometry/fiducials/',shape.fid.label{i} ,'/location'));
    end
    if isempty(coordsys)
      coordsys='dewar'
    end
    if strcmp(coordsys,'dewar')
      try
        tCCStoMegscanScs = h5read(filename,[strcat('/acquisitions/',char(string(acquisition))) '/ccs_to_scs_transform']);
        T = maketform('affine',tCCStoMegscanScs);
        shape.pos=tforminv(T,shape.pos(:,1),shape.pos(:,2),shape.pos(:,3));
        shape.fid.pos=tforminv(T,shape.fid.pos(:,1),shape.fid.pos(:,2),shape.fid.pos(:,3));
      catch
        error('No head to dewar transform available in hdf5 file');
      end
    end

  case 'ply'
    [vert, face] = read_ply(filename);
    shape.pos = [vert.x vert.y vert.z];
    if isfield(vert, 'red') && isfield(vert, 'green') && isfield(vert, 'blue')
      shape.color = double([vert.red vert.green vert.blue])/255;
    end
    switch size(face,2)
      case 3
        shape.tri = face;
      case 4
        shape.tet = face;
      case 8
        shape.hex = face;
    end

  case 'polhemus_fil'
    [shape.fid.pos, shape.pos, shape.fid.label] = read_polhemus_fil(filename, 0);

  case 'polhemus_pos'
    [shape.fid.pos, shape.pos, shape.fid.label] = read_ctf_pos(filename);

  case 'spmeeg_mat'
    tmp = load(filename);
    if isfield(tmp.D, 'fiducials') && ~isempty(tmp.D.fiducials)
      shape = tmp.D.fiducials;
    else
      ft_error('no headshape found in SPM EEG file');
    end

  case 'matlab'
    % determine which variables are contained in the file
    tmp = load(filename);

    if isfield(tmp, 'shape')
      shape = tmp.shape;
    elseif isfield(tmp, 'headshape')
      shape = tmp.headshape;
    elseif isfield(tmp, 'surface')
      shape = tmp.surface;
    elseif isfield(tmp, 'bnd')
      % the variable in the file is most likely a precomputed triangulation of some sort
      shape = tmp.bnd;
    elseif isfield(tmp, 'mesh')
      % the variable in the file is most likely a precomputed triangulation of some sort
      shape = tmp.mesh;
    elseif isfield(tmp, 'elec')
      % represent the electrodes as headshape
      tmp.elec        = ft_datatype_sens(tmp.elec);
      shape.fid.pos   = tmp.elec.chanpos;
      shape.fid.label = tmp.elec.label;
    elseif isfield(tmp, 'Vertices')
      % this applies to BrainStorm cortical meshes
      shape.pos = tmp.Vertices;
      % copy some optional fields over with a new name
      shape = copyfields(tmp, shape, {'Faces', 'Curvature', 'SulciMap'});
      shape = renamefields(shape, {'Faces', 'Curvature', 'SulciMap'}, {'tri', 'curv', 'sulc'});
    elseif numel(fieldnames(tmp))==1
      fn = fieldnames(tmp);
      shape = tmp.(fn{1});
      % check that it has vertices and triangles
      assert(isfield(shape, 'pos') && isfield(shape, 'tri'), 'no headshape found in MATLAB file')
    else
      ft_error('no headshape found in MATLAB file');
    end

  case {'freesurfer_triangle_binary', 'freesurfer_quadrangle'}
    % the freesurfer toolbox is required for this
    ft_hastoolbox('freesurfer', 1);

    [pos, tri] = read_surf(filename);

    if min(tri(:)) == 0
      % start counting from 1
      tri = tri + 1;
    end
    shape.pos = pos;
    shape.tri = tri;

    % for the left and right
    [path,name,ext] = fileparts(filename);

    if strcmp(ext, '.inflated') % does the shift only for inflated surface
      if strcmp(name, 'lh')
        % assume freesurfer inflated mesh in mm, mni space
        % move the mesh a bit to the left, to avoid overlap with the right
        % hemisphere
        shape.pos(:,1) = shape.pos(:,1) - max(shape.pos(:,1)) - 10;

      elseif strcmp(name, 'rh')
        % id.
        % move the mesh a bit to the right, to avoid overlap with the left
        % hemisphere
        shape.pos(:,1) = shape.pos(:,1) - min(shape.pos(:,1)) + 10;
      end
    end

    if exist(fullfile(path, [name,'.sulc']), 'file'), shape.sulc = read_curv(fullfile(path, [name,'.sulc'])); end
    if exist(fullfile(path, [name,'.curv']), 'file'), shape.curv = read_curv(fullfile(path, [name,'.curv'])); end
    if exist(fullfile(path, [name,'.area']), 'file'), shape.area = read_curv(fullfile(path, [name,'.area'])); end
    if exist(fullfile(path, [name,'.thickness']), 'file'), shape.thickness = read_curv(fullfile(path, [name,'.thickness'])); end

  case 'stl'
    [pos, tri, nrm] = read_stl(filename);
    shape.pos = pos;
    shape.tri = tri;

  case 'obj'
    ft_hastoolbox('wavefront', 1);
    % Only tested for structure.io .obj thus far
    [pos, tri, texture, textureIdx] = read_obj_new(filename);

    % check if the texture is defined per vertex, in which case the texture can be refined below
    if size(texture, 1)==size(pos, 1)
      texture_per_vert = true;
    else
      texture_per_vert = false;
    end

    % remove the triangles with 0's first
    allzeros = sum(tri==0,2)==3;
    tri(allzeros, :)        = [];
    textureIdx(allzeros, :) = [];

    % check whether all vertices belong to a triangle. If not, then prune the vertices and keep the faces consistent.
    utriIdx = unique(tri(:));
    remove  = setdiff((1:size(pos, 1))', utriIdx);
    if ~isempty(remove)
      [pos, tri] = remove_vertices(pos, tri, remove);
      if texture_per_vert
        % also remove the removed vertices from the texture
        texture(remove, :) = [];
      end
    end

    if hasimage
      % there is an image with color information

      if texture_per_vert
        % refine the mesh and texture mapping to increase the resolution
        for i=1:refine_
          [pos, tri, texture] = refine(pos, tri, 'banks', texture);
        end

        % do the texture to color mapping
        picture = imread(image);
        [sy, sx, sz] = size(picture);

        color = zeros(size(pos, 1), 3);
        for i = 1:size(pos, 1)
          x = floor((1-texture(i,2))*sx);
          y = 1+floor(texture(i,1)*sy);
          color(i,1:3) = picture(x,y,1:3);
        end

      else
        % do the texture to color mapping based on the textureIdx
        picture      = flip(imread(image),1);
        [sy, sx, sz] = size(picture);
        picture      = reshape(picture, sy*sx, sz);

        % make image 3D if grayscale
        if sz == 1
          picture = repmat(picture, 1, 3);
        end
        [dum, ix]  = unique(tri);
        textureIdx = textureIdx(ix);

        % get the indices into the image
        x     = abs(round(texture(:,1)*(sx-1)))+1;
        y     = abs(round(texture(:,2)*(sy-1)))+1;
        xy    = sub2ind([sy sx], y, x);
        sel   = xy(textureIdx);
        color = double(picture(sel,:))/255;
      end

    elseif size(pos, 2)==6
      % there is no separate image, but the vertices also contain RGB colors
      color = pos(:, 4:6);
      pos   = pos(:, 1:3);

    else
      % there is no color information
      color = [];
    end

    shape.pos = pos - repmat(mean(pos,1), [size(pos, 1),1]); % centering vertices
    shape.tri = tri;

    if ~isempty(color)
      if range(color(:)) > 1
        % color should be specified between 0 and 1
        shape.color = color./255;
      else
        shape.color = color;
      end
    end

  case 'vtk'
    [pos, tri, attr] = read_vtk(filename);
    shape.pos = pos;
    if ~isempty(tri)
      shape.tri = tri;
    end
    if ~isempty(attr)
      shape.data = attr;
    end

  case 'vtk_xml'
    data = read_vtk_xml(filename);
    shape.orig = data;
    shape.pos  = data.Points;
    if isfield(data, 'Lines')
      shape.line = data.Lines;
    end

  case 'mrtrix_tck'
    ft_hastoolbox('mrtrix', 1);
    shape = read_tck(filename);

  case 'trackvis_trk'
    shape = read_trk(filename);

  case 'off'
    [pos, plc] = read_off(filename);
    shape.pos  = pos;
    shape.tri  = plc;

  case 'mne_tri'
    % FIXME this should be implemented, consistent with ft_write_headshape
    keyboard

  case 'mne_pos'
    % FIXME this should be implemented, consistent with ft_write_headshape
    keyboard

  case 'netmeg'
    hdr = ft_read_header(filename);
    if isfield(hdr.orig, 'headshapedata')
      shape.pos = hdr.orig.Var.headshapedata;
    else
      ft_error('the NetMEG file "%s" does not contain headshape data', filename);
    end

  case 'vista'
    ft_hastoolbox('simbio', 1);
    [nodes,elements,labels] = read_vista_mesh(filename);
    shape.pos     = nodes;
    if size(elements,2)==8
      shape.hex     = elements;
    elseif size(elements,2)==4
      shape.tet = elements;
    else
      ft_error('unknown elements format')
    end
    % representation of data is compatible with ft_datatype_parcellation
    shape.tissue = zeros(size(labels));
    numlabels = size(unique(labels),1);
    shape.tissuelabel = {};
    for i = 1:numlabels
      ulabel = unique(labels);
      shape.tissue(labels == ulabel(i)) = i;
      shape.tissuelabel{i} = num2str(ulabel(i));
    end

  case 'tet'
    % the toolbox from Gabriel Peyre has a function for this
    ft_hastoolbox('toolbox_graph', 1);
    [vertexline, face] = read_tet(filename);
    %     'vertex' is a '3 x nb.vert' array specifying the position of the vertices.
    %     'face' is a '4 x nb.face' array specifying the connectivity of the tet mesh.
    shape.pos = vertexline';
    shape.tet = face';

  case 'tetgen_ele'
    % reads in the tetgen format and rearranges according to FT conventions
    % tetgen files also return a 'faces' field, which is not used here
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, f); % without the extension
    IMPORT = importdata([filename '.ele'],' ',1);
    shape.tet = IMPORT.data(:,2:5);
    if size(IMPORT.data,2)==6
      labels = IMPORT.data(:,6);
      % representation of tissue type is compatible with ft_datatype_parcellation
      numlabels    = size(unique(labels),1);
      ulabel       = unique(labels);
      shape.tissue = zeros(size(labels));
      shape.tissuelabel = {};
      for i = 1:numlabels
        shape.tissue(labels == ulabel(i)) = i;
        shape.tissuelabel{i} = num2str(ulabel(i));
      end
    end
    IMPORT = importdata([filename '.node'],' ',1);
    shape.pos = IMPORT.data(:,2:4);

  case 'brainsuite_dfs'
    % this requires the readdfs function from the BrainSuite MATLAB utilities
    ft_hastoolbox('brainsuite', 1);

    dfs = readdfs(filename);
    % these are expressed in MRI dimensions
    shape.pos  = dfs.vertices;
    shape.tri  = dfs.faces;
    shape.unit = 'unkown';

    % the filename is something like 2467264.right.mid.cortex.svreg.dfs
    % whereas the corresponding MRI is 2467264.nii and might be gzipped
    [p, f, x] = fileparts(filename);
    while ~isempty(x)
      [junk, f, x] = fileparts(f);
    end

    if exist(fullfile(p, [f '.nii']), 'file')
      fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
      mri = ft_read_mri(fullfile(p, [f '.nii']));
      transform = eye(4);
      transform(1:3,4) = mri.transform(1:3,4); % only use the translation
      shape.pos  = ft_warp_apply(transform, shape.pos);
      shape.unit = mri.unit;
    elseif exist(fullfile(p, [f '.nii.gz']), 'file')
      fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
      mri = ft_read_mri(fullfile(p, [f '.nii.gz']));
      transform = eye(4);
      transform(1:3,4) = mri.transform(1:3,4); % only use the translation
      shape.pos  = ft_warp_apply(transform, shape.pos);
      shape.unit = mri.unit;
    else
      ft_warning('could not find accompanying MRI file, returning vertices in voxel coordinates');
    end

  case 'brainvisa_mesh'
    % this requires the loadmesh function from the BrainVISA MATLAB utilities
    ft_hastoolbox('brainvisa', 1);
    [shape.pos, shape.tri, shape.nrm] = loadmesh(filename);
    shape.tri = shape.tri + 1; % they should be 1-offset, not 0-offset
    shape.unit = 'unkown';

    if exist([filename '.minf'], 'file')
      minffid = fopen_or_error([filename '.minf']);
      hdr=fgetl(minffid);
      tfm_idx = strfind(hdr,'''transformations'':') + 21;
      transform = sscanf(hdr(tfm_idx:end),'%f,',[4 4])';
      fclose(minffid);
      if ~isempty(transform)
        shape.pos = ft_warp_apply(transform, shape.pos);
        shape = rmfield(shape, 'unit'); % it will be determined later on, based on the size
      end
    end

    if isempty(transform)
      % the transformation was not present in the minf file, try to get it from the MRI

      % the filename is something like subject01_Rwhite_inflated_4d.mesh
      % and it is accompanied by subject01.nii
      [p, f, x] = fileparts(filename);
      f = tokenize(f, '_');
      f = f{1};

      if exist(fullfile(p, [f '.nii']), 'file')
        fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii']));
        mri = ft_read_mri(fullfile(p, [f '.nii']));
        shape.pos  = ft_warp_apply(mri.transform, shape.pos);
        shape.unit = mri.unit;
        transform = true; % used for feedback
      elseif exist(fullfile(p, [f '.nii.gz']), 'file')
        fprintf('reading accompanying MRI file "%s"\n', fullfile(p, [f '.nii.gz']));
        mri = ft_read_mri(fullfile(p, [f '.nii.gz']));
        shape.pos  = ft_warp_apply(mri.transform, shape.pos);
        shape.unit = mri.unit;
        transform = true; % used for feedback
      end
    end

    if isempty(transform)
      ft_warning('cound not determine the coordinate transformation, returning vertices in voxel coordinates');
    end

  case 'brainvoyager_srf'
    [pos, tri, srf] = read_bv_srf(filename);
    shape.pos = pos;
    shape.tri = tri;

    % FIXME add details from srf if possible
    % FIXME do transform
    % FIXME remove vertices that are not in a triangle
    % FIXME add unit

  case 'besa_sfp'
    [lab, pos] = read_besa_sfp(filename, 0);
    shape.pos = pos;

    % assume that all non-'headshape' points are fiducial markers
    hs = strmatch('headshape', lab);
    lab(hs) = [];
    pos(hs, :) = [];
    shape.fid.label = lab;
    shape.fid.pos = pos;

  case 'asa_elc'
    elec = ft_read_sens(filename);

    shape.fid.pos   = elec.chanpos;
    shape.fid.label = elec.label;

    npos = read_ini(filename, 'NumberHeadShapePoints=', '%d');
    if ~isempty(npos) && npos>0
      origunit = read_ini(filename, 'UnitHeadShapePoints', '%s', 1);
      pos = read_ini(filename, 'HeadShapePoints', '%f', npos, ':');
      pos = ft_scalingfactor(origunit, 'mm')*pos;

      shape.pos = pos;
    end

  case 'neuromag_mesh'
    fid = fopen_or_error(filename, 'rt');
    npos = fscanf(fid, '%d', 1);
    pos = fscanf(fid, '%f', [6 npos])';
    ntri = fscanf(fid, '%d', 1);
    tri = fscanf(fid, '%d', [3 ntri])';
    fclose(fid);

    shape.pos = pos(:,1:3); % vertex positions
    shape.nrm = pos(:,4:6); % vertex normals
    shape.tri = tri;

  case 'duneuro_dgf'
    lines = readlines(filename);
    lines = cellstr(lines);
    % remove comments
    sel = startsWith(lines, '#');
    lines(sel) = [];
    % remove empty lines
    sel = cellfun(@isempty, lines);
    lines(sel) = [];

    vertexline = find(strcmp(lines, 'Vertex'));
    cubeline   = find(strcmp(lines, 'Cube'));
    paramline  = find(startsWith(lines, 'parameters'));

    % parse this line to determine the number of parameters, and thereby the number of additional columns
    numparam = sscanf(lines{paramline}, 'parameters %d');

    npos = cubeline - vertexline - 1;
    shape.pos = nan(npos, 3);
    for i=1:npos
      shape.pos(i,:) = str2num(lines{vertexline+i});
    end

    nhex = length(lines)-paramline;
    shape.hex = nan(nhex, 8+numparam);
    for i=1:nhex
      shape.hex(i,:) = str2num(lines{paramline+i});
    end

    if numparam==1
      % assume that this is the tissue class
      shape.tissue = shape.hex(:,9);
      shape.tissue = shape.tissue + 1; % this should be one-offset
    end

    % remove the parameter columns
    shape.hex = shape.hex(:,1:8);
    shape.hex = shape.hex + 1; % this should be one-offset

  case {'gmsh_ascii' 'gmsh_binary'}
    % use the SimNIBS reader, this does not read all gmsh properties/tags
    ft_hastoolbox('simnibs', 1);
    shape = mesh_load_gmsh4(filename);
    shape = fixpos(shape);

    % remove empty fields
    fn = fieldnames(shape);
    for i=1:numel(fn)
      if isempty(shape.(fn{i}))
        shape = rmfield(shape, fn{i});
      end
    end

    hastri = isfield(shape, 'tri');
    hastet = isfield(shape, 'tet');
    
    if isempty(meshtype)
      if hastri && hastet
        ft_warning('mesh has both tri and tet, returning tet');
        meshtype = 'tet';
      elseif hastet
        meshtype = 'tet';
      elseif hastri
        meshtype = 'tri';
      end
    end

    if isfield(shape, 'triangle_regions')
      if strcmp(meshtype, 'tri')
        % use this as the tissue type when keeping triangles
        shape.tissue = shape.triangle_regions;
      end
      shape = rmfield(shape, 'triangle_regions');
    end

    if isfield(shape, 'tetrahedron_regions')
      if strcmp(meshtype, 'tet')
        % use this as the tissue type when keeping tetrahedrons
        shape.tissue = shape.tetrahedron_regions;
      end
      shape = rmfield(shape, 'tetrahedron_regions');
    end

    if isfield(shape, 'tissue')
      if all(shape.tissue(:)>999)
        % this is the case for surface meshes with triangles
        shape.tissue = shape.tissue-1000;
      end
      [labels, values, rgba] = simnibs_labels;
      shape.tissuelabel = labels;
      shape.rgba = rgba; % this is consistent with FT_READ_ATLAS
    end

  case 'gmsh_binary_v1'
    % use Jan-Mathijs' reader, this only works for binary files but does read all gmsh properties/tags
    [nodes, elements] = read_gmsh_binary(filename);
    shape.pos = nodes.nodes(nodes.indx, :);

    % this file format may contain a mixture of differently shaped elements
    fnames = fieldnames(elements);
    for k=1:numel(fnames)
      el = elements.(fnames{k});
      switch fnames{k}
        case 'lines'
          shape.line = el;
        case 'triangles'
          shape.tri = el;
        case 'tetrahedra'
          shape.tet = el;
        case 'hexahedra'
          shape.hex = el;
        case 'lines_tag'
          shape.tag_line = el;
        case 'triangles_tag'
          shape.tag_tri = el;
        case 'tetrahedra_tag'
          shape.tag_tet = el;
        case 'hexahedra_tag'
          shape.tag_hex = el;
        otherwise
          ft_warning('skipping element field %s', fnames{k});
      end
    end

  case {'neurojson_jmesh' 'neurojson_bmesh'}
    % see https://github.com/NeuroJSON/jmesh/blob/master/JMesh_specification.md
    ft_hastoolbox('jsonlab', 1);

    extraopt = jsonopt('jmeshopt', {}, varargin2struct(varargin{:}));
    if strcmp(fileformat, 'neurojson_bmesh')
      jmesh = loadbj(filename, extraopt{:});
    else
      jmesh = loadjson(filename, extraopt{:});
    end

    % jmesh metadata
    if(isfield(jmesh, encodevarname('_DataInfo_')))
      shape.info = jmesh.(encodevarname('_DataInfo_'));
    end

    % node data
    if(isfield(jmesh, 'MeshVertex3'))
      shape.pos  = jmesh.MeshVertex3;
    elseif(isfield(jmesh, 'MeshNode'))
      shape.pos  = jmesh.MeshNode;
    else
      ft_error('no vertex positions found');
    end

    % this applies when "Failed to decode embedded JData annotations, return raw JSON data"
    if ~isnumeric(shape.pos)
      ft_error('cannot read file "%s"', filename);
    end

    % extract node label if present
    if(isfield(shape, 'pos') && isstruct(shape.pos))
      if(isfield(shape.pos, 'Properties'))
        shape = copyfields(shape, shape.pos.Properties, {'Color','Normal','Size','Tag','Value','Texture'});
      end
      if(isfield(shape.pos, 'Data'))
        shape.pos = shape.pos.Data;
      end
    end

    % surface data
    if(isfield(jmesh, 'MeshTri3'))
      shape.tri  = jmesh.MeshTri3;
    elseif(isfield(jmesh, 'MeshSurf'))
      shape.tri  = jmesh.MeshSurf;
    end

    % extract surface label if present
    if(isfield(shape, 'tri') && isstruct(shape.tri))
      if(isfield(shape.tri, 'Properties'))
        shape = copyfields(shape.tri.Properties, shape, {'Color','Normal','Size','Tag','Value','Texture'});
      end
      if(isfield(shape.tri, 'Data'))
        shape.tri = shape.tri.Data;
      end
    end

    % tet element data
    if(isfield(jmesh, 'MeshTet4'))
      shape.tet = jmesh.MeshTet4;
    elseif(isfield(jmesh, 'MeshElem'))
      shape.tet = jmesh.MeshElem;
    end

    % extract tet label if present
    if(isfield(shape, 'tet') && isstruct(shape.tet))
      if(isfield(shape.tet, 'Properties'))
        shape = copyfields(shape.tet.Properties, shape, {'Color','Normal','Size','Tag','Value','Texture'});
      end
      if(isfield(shape.tet, 'Data'))
        shape.tet = shape.tet.Data;
      end
    end

    % hex element data
    if(isfield(jmesh, 'MeshHex8'))
      shape.hex = jmesh.MeshHex8;
    end

    % extract hex label if present
    if(isfield(shape, 'hex') && isstruct(shape.hex))
      if(isfield(shape.hex, 'Properties'))
        shape = copyfields(shape.hex.Properties, shape, {'Color','Normal','Size','Tag','Value','Texture'});
      end
      if(isfield(shape.hex, 'Data'))
        shape.hex = shape.hex.Data;
      end
    end

    % line segment data
    if(isfield(jmesh, 'MeshEdge'))
      shape.line  = jmesh.MeshEdge;
    end

    % extract line label if present
    if(isfield(shape, 'line') && isstruct(shape.line))
      if(isfield(shape.hex, 'Properties'))
        shape = copyfields(shape.line.Properties, shape, {'Color','Normal','Size','Tag','Value','Texture'});
      end
      if(isfield(shape.line, 'Data'))
        shape.line = shape.line.Data;
      end
    end

    % rename upper-case fields to lowe-case
    fn = intersect(fieldnames(shape), {'Color','Normal','Size','Tag','Value','Texture'});
    for i=1:numel(fn)
      tmp = shape.(fn{i});
      shape = rmfield(shape, fn{i});
      % assume that the number of pos/tri/tet/hex is larger than the number of values
      if size(tmp,2)>size(tmp,1)
        tmp = tmp';
      end
      shape.(lower(fn{i})) = tmp;
    end

  otherwise
    % try reading it from an electrode of volume conduction model file
    success = false;

    if ~success
      % try reading it as electrode positions and treat those as fiducials
      try
        elec = ft_read_sens(filename, 'senstype', 'eeg');
        if ~ft_senstype(elec, 'eeg')
          ft_error('headshape information can not be read from MEG gradiometer file');
        else
          shape.fid.pos   = elec.chanpos;
          shape.fid.label = elec.label;
          success = 1;
        end
      catch
        success = false;
      end % try
    end

    if ~success
      % try reading it as volume conductor
      % and treat the skin surface as headshape
      try
        headmodel = ft_read_headmodel(filename);
        if ~ft_headmodeltype(headmodel, 'bem')
          ft_error('skin surface can only be extracted from boundary element model');
        else
          if ~isfield(headmodel, 'skin')
            headmodel.skin = find_outermost_boundary(headmodel.bnd);
          end
          shape.pos = headmodel.bnd(headmodel.skin).pos;
          shape.tri = headmodel.bnd(headmodel.skin).tri; % also return the triangulation
          success = 1;
        end
      catch
        success = false;
      end % try
    end

    if ~success
      ft_error('unknown fileformat "%s" for head shape information', fileformat);
    end
end % switch fileformat

if isfield(shape, 'label')
  % ensure that it is a column
  shape.label = shape.label(:);
end

if isfield(shape, 'fid') && isfield(shape.fid, 'label')
  % ensure that it is a column
  shape.fid.label = shape.fid.label(:);
end

% this will add the units to the head shape and optionally convert
if ~isempty(unit)
  shape = ft_convert_units(shape, unit);
else
  try
    % ft_determine_units will fail for triangle-only gifties.
    shape = ft_determine_units(shape);
  end
end

% ensure that vertex positions are given in pos, not in pnt
shape = fixpos(shape);

% ensure that the numerical arrays are represented in double precision and not as integers
shape = ft_struct2double(shape);

% determine which types of meshes to keep if there are multiple
hastri = isfield(shape, 'tri');
hastet = isfield(shape, 'tet');
hashex = isfield(shape, 'hex');

if (hastri+hastet+hashex)>1
  % there are multiple types of meshes
  % this also allows for specifications like 'tri+tet'
  if  hastri && ~strcmp(meshtype, 'tri')
    if isempty(meshtype)
      ft_warning('removing surface mesh (tri), use the ''meshtype'' option to keep it')
    end
    shape = removefields(shape, 'tri');
  elseif hastet && ~strcmp(meshtype, 'tet')
    if isempty(meshtype)
      ft_warning('removing tetrahedral mesh (tet), use the ''meshtype'' option to keep it')
    end
    shape = removefields(shape, 'tet');
  elseif hashex && ~strcmp(meshtype, 'hex')
    if isempty(meshtype)
      ft_warning('removing hexahedral mesh (hex), use the ''meshtype'' option to keep it')
    end
    shape = removefields(shape, 'hex');
  end

end % if multiple types of meshes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [labels, values, rgba] = simnibs_labels
% the following is from the final_tissues_LUT.txt and the .msh.opt files
% the values should be offset with 1000 in case of triangles

values = [
  1
  2
  3
  4
  5
  6
  7
  8
  9
  100
  500
  ];

labels = {
  'WM'
  'GM'
  'CSF'
  'Bone'
  'Scalp'
  'Eye_balls'
  'Compact_bone'
  'Spongy_bone'
  'Blood'
  'Muscle'
  'Electrode'
  'Saline_or_gel'
  };

rgba = [
  230 230 230 255
  129 129 129 255
  104 163 255 255
  255 239 179 255
  255 166 133 255
  255 240 0 255
  255 239 179 255
  255 138 57 255
  0 65 142 255
  0 118 14 255
  37 79 255 255
  103 255 226 255
  ];

