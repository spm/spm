function [shape] = read_headshape(filename, varargin)

% READ_HEADSHAPE reads the fiducials and/or the measured headshape
% from a variety of files (like CTF and Polhemus). The headshape and
% fiducials can for example be used for coregistration.
%
% Use as
%   [shape] = read_headshape(filename)
%
% See also READ_VOL, READ_SENS

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: read_headshape.m,v $
% Revision 1.3  2008/04/16 08:04:03  roboos
% allow headshape to be extracted from BEM volume conduction model
%
% Revision 1.2  2008/04/14 20:52:11  roboos
% ensure consistent output for all  file formats (thanks to Vladimir)
% added convert_units
%
% Revision 1.1  2008/04/11 12:04:55  roboos
% new impoementation, required for clean interface towards SPM
%

% test whether the file exists
if ~exist(filename)
  error(sprintf('file ''%s'' does not exist', filename));
end

% get the options
fileformat = keyval('fileformat',  varargin);

if isempty(fileformat)
  fileformat = filetype(filename);
end

% start with an empty structure
shape           = [];
shape.pnt       = [];
shape.fid.pnt   = [];
shape.fid.label = {};

switch fileformat
  case {'ctf_ds', 'ctf_hc', 'ctf_meg4', 'ctf_res4'}
    [p, f, x] = fileparts(filename);
    if strcmp(fileformat, 'ctf_ds')
      filename = fullfile(p, f, [f '.hc']);
    elseif strcmp(fileformat, 'ctf_meg4')
      filename = fullfile(p, [f '.hc']);
    elseif strcmp(fileformat, 'ctf_res4')
      filename = fullfile(p, [f '.hc']);
    end

    orig = read_ctf_hc(filename);
    shape.fid.pnt = cell2mat(struct2cell(orig.head));
    shape.fid.label = fieldnames(orig.head);

  case 'ctf_shape'
    orig = read_ctf_shape(filename);
    shape.pnt = orig.pnt;
    shape.fid.label = {'NASION', 'LEFT_EAR', 'RIGHT_EAR'};
    for i = 1:numel(shape.fid.label)
      shape.fid.pnt = cat(1, shape.fid.pnt, ...
        getfield(orig.MRI_Info, shape.fid.label{i}));
    end

  case '4d_hs'
    shape.pnt = read_bti_hs(filename);

  case 'polhemus_fil'
    [shape.fid.pnt, shape.pnt, shape.fid.label] = read_polhemus_fil(filename, 0);

  case 'matlab'
    tmp = load(filename);
    if isfield(tmp, 'shape')
      shape = tmp.shape;
    elseif isfield(tmp, 'elec')
        shape.fid.pnt   = tmp.elec.pnt;
        shape.fid.label = tmp.elec.label;
    else
      error('no headshape found in Matlab file');
    end

  otherwise

    success = 0;
    if ~success
      % try reading it as electrode positions
      % and treat those as fiducials
      try
        elec = read_sens(filename);
        if ~senstype(elec, 'eeg')
          error('headshape information can not be read from MEG gradiometer file');
        else
          shape.fid.pnt   = elec.pnt;
          shape.fid.label = elec.label;
          success = 1;
        end
      end
    end

    if ~success
      % try reading it as volume conductor
      % and treat the skin surface as headshape
      try
        vol = read_vol(filename);
        if ~voltype(vol, 'bem')
          error('skin surface can only be extracted from boundary element model');
        else
          if ~isfield(vol, 'skin')
            vol.skin = find_outermost_boundary(vol.bnd);
          end
          shape.pnt = vol.bnd(vol.skin).pnt;
          shape.tri = vol.bnd(vol.skin).tri; % also return the triangulation
          success = 1;
        end
      end
    end

    if ~success
      error('unknown fileformat for head shape information');
    end
end

% this will add the units to the head shape
shape = convert_units(shape);
