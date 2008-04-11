function [shape] = read_headshape(filename, varargin)

% READ_HEADSHAPE reads the fiducials and/or the measured headshape
% from a variety of files (like CTF and Polhemus). The headshape  and
% fiducials can for example be used for coregistration.
%
% Use as
%   [shape] = read_headshape(filename)
%
% See also READ_VOL, READ_SENS

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: read_headshape.m,v $
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

switch fileformat
  case 'ctf_ds'
    [p, f, x] = fileparts(filename);
    filename = fullfile(p, f, [f 'hs']);
    shape = read_ctf_hs(filename);

  case 'ctf_hc'
    shape = read_ctf_hc(filename);

  case 'ctf_shape'
    shape = read_ctf_shape(filename);

  case '4d_hs'
    shape.pnt = read_bti_hs(filename);

  case 'polhemus_fil'
    [fid, sens] = read_polhemus_fil(filename, 0, 0);
    shape.fid  = fid;
    shape.sens = sens;

  case 'matlab'
    tmp = load(filename);
    if isfield(tmp, 'shape')
      shape = tmp.shape;
    else
      error('no headshape found in Matlab file');
    end

  otherwise
    % try reading it as electrode positions
    try
      elec = read_sens(filename);
      if ~senstype(elec, 'eeg')
        error('headshape information can not be read from MEG gradiometer file');
      else
        shape.pnt = elec.pnt;
      end
    catch
      disp(lasterr);
      error('unknown fileformat for head shape information');
    end
end

