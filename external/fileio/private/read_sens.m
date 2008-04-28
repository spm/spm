function [sens] = read_sens(filename, varargin)

% READ_SENS read sensor positions from various manufacturer specific files. 
% Currently supported are ASA, BESA, Polhemus and Matlab for EEG 
% electrodes and CTF and Neuromag for MEG gradiometers.
%
% Use as
%   grad = read_sens(filename, ...)  % for gradiometers
%   elec = read_sens(filename, ...)  % for electrodes
%
% Additional options should be specified in key-value pairs and can be
%   'fileformat'   string
%
% An electrode definition contain the following fields
%   elec.pnt     Nx3 matrix with carthesian (x,y,z) coordinates of each electrodes
%   elec.label   cell-array of length N with the label of each electrode
%
% A gradiometer definition generally consists of multiple coils per
% channel, e.g.two coils for a 1st order gradiometer in which the
% orientation of the coils is opposite. Each coil is described
% separately and a large "tra" matrix (can be sparse) has to be
% given that defines how the forward computed field is combined over
% the coils to generate the output of each channel. The gradiometer
% definition constsis of the following fields
%   grad.pnt     Mx3 matrix with the position of each coil
%   grad.ori     Mx3 matrix with the orientation of each coil
%   grad.tra     NxM matrix with the weight of each coil into each channel
%   grad.label   cell-array of length N with the label of each of the channels
%
% See also TRANSFORM_SENS, PREPARE_VOL_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2005-2008, Robert Oostenveld
%
% $Log: read_sens.m,v $
% Revision 1.8  2008/04/14 20:51:36  roboos
% added convert_units
%
% Revision 1.7  2008/04/11 16:17:22  roboos
% added polhemus_fil
%
% Revision 1.6  2008/03/20 13:43:14  roboos
% added support for besa_pos
%
% Revision 1.5  2008/03/18 12:34:30  roboos
% fixed bug: added varargin to input arguments, thanks to Juan
%
% Revision 1.4  2008/03/06 09:27:54  roboos
% updated documentation
%
% Revision 1.3  2008/03/05 11:06:11  roboos
% test the presence of the fileio toolbox, needed when this function is included in forwinv
%
% Revision 1.2  2008/03/05 10:54:05  roboos
% added optional argument for fileformat
% some documentation changes
%
% Revision 1.1  2008/01/28 20:10:11  roboos
% new functions based on existing fieldtrip code
%

% test whether the file exists
if ~exist(filename)
  error(sprintf('file ''%s'' does not exist', filename));
end

% get the options
fileformat = keyval('fileformat',  varargin);

% determine the filetype
if isempty(fileformat)
  fileformat = filetype(filename);
end

switch fileformat

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % read the content from various files that contain EEG electrode positions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'asa_elc'
    sens = read_asa_elc(filename);
  
  case 'polhemus_pos'
    sens = read_brainvision_pos(filename);

  case 'besa_pos'
    tmp = importdata(filename);
    if ~isnumeric(tmp)
      error('unexpected file format for fileformat=besa_pos')
    end
    [nchan,nrow] = size(tmp);
    if nrow==3
      sens.pnt = tmp;
    elseif nrow==9
      pnt1 = tmp(:,1:3);  % bottom coil
      pnt2 = tmp(:,4:6);  % top coil
      ori  = tmp(:,7:9);  % orientation of bottom coil
      sens.pnt = [pnt1; pnt2];
      sens.ori = [ori; ori];
      sens.tra = [eye(nchan) -eye(nchan)];
    else
      error('unexpected file format for fileformat=besa_pos')
    end
    [p, f, x] = fileparts(filename);
    elpfile = fullfile(p, [f '.elp']);
    elafile = fullfile(p, [f '.ela']);
    if exist(elpfile, 'file')
      warning(sprintf('reading channel labels from %s', elpfile));
      % read the channel names from the accompanying ELP file
      lbl = importdata(elpfile);
      sens.label = strrep(lbl.textdata(:,2) ,'''', '');
    elseif exist(elafile, 'file')
      warning(sprintf('reading channel labels from %s', elafile));
      % read the channel names from the accompanying ELA file
      lbl = importdata(elafile);
      lbl = strrep(lbl, 'MEG ', ''); % remove the channel type
      lbl = strrep(lbl, 'EEG ', ''); % remove the channel type
      sens.label = lbl;
    else
      % the file does not have channel labels in it
      warning('creating fake channel names for besa_pos');
      for i=1:nchan
        sens.label{i} = sprintf('%03d', i);
      end
    end

  case 'besa_sfp'
    tmp        = importdata(filename);
    sens.label = tmp.textdata;
    sens.pnt   = tmp.data;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % gradiometer information is always stored in the header of the MEG dataset
  % hence uses the standard fieldtrip/fileio read_header function
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case {'ctf_ds', 'ctf_res4', 'neuromag_fif', 'yokogawa_ave', 'yokogawa_con', 'yokogawa_raw'}
    % check the availability of the required low-level toolbox
    % this is required because the read_sens function is also on itself included in the forwinv toolbox
    hastoolbox('fileio');
    hdr = read_header(filename);
    sens = hdr.grad;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % these are created at the FIL in London with a polhemus tracker
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case 'polhemus_fil'
    [sens.fid, sens.pnt] = read_polhemus_fil(filename, 0, 0);

    % the file does not have channel labels in it
    warning('creating fake channel names for polhemus_fil');
    for i=1:nchan
      sens.label{i} = sprintf('%03d', i);
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % matlab files can contain either electrodes or gradiometers
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  case 'matlab'
    matfile = filename;   % this solves a problem with the matlab compiler v3
    warning('off', 'MATLAB:load:variableNotFound');
    tmp = load(matfile, 'elec', 'grad', 'sens', 'elc');
    warning('on', 'MATLAB:load:variableNotFound');
    if isfield(tmp, 'grad')
      sens = getfield(tmp, 'grad');
    elseif isfield(tmp, 'elec')
      sens = getfield(tmp, 'elec');
    elseif isfield(tmp, 'sens')
      sens = getfield(tmp, 'sens');
    elseif isfield(tmp, 'elc')
      sens = getfield(tmp, 'elc');
    else
      error('no electrodes or gradiometers found in Matlab file');
    end

  otherwise
    error('unknown fileformat for electrodes or gradiometers');
end

if senstype(sens, 'eeg')
  % only keep positions and labels in case of EEG electrodes
  dum  = sens;
  sens = [];
  sens.pnt   = dum.pnt;
  sens.label = dum.label;
end

% this will add the units to the sensor array
sens = convert_units(sens);
