function [vol, sens] = prepare_vol_sens(vol, sens, varargin)

% PREPARE_VOL_SENS does some bookkeeping to ensures that the volume
% conductor model and the sensor array are appropriate. Furthermore it
% takes care of pre-computations that can be done efficiently prior to the
% leadfield calculations.
%
% The prepare_vol_sens has different roles for EEG and for MEG and % for
% the different volume conductor models. it for example projects the 3D
% electrode positions onto the skin compartment of the volume conductor to
% ensure that they do not float above the surface with a few mm due to
% mis-alignment. Or for example for MEG gradiometer sensors with a
% multisphere volume conductor it ensured that each coil of the gradiometer
% array is associated with a sphere.
%
% Use as
%   [vol, sens] = prepare_vol_sens(vol, sens, ...)
%
% The vol structure represents a volume conductor model, its contents
% depend on the type of model. The sens structure represents a sensor
% arary, i.e. EEG electrodes or MEG gradiometers.
%
% Additional options should be specified in key-value pairs.
%
% See also READ_VOL, READ_SENS, TRANSFORM_VOL, TRANSFORM_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2004-2008, Robert Oostenveld
%
% $Log: prepare_vol_sens.m,v $
% Revision 1.1  2008/03/06 09:30:36  roboos
% Created skeleton implementation according to how it should be for the forwinv toolbox, i.e. fieldtrip independent, so that it can be included in spm8.
% The functionality should be moved from the existing fieldtrip/private/prepare_vol_sens.m function into this new function.
%

% get the options
% fileformat = keyval('fileformat',  varargin);

% determine whether the input contains EEG or MEG seosors
iseeg = senstype(sens, 'eeg');
ismeg = senstype(sens, 'meg');

if ismeg && iseeg
  % this is something that could be implemented relatively easily
  error('simultaneous EEG and MEG not supported');

elseif ~ismeg && ~iseeg
  error('the input does not look like EEG, nor like MEG');

elseif ismeg
  switch voltype(vol)
    case 'multisphere'
    case 'neuromag'
    case 'nolte'
    case 'infinite'
    case 'multisphere'
    case {'singlesphere', 'concentric'}
    otherwise
      error('unsupported volume conductor model for MEG');
  end

elseif iseeg
  switch voltype(vol)
    case 'singlesphere'
    case 'infinite'
    case {'singlesphere', 'concentric'}
    case 'bem'
    otherwise
      error('unsupported volume conductor model for EEG');
  end

end % if iseeg or ismeg

