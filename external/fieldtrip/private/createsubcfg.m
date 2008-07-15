function [topcfg] = createsubcfg(topcfg, subname);

% CREATESUBCFG collects the optional arguments for some of the low-level
% functions and puts them in a seperate substructure This function is to
% ensure backward compatibility of end-user scripts, fieldtrip functions
% and documentation that do not use the nested detailled configuration
% but that use a flat configuration.
%
% Use as
%   cfg = createsubcfg(cfg, 'something')
% which will put all fields related to 'something' into cfg.something
%
% See also KEYVAL2CFG, CFG2KEYVAL

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: createsubcfg.m,v $
% Revision 1.9  2008/07/15 19:34:13  roboos
% removed inwardshift from cfg.grid, will have to be reconsidered later
%
% Revision 1.8  2008/07/15 18:18:45  roboos
% added subfield grid
%
% Revision 1.7  2008/07/11 13:18:40  roboos
% removed all lnfilter references, added error to preprocessing and preproc
%
% Revision 1.6  2006/10/12 09:11:43  roboos
% split lcmv and dics, changed keepcsd into keepcov for lcmv
%
% Revision 1.5  2006/10/05 09:58:11  roboos
% removed powmethod for pcc
%
% Revision 1.4  2006/10/02 07:07:04  jansch
% added realfilter-option in beamformer-cfg
%
% Revision 1.3  2006/08/31 08:47:18  roboos
% added filter direction for preproc
%
% Revision 1.2  2006/06/22 12:29:42  roboos
% changed numcomponents into numcomponent for music
%
% Revision 1.1  2006/06/14 11:49:08  roboos
% new implementation, replaces construct_optarg in sourceanalysis and is now also used for cfg.preproc
%

if isfield(topcfg, subname)
  % get the options that are already specified in the substructure
  subcfg = getfield(topcfg, subname);
else
  % start with an empty substructure
  subcfg = [];
end

% add all other relevant options to the substructure
switch subname
  case 'preproc'
    fieldname = {
      'reref'
      'refchannel'
      'implicitref'
      'detrend'
      'blc'
      'blcwindow'
      'bpfiltdir'
      'bpfilter'
      'bpfiltord'
      'bpfilttype'
      'bpfreq'
      'bsfiltdir'
      'bsfilter'
      'bsfiltord'
      'bsfilttype'
      'bsfreq'
      'dftfilter'
      'dftfreq'
      'hpfiltdir'
      'hpfilter'
      'hpfiltord'
      'hpfilttype'
      'hpfreq'
      'lpfiltdir'
      'lpfilter'
      'lpfiltord'
      'lpfilttype'
      'lpfreq'
      'medianfilter'
      'medianfiltord'
      'hilbert'
      'derivative'
      'rectify'
      'boxcar'
      'absdiff'
      };

  case 'grid'
    fieldname = {
      'xgrid'
      'ygrid'
      'zgrid'
      'resolution'
      'filter'
      'leadfield'
      'inside'
      'outside'
      'pos'
      'dim'
      'tight'
      };

  case 'dics'
    fieldname = {
      'feedback'
      'keepfilter'
      'keepmom'
      'lambda'
      'normalize'
      'powmethod'
      'projectnoise'
      'reducerank'
      'keepcsd'
      'realfilter'
      };

  case 'lcmv'
    fieldname = {
      'feedback'
      'keepfilter'
      'keepmom'
      'lambda'
      'normalize'
      'powmethod'
      'projectnoise'
      'reducerank'
      'keepcov'
      };

  case 'pcc'
    fieldname = {
      'feedback'
      'keepfilter'
      'keepmom'
      'lambda'
      'normalize'
%     'powmethod'
      'projectnoise'
      'reducerank'
      'keepcsd'
      'realfilter'
      };

  case {'mne', 'loreta', 'rv'}
    fieldname = {
      'feedback'
      };

  case 'music'
    fieldname = {
      'feedback'
      'numcomponent'
      };

  otherwise
    error('unexpected name of the subfunction');
    fieldname = {};

end % switch subname

for i=1:length(fieldname)
  if ~isfield(subcfg, fieldname{i}) && isfield(topcfg, fieldname{i})
    subcfg = setfield(subcfg, fieldname{i}, getfield(topcfg, fieldname{i}));   % set it in the subconfiguration
    topcfg = rmfield(topcfg, fieldname{i});                                    % remove it in the main configuration
  end
end

% copy the substructure back into the main configuration structure
topcfg = setfield(topcfg, subname, subcfg);

