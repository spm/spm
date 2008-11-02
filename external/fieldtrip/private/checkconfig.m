function [cfg] = checkconfig(cfg, varargin)

% CHECKCONFIG checks the input cfg of the main FieldTrip functions.
%
% 1: It checks whether the cfg contains all the required options, it gives
% a warning when renamed or deprecated options are used, and it makes sure
% no forbidden options are used. If necessary and possible, this function
% will adjust the cfg to the input requirements. If the input cfg does NOT
% correspond to the requirements, this function gives an elaborate warning
% message.
%
% 2: It controls the relevant cfg options that are being passed on to other
% functions, by putting them into substructures or converting them into the
% required format.
%
% 3: It controls the output cfg (data.cfg) such that it only contains
% relevant and used fields. This part of the functionality is still under
% construction!
%
%
% Use as
%   [cfg] = checkconfig(cfg, ...)
%
% Optional input arguments should be specified as key-value pairs and can include
%   renamed           = {'old',  'new'}        % list the old and new option
%   renamedval        = {'opt',  'old', 'new'} % list option and old and new value
%   required          = {'opt1', 'opt2', etc.} % list the required options
%   deprecated        = {'opt1', 'opt2', etc.} % list the deprecated options
%   unused            = {'opt1', 'opt2', etc.} % list the unused options, these will be removed and a warning is issued
%   forbidden         = {'opt1', 'opt2', etc.} % list the forbidden options, these result in an error
%   createsubcfg      = {'subname', etc.}      % list the names of the subcfg
%   dataset2files     = 'yes', 'no'            % converts dataset into headerfile and datafile
%   trackconfig       = 'on', 'off'            % turn config tracking on/off
%
% See also CHECKDATA

% Copyright (C) 2007-2008, Robert Oostenveld, Saskia Haegens
%
% $Log: checkconfig.m,v $
% Revision 1.14  2008/11/02 10:56:34  roboos
% explicit sharing of code for dataset2files with fileio read_header/data
%
% Revision 1.13  2008/10/28 19:20:18  sashae
% configtracking can be turned on/off with trackconfig.
% user can request report on unused options with cfg.trackconfig='report', or get cleaned cfg with cfg.trackconfig='cleanup'.
%
% Revision 1.12  2008/10/13 13:38:59  sashae
% change in dataset2files code: empty dataset/headerfile/datafile fields are removed
%
% Revision 1.11  2008/10/13 12:41:33  jansch
% added projectmom for lcmv (in contrast to pcc the projection is done within
% beamformer_lcmv, instead of in sourcedescriptives). Probably this should
% change back when a clean version of sourcedescriptives is implemented.
%
% Revision 1.10  2008/10/10 12:33:04  sashae
% incorporated dataset2files
%
% Revision 1.9  2008/10/10 10:50:34  sashae
% updated documentation
%
% Revision 1.8  2008/10/02 14:06:21  roboos
% get the fields from ft_default and add them to the cfg structure
% implemented cfg.checkconfig=silent/loose/pedantic (default is in ft_default, i.e. fieldtripdefs function)
%
% Revision 1.7  2008/10/02 12:35:14  roboos
% added option "unused", renamed tracking to configtracking
%
% Revision 1.6  2008/10/01 15:45:29  sashae
% incorporated createsubcfg
%
% Revision 1.5  2008/09/30 13:05:19  roboos
% adedd first version of configuration tracking for testing
%
% Revision 1.4  2008/09/23 12:05:33  sashae
% some small changes; checkconfig can now handle empty cfgs
%
% Revision 1.3  2008/09/18 10:01:57  sashae
% added 'renamedval' which checks/adjusts renamed values
%
% Revision 1.2  2008/09/18 08:33:48  sashae
% new version: checks required, renamed, deprecated and forbidden configuration options,
% adjusts where possible and gives warning/error messages. to be implemented in all main
% fieldtrip functions, comparable to checkdata
%
% Revision 1.1  2008/07/08 15:39:22  roboos
% initial version for Saskia to work on
%

if isempty(cfg)
  cfg = struct; % ensure that it is an empty struct, not empty double
end

global ft_default
if isempty(ft_default)
  ft_default = struct;
end
fieldsused = fieldnames(ft_default);
for i=1:length(fieldsused)
  fn = fieldsused{i};
  if ~isfield(cfg, fn),
    cfg.(fn) = ft_default.(fn);
  end
end

renamed         = keyval('renamed',         varargin);
renamedval      = keyval('renamedval',      varargin);
required        = keyval('required',        varargin);
deprecated      = keyval('deprecated',      varargin);
unused          = keyval('unused',          varargin);
forbidden       = keyval('forbidden',       varargin);
createsubcfg    = keyval('createsubcfg',    varargin);
dataset2files   = keyval('dataset2files',   varargin);
trackconfig     = keyval('trackconfig',     varargin);
if isempty(trackconfig)
  if strcmp(cfg.trackconfig, 'report') || strcmp(cfg.trackconfig, 'cleanup')
    trackconfig = 'on'; % turn on configtracking if user requests report/cleanup
  else
    trackconfig = 'no';
  end
end

% these should be cell arrays and not strings
if ischar(required),   required   = {required};   end
if ischar(deprecated), deprecated = {deprecated}; end
if ischar(unused),     unused     = {unused};     end
if ischar(forbidden),  forbidden  = {forbidden};  end

if isfield(cfg, 'checkconfig')
  silent   = strcmp(cfg.checkconfig, 'silent');
  loose    = strcmp(cfg.checkconfig, 'loose');
  pedantic = strcmp(cfg.checkconfig, 'pedantic');
else
  silent   = false;
  loose    = true;
  pedantic = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename old to new options, give warning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(renamed)
  fieldsused = fieldnames(cfg);
  if any(strcmp(renamed{1}, fieldsused))
    cfg = setfield(cfg, renamed{2}, (getfield(cfg, renamed{1})));
    cfg = rmfield(cfg, renamed{1});
    warning(sprintf('use cfg.%s instead of cfg.%s', renamed{2}, renamed{1}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename old to new value, give warning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(renamedval) && isfield(cfg, renamedval{1})
  if strcmpi(getfield(cfg, renamedval{1}), renamedval{2})
    cfg = setfield(cfg, renamedval{1}, renamedval{3});
    warning(sprintf('use cfg.%s = %s instead of cfg.%s = %s', renamedval{1}, renamedval{3}, renamedval{1}, renamedval{2}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for required fields, give error when missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(required)
  fieldsused = fieldnames(cfg);
  [c, ia, ib] = setxor(required, fieldsused);
  if ~isempty(ia)
    error(sprintf('The field cfg.%s is required\n', required{ia}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for deprecated fields, give warning when present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(deprecated)
  fieldsused = fieldnames(cfg);
  if any(ismember(deprecated, fieldsused))
    if silent
      % don't mention it
    elseif loose
      warning(sprintf('The option cfg.%s is deprecated, support is no longer guaranteed\n', deprecated{ismember(deprecated, fieldsused)}));
    elseif pedantic
      error(sprintf('The option cfg.%s is deprecated, support is no longer guaranteed\n', deprecated{ismember(deprecated, fieldsused)}));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for unused fields, give warning when present and remove them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(unused)
  fieldsused = fieldnames(cfg);
  if any(ismember(unused, fieldsused))
    if silent
      % don't mention it
      cfg = rmfield(cfg, unused(ismember(unused, fieldsused)));
    elseif loose
      warning(sprintf('The field cfg.%s is unused, it will be removed from your configuration\n', unused{ismember(unused, fieldsused)}));
      cfg = rmfield(cfg, unused(ismember(unused, fieldsused)));
    elseif pedantic
      error(sprintf('The field cfg.%s is unused\n', unused{ismember(unused, fieldsused)}));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for forbidden fields, give error when present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(forbidden)
  fieldsused = fieldnames(cfg);
  if any(ismember(forbidden, fieldsused))
    if silent
      % don't mention it
      cfg = rmfield(cfg, forbidden(ismember(forbidden, fieldsused)));
    elseif loose
      warning(sprintf('The field cfg.%s is forbidden, it will be removed from your configuration\n', forbidden{ismember(forbidden, fieldsused)}));
      cfg = rmfield(cfg, forbidden(ismember(forbidden, fieldsused)));
    elseif pedantic
      error(sprintf('The field cfg.%s is forbidden\n', forbidden{ismember(forbidden, fieldsused)}));
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createsubcfg
%
% This collects the optional arguments for some of the low-level
% functions and puts them in a separate substructure. This function is to
% ensure backward compatibility of end-user scripts, fieldtrip functions
% and documentation that do not use the nested detailled configuration
% but that use a flat configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(createsubcfg)
  for j=1:length(createsubcfg)
    subname = createsubcfg{j};

    if isfield(cfg, subname)
      % get the options that are already specified in the substructure
      subcfg = getfield(cfg, subname);
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
          'projectmom'
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
          %'powmethod'
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
      if ~isfield(subcfg, fieldname{i}) && isfield(cfg, fieldname{i})
        subcfg = setfield(subcfg, fieldname{i}, getfield(cfg, fieldname{i}));  % set it in the subconfiguration
        cfg = rmfield(cfg, fieldname{i});                                      % remove it from the main configuration
      end
    end

    % copy the substructure back into the main configuration structure
    cfg = setfield(cfg, subname, subcfg);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataset2files
%
% Converts cfg.dataset into cfg.headerfile and cfg.datafile if neccessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(dataset2files) && strcmp(dataset2files, 'yes')

  % start with empty fields if they are not present
  if ~isfield(cfg, 'dataset')
    cfg.dataset = [];
  end
  if ~isfield(cfg, 'datafile')
    cfg.datafile = [];
  end
  if ~isfield(cfg, 'headerfile')
    cfg.headerfile = [];
  end

  if ~isempty(cfg.dataset)
    if strcmp(cfg.dataset, 'gui');
      d = uigetdir;
      if d==0
        [f, p] = uigetfile;
        if f==0
          error('You should select a dataset file or directory');
        else
          d = fullfile(p, f);
        end
      end
      cfg.dataset = d;
    end

    % the following code is shared with fileio read_header/read_data
    % therefore the three local variables are used outside of the cfg
    filename   = cfg.dataset;
    datafile   = [];
    headerfile = [];
    switch filetype(filename)
      case '4d_pdf'
        datafile   = filename;
        headerfile = [datafile '.m4d'];
        sensorfile = [datafile '.xyz'];
      case {'4d_m4d', '4d_xyz'}
        datafile   = filename(1:(end-4)); % remove the extension
        headerfile = [datafile '.m4d'];
        sensorfile = [datafile '.xyz'];
      case '4d'
        [path, file, ext] = fileparts(filename);
        datafile   = fullfile(path, [file,ext]);
        headerfile = fullfile(path, [file,ext]);
        configfile = fullfile(path, 'config');
      case {'ctf_ds', 'ctf_old'}
        % convert CTF filename into filenames
        [path, file, ext] = fileparts(filename);
        if isempty(path) && isempty(file)
          % this means that the dataset was specified as the present working directory, i.e. only with '.'
          filename = pwd;
          [path, file, ext] = fileparts(filename);
        end
        headerfile = fullfile(filename, [file '.res4']);
        datafile   = fullfile(filename, [file '.meg4']);
        if length(path)>3 && strcmp(path(end-2:end), '.ds')
          filename = path; % this is the *.ds directory
        end
      case 'ctf_meg4'
        [path, file, ext] = fileparts(filename);
        if isempty(path)
          path = pwd;
        end
        headerfile = fullfile(path, [file '.res4']);
        datafile   = fullfile(path, [file '.meg4']);
        if length(path)>3 && strcmp(path(end-2:end), '.ds')
          filename = path; % this is the *.ds directory
        end
      case 'ctf_res4'
        [path, file, ext] = fileparts(filename);
        if isempty(path)
          path = pwd;
        end
        headerfile = fullfile(path, [file '.res4']);
        datafile   = fullfile(path, [file '.meg4']);
        if length(path)>3 && strcmp(path(end-2:end), '.ds')
          filename = path; % this is the *.ds directory
        end
      case 'brainvision_vhdr'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        if exist(fullfile(path, [file '.eeg']))
          datafile   = fullfile(path, [file '.eeg']);
        elseif exist(fullfile(path, [file '.seg']))
          datafile   = fullfile(path, [file '.seg']);
        elseif exist(fullfile(path, [file '.dat']))
          datafile   = fullfile(path, [file '.dat']);
        end
      case 'brainvision_eeg'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        datafile   = fullfile(path, [file '.eeg']);
      case 'brainvision_seg'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        datafile   = fullfile(path, [file '.seg']);
      case 'brainvision_dat'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.vhdr']);
        datafile   = fullfile(path, [file '.dat']);
      case 'fcdc_matbin'
        [path, file, ext] = fileparts(filename);
        headerfile = fullfile(path, [file '.mat']);
        datafile   = fullfile(path, [file '.bin']);
      otherwise
        % convert filename into filenames, assume that the header and data are the same
        datafile   = filename;
        headerfile = filename;
    end
    % end sharing with fileio read_header/read_data
    % put everything back into the cfg
    cfg.dataset    = filename;
    cfg.datafile   = datafile;
    cfg.headerfile = headerfile;

  elseif ~isempty(cfg.datafile) && isempty(cfg.headerfile);
    % assume that the datafile also contains the header
    cfg.headerfile = cfg.datafile;
  elseif isempty(cfg.datafile) && ~isempty(cfg.headerfile);
    % assume that the headerfile also contains the data
    cfg.datafile = cfg.headerfile;
  end
  % remove empty fields (otherwise a subsequent check on required fields doesn't make any sense)
  if isempty(cfg.dataset),    cfg=rmfield(cfg, 'dataset');    end
  if isempty(cfg.headerfile), cfg=rmfield(cfg, 'headerfile'); end
  if isempty(cfg.datafile),   cfg=rmfield(cfg, 'datafile');   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configtracking
%
% switch configuration tracking on/off
% FIXME dit is slechts een eerste versie "voor de leuk"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  switch trackconfig
    case 'on' % turn on configtracking
      if isa(cfg, 'struct')
        cfg = configuration(cfg);
      elseif isa(cfg, 'configuration')
        % configtracking is already turned on, do nothing
      else
        warning('cannot convert "%s" to a configuration object for tracking', class(cfg));
      end

    case 'off' % turn off configtracking
      if strcmp(cfg.trackconfig, 'report') % give feedback if requested
        access(cfg);
      end

      if strcmp(cfg.trackconfig, 'cleanup') % clean up cfg if requested
        cfg=cleancfg(cfg);
      end

      cfg = struct(cfg); % convert back to struct

    otherwise
      % do nothing
  end
catch
  disp(lasterr);
end
