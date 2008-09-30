function [cfg] = checkconfig(cfg, varargin)

% CHECKCONFIG checks the input cfg of the main FieldTrip functions. It
% checks whether the cfg contains all the required options, it gives a
% warning when renamed or deprecated options are used, and it makes sure
% no forbidden options are used. If necessary and possible, this function
% will adjust the cfg to the input requirements.
%
% If the input cfg does NOT correspond to the requirements, this function
% is supposed to give an elaborate warning message and if applicable point
% the user to external documentation.
%
% Use as
%   [cfg] = checkconfig(cfg, ...)
%
% Optional input arguments should be specified as key-value pairs and can include
%   renamed           = {'old',  'new'}        % list the old and new option
%   renamedval        = {'opt',  'old', 'new'} % list option and old and new value
%   required          = {'opt1', 'opt2', etc.} % list the required options
%   deprecated        = {'opt1', 'opt2', etc.} % list the deprecated options
%   forbidden         = {'opt1', 'opt2', etc.} % list the forbidden options
%   tracking          = 'on', 'off' or 'report'
%
% See also CHECKDATA

% Copyright (C) 2007-2008, Robert Oostenveld
%
% $Log: checkconfig.m,v $
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

renamed         = keyval('renamed',     varargin);
renamedval      = keyval('renamedval',  varargin);
required        = keyval('required',    varargin);
deprecated      = keyval('deprecated',  varargin);
forbidden       = keyval('forbidden',   varargin);
tracking        = keyval('tracking',    varargin); if isempty(tracking), tracking = 'no'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename old to new options, give warning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(renamed)
  try fieldsused = fieldnames(cfg); catch fieldsused={}; end
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
  try fieldsused = fieldnames(cfg); catch fieldsused={}; end
  [c, ia, ib] = setxor(required, fieldsused);
  if ~isempty(ia)
    error(sprintf('The field cfg.%s is required\n', required{ia}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for deprecated fields, give warning when present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(deprecated)
  try fieldsused = fieldnames(cfg); catch fieldsused={}; end
  if any(ismember(deprecated, fieldsused))
    warning(sprintf('the option cfg.%s is deprecated, support is no longer guaranteed\n', deprecated{[ismember(deprecated, fieldsused)]}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for forbidden fields, give error when present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(forbidden)
  try fieldsused = fieldnames(cfg); catch fieldsused={}; end
  if any(ismember(forbidden, fieldsused))
    error(sprintf('The field cfg.%s is forbidden\n', forbidden{[ismember(forbidden, fieldsused)]}));
  end
end

try
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % switch configuration tracking on/off
  % FIXME dit is slechts een eerste versie "voor de leuk"
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch tracking
    case 'on'
      if isa(cfg, 'struct')
        cfg = configuration(cfg);
      else
        warning('cannot convert "%s" to a configuration object for tracking', class(cfg));
      end

    case 'report'
      access(cfg);

    case 'off'
      cfg = struct(cfg);

    otherwise
      % do nothing
  end
catch
  disp(lasterr);
end
