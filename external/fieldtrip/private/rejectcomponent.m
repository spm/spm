function [data] = rejectcomponent(cfg, comp, data)

% REJECTCOMPONENT reconstructs an ICA decomposition after removing 
% the independent components that contain the artifacts.  This function 
% does not automatically detect the artifact components, you will 
% have to do that yourself.
%
% Use as
%    [data] = rejectcomponent(cfg, comp)
% or
%    [data] = rejectcomponent(cfg, comp, data)
% where comp is the result of COMPONENTANALYSIS and data is the output
% of PREFPROCESSING.
% 
% The configuration should contain
%   cfg.component = list of components to remove, e.g. [1 4 7]
% 
% If you have performed a decomposition in which all components were
% computed, the data can be reconstructed from the components, leaving
% the artefact components out. In that case you do not have to specify the data itself. 
%
% If you have performed a decomposition of the data resuling in only a few
% components, the data cannot be fully reconstructed. By also giving the
% data as input, the artefact components can be filtered out of the data.
%
% See also COMPONENTANALYSIS, PREFPROCESSING

% Undocumented local options:
% cfg.channel
% cfg.previous
% cfg.version

% Copyright (C) 2005, Robert Oostenveld
% 
% $Log: rejectcomponent.m,v $
% Revision 1.6  2007/03/04 14:45:28  chrhes
% fixed a small bug by adding the third input argument 'data' to the function
% interface (declaration) in the first line of the file. Also added an error
% message as a reminder that component removal from incomplete decompositions
% is not yet implemented: it is not clear how this should work, probably the
% optimal approach is to estimate the component waveforms using MMSE, but the
% question is whether to use a global cross-covariance estimate or a local
% (trial based) estimate to compute the adaptive spatial filter.
%
% Revision 1.5  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.4  2006/01/30 17:00:46  roboos
% added initial support for incomplete decompositions (requires input data)
%
% Revision 1.3  2005/10/27 09:37:55  jansch
% fixed bug in assignment of labels
%
% Revision 1.2  2005/10/06 07:11:09  roboos
% changed from DOS to UNIX ascii format
%
% Revision 1.1  2005/09/29 16:23:24  jansch
% first implementation
%

if ~isfield(cfg, 'component'), cfg.component = []; end

ntrials = length(comp.trial);
ncomps  = length(comp.label);

if min(cfg.component)<1
  error('you cannot remove components that are not present in the data');
end

if max(cfg.component)>ncomps
  error('you cannot remove components that are not present in the data');
end

% determine the number of channels and the number of components
Nchan = length(comp.topolabel);
Ncomp = length(comp.label);

if Nchan==Ncomp
  if nargin>2
    warning('not using input data for reconstruction');
  end
  % set the rejected compponent amplitudes to zero 
  comp.topo(:,cfg.component) = 0;
  fprintf('removing %d components\n', length(cfg.component)); 
  fprintf('keeping %d components\n',  ncomps-length(cfg.component));
  % recontruct the data from the independent components 
  for i=1:ntrials
    data.trial{i} = comp.topo * comp.trial{i}; 
  end
  data.fsample = comp.fsample;
  data.time    = comp.time;  
  data.label   = comp.cfg.channel; 
  try, data.grad    = comp.grad;    end
else
  if nargin<3
    error('incomplete component decomposition, this also requires input data');
  end
  [seldat, selcomp] = match_str(data.label, comp.topolabel);
  comp = comp.topo(selcomp,cfg.component);
  % FIXME, continue implementation here...
  error('FIXME: this section of code is still under construction');
end 

% add the version details of this function call to the configuration 
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath'); 
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: rejectcomponent.m,v 1.6 2007/03/04 14:45:28 chrhes Exp $';
% remember the configuration details of the input data 
try, cfg.previous = comp.cfg; end
% keep the configuration in the output
data.cfg = cfg;

