function [data] = rejectcomponent(cfg, comp, data)

% REJECTCOMPONENT backprojects an ICA decomposition to the channel
% level after removing the independent components that contain the
% artifacts. This function does not automatically detect the artifact
% components, you will have to do that yourself.
%
% Use as
%    [data] = rejectcomponent(cfg, comp)
% where the input comp is the result of COMPONENTANALYSIS. The output
% data will have the same format as the output of PREFPROCESSING.
% 
% The configuration should contain
%   cfg.component = list of components to remove, e.g. [1 4 7]
% 
% See also COMPONENTANALYSIS, PREFPROCESSING

% Copyright (C) 2005-2009, Robert Oostenveld
% 
% $Log: rejectcomponent.m,v $
% Revision 1.9  2009/03/23 20:03:50  roboos
% removed obsolete and non-functional code
%
% Revision 1.8  2008/11/04 20:20:34  roboos
% added optional baseline correction
%
% Revision 1.7  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
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

fieldtripdefs

if ~isfield(cfg, 'component'), cfg.component = [];      end

ntrials = length(comp.trial);
ncomps  = length(comp.label);

if min(cfg.component)<1
  error('you cannot remove components that are not present in the data');
end

if max(cfg.component)>ncomps
  error('you cannot remove components that are not present in the data');
end

% set the rejected compponent amplitudes to zero 
fprintf('removing %d components\n', length(cfg.component)); 
fprintf('keeping %d components\n',  ncomps-length(cfg.component));
comp.topo(:,cfg.component) = 0;

% recontruct the data from the independent components 
data = [];
for i=1:ntrials
  data.trial{i} = comp.topo * comp.trial{i}; 
end
data.fsample   = comp.fsample;
data.time      = comp.time;  
data.label     = comp.topolabel;
try, data.grad = comp.grad; end

% add the version details of this function call to the configuration 
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath'); 
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: rejectcomponent.m,v 1.9 2009/03/23 20:03:50 roboos Exp $';
% remember the configuration details of the input data 
try, cfg.previous = comp.cfg; end
% keep the configuration in the output
data.cfg = cfg;

