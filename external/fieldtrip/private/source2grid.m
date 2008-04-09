function [grid] = source2grid(source)

% SOURCE2GRID removes the fields from a source structure that are
% not neccessary to reuse the dipole grid in another source estimation.
%
% Use as
%   [grid] = source2grid(source);
%
% The resulting grid can be used in the configuration of another
% run of SOURCANALYSIS.
%
% See also SOURCE2SPARSE, SOURCE2FULL

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: source2grid.m,v $
% Revision 1.3  2006/05/17 08:41:18  roboos
% also keep the filters if present
%
% Revision 1.2  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%
% Revision 1.1  2004/08/03 09:06:19  roboos
% initial implementation of these helper functions for beamformer sourceanalysis
%

grid.xgrid   = source.xgrid;
grid.ygrid   = source.ygrid;
grid.zgrid   = source.zgrid;
grid.pos     = source.pos;
grid.dim     = source.dim;
grid.inside  = source.inside;
grid.outside = source.outside;
if issubfield(source, 'filter')
  grid.filter = source.filter;
elseif issubfield(source, 'avg.filter')
  grid.filter = source.avg.filter;
end
if isfield(source, 'leadfield')
  grid.leadfield = source.leadfield;
end
