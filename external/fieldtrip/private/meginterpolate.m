function [data] = meginterpolate(cfg, data)

% MEGINTERPOLATE is deprecated, please use MEGREALIGN, MEGPLANAR and/or MEGREPAIR

% Copyright (C) 2003-2006, Robert Oostenveld, F.C. Donders Centre
%
% $Log: meginterpolate.m,v $
% Revision 1.21  2006/10/03 12:59:39  roboos
% updated documentation, explicitely state that it is deprecated, give warnings
%

% set the default configuration 
if ~isfield(cfg, 'realign'),       cfg.realign = 'no';            end
if ~isfield(cfg, 'repair'),        cfg.repair = 'no';             end
if ~isfield(cfg, 'planar'),        cfg.planar = 'no';             end

if strcmp(cfg.repair, 'yes')
  warning('MEGINTERPOLATE is deprecated, please use MEGREPAIR')
  data = megrepair(cfg, data);
else
  fprintf('not repairing bad channels\n');
end

if strcmp(cfg.realign, 'yes')
  warning('MEGINTERPOLATE is deprecated, please use MEGREALIGN')
  data = megrealign(cfg, data);
else
  fprintf('not realigning to template\n');
end

if strcmp(cfg.planar, 'yes')
  warning('MEGINTERPOLATE is deprecated, please use MEGPLANAR')
  data = megplanar(cfg, data);
else
  fprintf('not computing planar gradient\n');
end
