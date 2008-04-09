function check_cfg_required(cfg, a)

% CHECK_CFG_REQUIRED will show an error if any of the fieldnames is not
% present in the configuration structure
%
% Use as
%   check_cfg_required(cfg, fieldnames)
% where fieldnames should be a cell-array containing the fields that
% are required.

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: check_cfg_required.m,v $
% Revision 1.1  2004/11/19 12:14:34  roboos
% new implementation
%

msgid = 'FIELDTRIP:configurationRequired';

if ~iscell(a)
  a = {a};
end

b = fieldnames(cfg);
[c, ia, ib] = setxor(a, b);

% ensure that they are row vectors
ia = ia(:)';
ib = ib(:)';

if length(ia)>0
  if str2num(version('-release'))>12.1
    % for Matlab versions 6.5 and later
    error(msgid, 'The field cfg.%s is required\n', a{1});
  else
    % backward compatible for Matlab 6.1
    error(sprintf('The field cfg.%s is required\n', a{1}));
  end
end

