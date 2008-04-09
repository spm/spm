function check_cfg_unused(cfg, a)


% CHECK_CFG_UNUSED will show a warning if the configuration structure
% contains fields that are not used
%
% Use as
%   check_cfg_unused(cfg, fieldnames)
% where fieldnames should be a cell-array with the fields that are
% used.

% Copyright (C) 2004, Robert Oostenveld
% 
% $Log: check_cfg_unused.m,v $
% Revision 1.1  2004/11/19 12:14:34  roboos
% new implementation
% 

% Check the configuration structure for unused fields

msgid = 'FIELDTRIP:configurationUnused';

if ~iscell(a)
  a = {a};
end

b = fieldnames(cfg);
[c, ia, ib] = setxor(a, b);

% ensure that they are row vectors
ia = ia(:)';
ib = ib(:)';

if str2num(version('-release'))>12.1
  % for Matlab versions 6.5 and later
  if length(ib)>0
    warning off backtrace
    warning off verbose
    for i=ib(1:(end-1))
      warning(msgid, 'The field cfg.%s is unused\n', b{i});
    end
    warning on backtrace
    warning on verbose
    warning(msgid, 'The field cfg.%s is unused\n', b{end});
  end
else
  % backward compatible for Matlab 6.1
  for i=ib(1:end)
    warning(sprintf('The field cfg.%s is unused\n', b{i}));
  end
end

