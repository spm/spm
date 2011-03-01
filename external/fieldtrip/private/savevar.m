function savevar(filename, varname, value)

% SAVEVAR is a helper function for cfg.outputfile

% Copyright (C) 2010, Robert Oostenveld
%
% $Id: savevar.m 3014 2011-03-01 15:36:42Z eelspa $

fprintf('writing ''%s'' to file ''%s''\n', varname, filename);

eval(sprintf('%s = value;', varname));

s = whos(varname);

% if variable < ~500 MB, store it in old (uncompressed) format, which is
% faster
if (s.bytes < 500000000)
  save(filename, varname, '-v6');
else
  save(filename, varname, '-v7.3');
end
