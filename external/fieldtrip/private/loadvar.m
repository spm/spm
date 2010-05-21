function value = loadvar(filename, varname)

% LOADVAR is a helper function for cfg.inputfile

% Copyright (C) 2010, Robert Oostenveld
%
% $Id: loadvar.m 1083 2010-05-18 09:27:10Z roboos $

fprintf('reading ''%s'' from file ''%s''\n', varname, filename);

var = whos('-file', filename);

if length(var)==1
  filecontent = load(filename); % read the one variable in the file, regardless of how it is called
  value       = filecontent.(var.name);
  clear filecontent
else
  filecontent = load(filename, varname);
  value       = filecontent.(varname);  % read the variable named according to the input specification
  clear filecontent
end

