function value = loadvar(filename, varname)

% LOADVAR is a helper function for cfg.inputfile

% Copyright (C) 2010, Robert Oostenveld
%
% $Id: loadvar.m 1358 2010-07-06 08:34:26Z roboos $

if nargin<2
  fprintf('reading variable from file ''%s''\n', filename);
else
  fprintf('reading ''%s'' from file ''%s''\n', varname, filename);
end

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

