function value = loadvar(filename, varname)

% LOADVAR is a helper function for cfg.inputfile

% Copyright (C) 2010, Robert Oostenveld
%
% $Id: loadvar.m 7123 2012-12-06 21:21:38Z roboos $

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

