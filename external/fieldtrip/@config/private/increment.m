function [varargout] = increment(varargin)

% INCREMENT Increments the input by one, using pass-by-reference.
%
% Example
%  a = 1;
%  increment(a);
%  disp(a);

% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
fname  = mfilename('fullpath');
mexsrc = [fname '.c'];
[mexdir, mexname] = fileparts(fname);

try
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  mex(mexsrc);
  cd(pwdir);
catch
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
end

