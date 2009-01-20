function [varargout] = reset(varargin)

% RESET Set the value of the input to zero, using pass-by-reference.
%
% Example
%  a = 1;
%  reset(a);
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

