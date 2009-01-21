function [varargout] = funname(varargin)

% READ_EEP_AVR reads averaged EEG data from an EEProbe *.avr file
% and returns a structure containing the header and data information.
%
% Use as
%   eeg = read_eep_avr(filename)
%
% See also READ_EEP_CNT, READ_EEP_TRG, READ_EEP_REJ

% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  mex(mexsrc);
  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end

