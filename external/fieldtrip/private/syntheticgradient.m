function [data] = syntheticgradient(cfg, data);

% SYNTHETICGRADIENT applies CTFs higher-order synthetic gradients 
% to preprocessed or averaged data and to the corresponding 
% gradiometer definition.
%
% Use as
%   [data] = syntheticgradient(cfg, data);
% where data should come from PREPROCESSING or TIMELOCKANALYSIS
% and the configuration should contain
%
% cfg.input  = 'none', 'G1BR', 'G2BR' or 'G3BR' specifies the balancing
%              of the input data (default is 'none')
%
% cfg.output = 'none', 'G1BR', 'G2BR' or 'G3BR' specifies the balancing
%              of the output data (default is 'none')
%
% Currently it is not possible to determine the balancing of the data from
% the preprocessed data or from the CTF *.res4 header. Therefore it is
% important to specify the balancing of the input data correctly.

% Copyright (C) 2004-2006, Robert Oostenveld
%
% $Log: syntheticgradient.m,v $
% Revision 1.8  2006/03/29 08:40:26  roboos
% fixed documentation
%
% Revision 1.7  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.6  2005/05/17 17:50:39  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.5  2004/07/02 13:48:36  roboos
% detected a bug in the grad computation and added an error messgae there
% the bug is still not fixed
%
% Revision 1.4  2004/07/02 12:37:52  roboos
% removed the default for the input and output
% added support for averaged data
%
% Revision 1.3  2004/07/02 11:51:17  roboos
% fixed some small bugs
%
% Revision 1.2  2004/07/02 11:40:26  roboos
% fixed bug in switch-otherwise
%
% Revision 1.1  2004/07/02 11:32:50  roboos
% new implementation, requires external coefficients file for the transformation
%

% these MUST be specified, since they cannot be determined from the data itself
if ~isfield(cfg, 'input')
  error('you must specify the balancing of the input data');
end
if ~isfield(cfg, 'output')
  error('you must specify the balancing of the output data');
end

% read the coefficients from the FCDC specific file
load ctf_coef

if isfield(data, 'avg')
  % treat the average just as if it was a single trial
  data.trial = {};
  data.trial{1} = data.avg;
  isavgdata = 1;
else
  isavgdata = 0;
end

% find the matching channels, channel order corresponds to coef
[selcoef, seldat] = match_str(coef.label, data.label);

if length(seldat)<length(coef.label)
  % in principle this could be less restrictive, but that requires more bookkeeping
  error('not all required channels were found in the data');
end

% find the matching channels, channel order corresponds to coef
[selcoef, selgrad] = match_str(coef.label, data.grad.label);

if length(selgrad)<length(coef.label)
  % in principle this could be less restrictive, but that requires more bookkeeping
  error('not all required channels were found in the gradiometer definition');
end

% translate the input type to a number
switch cfg.input
  case 'none'
    input = 0;
  case 'G1BR'
    input = 1;
  case 'G2BR'
    input = 2;
  case 'G3BR'
    input = 3;
  otherwise
    error('unknown balancing for input data');
end

% translate the output type to a number
switch cfg.output
  case 'none'
    output = 0;
  case 'G1BR'
    output = 1;
  case 'G2BR'
    output = 2;
  case 'G3BR'
    output = 3;
  otherwise
    error('unknown balancing for output data');
end

% compute the linear transformation matrix
if input==output
  error('you should select a different balancing for the input and output data');
elseif input==0 && output==1
  T = coef.t01;
elseif input==0 && output==2
  T = coef.t02;
elseif input==0 && output==3
  T = coef.t03;
elseif input==1 && output==0
    T =            inv(coef.t01);
elseif input==1 && output==2
    T = coef.t02 * inv(coef.t01);
elseif input==1 && output==3
    T = coef.t02 * inv(coef.t01);
elseif input==2 && output==0
    T =            inv(coef.t02);
elseif input==2 && output==1
    T = coef.t01 * inv(coef.t02);
elseif input==2 && output==3
    T = coef.t02 * inv(coef.t02);
elseif input==3 && output==0
    T =            inv(coef.t03);
elseif input==3 && output==1
    T = coef.t01 * inv(coef.t03);
elseif input==3 && output==2
    T = coef.t02 * inv(coef.t03);
end

% apply the transformatoin on the data
Ntrials = length(data.trial);
for i=1:Ntrials
  fprintf('processing trial %d from %d\n', i, Ntrials);
  data.trial{i}(seldat,:) = T*data.trial{i}(seldat,:);
end

% apply the transformation on the gradiometer definition
if input~=0
  error('there is a known bug in the code, sorry for that: only supoorted is 0th order input');
end
data.grad.tra(selgrad,:) = T*data.grad.tra(selgrad,:);

if isavgdata
  % reformat the output back into an average-like structure
  data.avg = data.trial{1};
  data = rmfield(data, 'trial');
  % remove those fields that were not updated
  try, data = rmfield(data, 'var'); end
  try, data = rmfield(data, 'cov'); end
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: syntheticgradient.m,v 1.8 2006/03/29 08:40:26 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
data.cfg = cfg;

