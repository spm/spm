function [cfg, artifact] = artifact_file(cfg);

% ARTIFACT_FILE reads rejection marks from a file
% 
% See also REJECTARTIFACT

% Undocumented local options:
% cfg.artfctdef
% cfg.headerfile
% cfg.rejectfile
% cfg.trl

% Copyright (C) 2003-2006, Robert Oostenveld
%
% $Log: artifact_file.m,v $
% Revision 1.10  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.9  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.8  2006/01/31 13:49:29  jansch
% included dataset2files to ensure the presence of cfg.headerfile or cfg.datafile
% whenever needed
%
% Revision 1.7  2005/12/20 08:36:47  roboos
% add the artifact Nx2 matrix to the output configuration
% changed some indentation and white space, renamed a few variables
%
% Revision 1.6  2005/06/29 12:42:00  roboos
% added version to the output configuration
%
% Revision 1.5  2005/05/17 17:50:36  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.4  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%
% Revision 1.3  2004/08/26 16:01:34  roboos
% added support for brainvision_marker
%
% Revision 1.2  2003/12/11 09:35:27  roberto
% fixed bug in "if isfield()"
%
% Revision 1.1  2003/12/08 12:34:27  roberto
% initial version, according to pluggable artifact rejection
%

if isfield(cfg, 'rejectfile') && ~strcmp(cfg.rejectfile, 'no')
  cfg = dataset2files(cfg);
  hdr = read_fcdc_header(cfg.headerfile);
  if filetype(cfg.rejectfile, 'eep_rej')
    artifact = read_eep_rej(cfg.rejectfile);
  elseif filetype(cfg.rejectfile, 'brainvision_marker')
    artifact = read_brainvision_marker(cfg.rejectfile);
  else
    error(sprintf('unrecognized filetype for rejection file ''%s''', cfg.rejectfile));
  end
  % convert rejection marks from miliseconds into seconds
  artifact = artifact/1000;
  % convert rejection marks into samples
  artifact = round(artifact * hdr.Fs);
  fprintf('%d rejection marks read\n', size(artifact,1));
else
  artifact = [];
  fprintf('no rejection marks read\n');
end

% remember the details that were used here
cfg.artfctdef.file.trl      = cfg.trl;
cfg.artfctdef.file.artifact = artifact;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: artifact_file.m,v 1.10 2006/05/03 08:12:51 ingnie Exp $';

