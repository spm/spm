function varargout = spm_select(varargin)
% File selector
% This functionality has moved to cfg_getfile, the matlabbatch file
% selector. See 
%                help cfg_getfile
% for a description of the file selection interface. spm_select will
% just pass on any inputs to cfg_getfile and return all of its output
% arguments.
% cfg_getfile has all functionality of spm_select, except handling of
% virtual files. Management, filtering, selection of virtual files has
% been very slow for large numbers of files. Also, once virtual files
% were selected, they did not change if the inputs changed from which
% they were derived.
% To pass output items from one batch module to inputs of another one,
% matlabbatch now uses dependencies (see help on cfg_dep). These
% dependencies allow to pass any output variable to another
% modules at run time.
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_select.m 1143 2008-02-07 19:33:33Z spm $

if ~(exist('cfg_getfile') == 2)
    addpath(fullfile(spm('dir'),'matlabbatch'));
end;
[varargout{1:nargout}] = cfg_getfile(varargin{:});
