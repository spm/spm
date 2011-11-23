function varargout = pm_get_defaults(defstr, varargin)

% FORMAT defval = pm_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT spm_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit spm_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton
% $Id: pm_get_defaults.m 4572 2011-11-23 17:35:10Z chloe $
% ---------------------------------------------------------------------
global pm_def
if isempty(pm_def)
    pm_defaults;
    pm_def.sessname='session';
    pm_def.pedir=2;
elseif ~isfield(pm_def,'pedir');
    pm_def.pedir=2;
elseif ~isfield(pm_def,'sessname');
    pm_def.sessname='session';
end

[default_file_path, tmpname] = fileparts(mfilename('fullpath'));
pm_def.defaultsfilename{1} = sprintf('%s%s%s',default_file_path,filesep,'pm_defaults.m');
if ischar(pm_def.MFLAGS.TEMPLATE)
    pm_def.MFLAGS.TEMPLATE = {pm_def.MFLAGS.TEMPLATE};
end
% tmpTEMPLATE = pm_def.MFLAGS.TEMPLATE;
% pm_def.MFLAGS.TEMPLATE=[];
% pm_def.MFLAGS.TEMPLATE{1}=tmpTEMPLATE;

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(pm_def, subs);
else
    pm_def = subsasgn(pm_def, subs, varargin{1});
end

