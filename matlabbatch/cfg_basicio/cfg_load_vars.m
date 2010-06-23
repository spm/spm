function varargout = cfg_load_vars(cmd, varargin)
% Load a .mat file, and return its contents via output dependencies.
% varargout = cfg_load_vars(cmd, varargin)
% where cmd is one of
% 'run'      - out = cfg_load_vars('run', job)
%              Run a job, and return its output argument
% 'vout'     - dep = cfg_load_vars('vout', job)
%              Create a virtual output for each requested variable. If
%              "all variables" are requested, only one output will be
%              generated.
% 'check'    - str = cfg_load_vars('check', subcmd, subjob)
%              'isvarname' - check whether the entered string is a valid
%                            MATLAB variable name. This does not check
%                            whether the variable is present in the .mat file.
% 'defaults' - defval = cfg_load_vars('defaults', key)
%              No defaults.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_load_vars.m 3944 2010-06-23 08:53:40Z volkmar $

rev = '$Rev: 3944 $'; %#ok

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job = local_getjob(varargin{1});
            % do computation, return results in variable out
            var = load(job.matname{1});
            if isfield(job.loadvars,'allvars')
                out{1} = var;
            else
                out = cell(size(job.loadvars.varname));
                for k = 1:numel(job.loadvars.varname)
                    try
                        out{k} = var.(job.loadvars.varname{k});
                    catch
                        error(['Variable ''%s'' could not be loaded from ' ...
                               'file ''%s''.'], job.loadvars.varname{k}, ...
                              job.matname{1});
                    end
                end
            end
            if nargout > 0
                varargout{1} = out;
            end
        case 'vout'
            job = local_getjob(varargin{1});
            % initialise empty cfg_dep array
            if isfield(job.loadvars,'allvars')
                dep = cfg_dep;
                dep(1).sname = 'Loaded Variables (struct)';
                dep(1).src_output = substruct('{}',{1});
                dep(1).tgt_spec   = cfg_findspec({{'class','cfg_entry'}});
            else
                for k = 1:numel(job.loadvars.varname)
                    dep(k) = cfg_dep;
                    if ~isempty(job.loadvars.varname{k}) && ...
                            ischar(job.loadvars.varname{k}) && ~ ...
                            strcmpi(job.loadvars.varname{k}, ...
                                    '<UNDEFINED>')
                        dep(k).sname = sprintf('Loaded Variable ''%s''', ...
                                               job.loadvars.varname{k});
                    else
                        dep(k).sname = sprintf('Loaded Variable #%d', k);
                    end
                    dep(k).src_output = substruct('{}',{k});
                    dep(k).tgt_spec   = cfg_findspec({{'class','cfg_entry'}});
                end
            end
            % determine outputs, return cfg_dep array in variable dep
            varargout{1} = dep;
        case 'check'
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    case 'isvarname'
                        if ~isvarname(subjob)
                            str = sprintf(['''%s'' is not a valid MATLAB ' ...
                                           'variable name'], subjob);
                        end
                    otherwise
                        cfg_message('unknown:check', ...
                            'Unknown check subcmd ''%s''.', subcmd);
                end
                varargout{1} = str;
            else
                cfg_message('ischar:check', 'Subcmd must be a string.');
            end
        case 'defaults'
            if nargin == 2
                varargout{1} = local_defs(varargin{1});
            else
                local_defs(varargin{1:2});
            end
        otherwise
            cfg_message('unknown:cmd', 'Unknown command ''%s''.', cmd);
    end
else
    cfg_message('ischar:cmd', 'Cmd must be a string.');
end

function varargout = local_defs(defstr, defval)
persistent defs;
if isempty(defs)
    % initialise defaults
end
if ischar(defstr)
    % construct subscript reference struct from dot delimited tag string
    tags = textscan(defstr,'%s', 'delimiter','.');
    subs = struct('type','.','subs',tags{1}');
    try
        cdefval = subsref(local_def, subs);
    catch
        cdefval = [];
        cfg_message('defaults:noval', ...
            'No matching defaults value ''%s'' found.', defstr);
    end
    if nargin == 1
        varargout{1} = cdefval;
    else
        defs = subsasgn(defs, subs, defval);
    end
else
    cfg_message('ischar:defstr', 'Defaults key must be a string.');
end

function job = local_getjob(job)
if ~isstruct(job)
    cfg_message('isstruct:job', 'Job must be a struct.');
end