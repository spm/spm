function varargout = cfg_run_call_matlab(cmd, varargin)
% A generic interface to call any MATLAB function through the batch system
% and make its output arguments available as dependencies.
% varargout = cfg_run_call_matlab(cmd, varargin)
% where cmd is one of
% 'run'      - out = cfg_run_call_matlab('run', job)
%              Run the function, and return the specified output arguments
% 'vout'     - dep = cfg_run_call_matlab('vout', job)
%              Return dependencies as specified via the output cfg_repeat.
% 'check'    - str = cfg_run_call_matlab('check', subcmd, subjob)
%              Examine a part of a fully filled job structure. Return an empty
%              string if everything is ok, or a string describing the check
%              error. subcmd should be a string that identifies the part of
%              the configuration to be checked.
% 'defaults' - defval = cfg_run_call_matlab('defaults', key)
%              Retrieve defaults value. key must be a sequence of dot
%              delimited field names into the internal def struct which is
%              kept in function local_def. An error is returned if no
%              matching field is found.
%              cfg_run_call_matlab('defaults', key, newval)
%              Set the specified field in the internal def struct to a new
%              value.
% Application specific code needs to be inserted at the following places:
% 'run'      - main switch statement: code to compute the results, based on
%              a filled job
% 'vout'     - main switch statement: code to compute cfg_dep array, based
%              on a job structure that has all leafs, but not necessarily
%              any values filled in
% 'check'    - create and populate switch subcmd switchyard
% 'defaults' - modify initialisation of defaults in subfunction local_defs
% Callbacks can be constructed using anonymous function handles like this:
% 'run'      - @(job)cfg_run_call_matlab('run', job)
% 'vout'     - @(job)cfg_run_call_matlab('vout', job)
% 'check'    - @(job)cfg_run_call_matlab('check', 'subcmd', job)
% 'defaults' - @(val)cfg_run_call_matlab('defaults', 'defstr', val{:})
%              Note the list expansion val{:} - this is used to emulate a
%              varargin call in this function handle.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_call_matlab.m 3810 2010-04-07 12:42:32Z volkmar $

rev = '$Rev: 3810 $'; %#ok

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job = local_getjob(varargin{1});
            % do computation, return results in variable out
            in  = cell(size(job.inputs));
            for k = 1:numel(in)
                in{k} = job.inputs{k}.(char(fieldnames(job.inputs{k})));
            end
            out.outputs = cell(size(job.outputs));
            [out.outputs{:}] = feval(job.fun, in{:});
            % make sure output filenames are cellstr arrays, not char
            % arrays
            for k = 1:numel(out.outputs)
                if isfield(job.outputs{k},'filter') && isa(out.outputs{k},'char')
                    out.outputs{k} = cellstr(out.outputs{k});
                end
            end
            if nargout > 0
                varargout{1} = out;
            end
        case 'vout'
            job = local_getjob(varargin{1});
            % initialise empty cfg_dep array
            dep = cfg_dep;
            dep = dep(false);
            % determine outputs, return cfg_dep array in variable dep
            for k = 1:numel(job.outputs)
                dep(k)            = cfg_dep;
                dep(k).sname      = sprintf('Call MATLAB: output %d - %s %s', k, char(fieldnames(job.outputs{k})), char(fieldnames(job.outputs{k}.(char(fieldnames(job.outputs{k}))))));
                dep(k).src_output = substruct('.','outputs','{}',{k});
                dep(k).tgt_spec   = cfg_findspec({{'strtype','e', char(fieldnames(job.outputs{k})), char(fieldnames(job.outputs{k}.(char(fieldnames(job.outputs{k})))))}});
            end
            varargout{1} = dep;
        case 'check'
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    % implement checks, return status string in variable str
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