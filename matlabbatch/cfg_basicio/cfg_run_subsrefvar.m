function varargout = cfg_run_subsrefvar(cmd, varargin)
% Template function to implement callbacks for an cfg_exbranch. The calling
% syntax is
% varargout = cfg_run_subsrefvar(cmd, varargin)
% where cmd is one of
% 'run'      - out = cfg_run_subsrefvar('run', job)
%              Run a job, and return its output argument
% 'vout'     - dep = cfg_run_subsrefvar('vout', job)
%              Examine a job structure with all leafs present and return an
%              array of cfg_dep objects.
% 'check'    - str = cfg_run_subsrefvar('check', subcmd, subjob)
%              Examine a part of a fully filled job structure. Return an empty
%              string if everything is ok, or a string describing the check
%              error. subcmd should be a string that identifies the part of
%              the configuration to be checked.
% 'defaults' - defval = cfg_run_subsrefvar('defaults', key)
%              Retrieve defaults value. key must be a sequence of dot
%              delimited field names into the internal def struct which is
%              kept in function local_def. An error is returned if no
%              matching field is found.
%              cfg_run_subsrefvar('defaults', key, newval)
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
% 'run'      - @(job)cfg_run_subsrefvar('run', job)
% 'vout'     - @(job)cfg_run_subsrefvar('vout', job)
% 'check'    - @(job)cfg_run_subsrefvar('check', 'subcmd', job)
% 'defaults' - @(val)cfg_run_subsrefvar('defaults', 'defstr', val{:})
%              Note the list expansion val{:} - this is used to emulate a
%              varargin call in this function handle.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_subsrefvar.m 3948 2010-06-25 09:48:03Z volkmar $

rev = '$Rev: 3948 $'; %#ok

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job = local_getjob(varargin{1});
            % do computation, return results in variable out
            subs = local_getsubs(job, true);
            out.output = subsref(job.input, subs);
            if nargout > 0
                varargout{1} = out;
            end
        case 'vout'
            job = local_getjob(varargin{1});
            subscode = char(gencode_substruct(local_getsubs(job, false)));
            dep = cfg_dep;
            if isempty(subscode)
                dep.sname = 'Referenced part of variable';
            else
                dep.sname = sprintf('val%s', subscode);
            end
            dep.src_output = substruct('.','output');
            if isequal(job.tgt_spec,'<UNKNOWN>')
                dep.tgt_spec = cfg_findspec({{'strtype','e'}});
            else
                fn = fieldnames(job.tgt_spec);
                dep.tgt_spec = job.tgt_spec.(fn{1});
            end
            varargout{1} = dep;
        case 'check'
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    % implement checks, return status string in variable str
                    case 'subsind'
                        if (ischar(subjob) && isequal(subjob, ':')) || ...
                                (isnumeric(subjob) && isequal(subjob, round(subjob)) && all(subjob > 0))
                            str = '';
                        else
                            str = 'Subscript index must be either a vector of natural numbers or the character '':''.';
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

function subs = local_getsubs(job, cwarn)
% generate subscript structure
% generate warning for cell references if requested
subs = struct('type',{},'subs',{});
for k = 1:numel(job.subsreference)
    switch char(fieldnames(job.subsreference{k}))
        case 'subsfield'
            subs(k).type = '.';
            subs(k).subs = job.subsreference{k}.subsfield;
        case 'subsinda'
            subs(k).type = '()';
            subs(k).subs = job.subsreference{k}.subsinda;
        case 'subsindc'
            subs(k).type = '{}';
            subs(k).subs = job.subsreference{k}.subsindc;
            if cwarn && any(cellfun(@(x)(isequal(x,':')||numel(x)>1), ...
                    job.subsreference{k}.subsindc))
                cfg_message('cfg_basicio:subsrefvar', ...
                    'Trying to access multiple cell elements - only returning first one.');
            end
    end
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