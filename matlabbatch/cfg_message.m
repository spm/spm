function cfg_message(varargin)

% function cfg_message(msgid, msgfmt, varargin)
% Display a message. The message identifier msgid will be looked up in a
% message database to decide how to treat this message. This database is
% a struct array with fields:
% .identifier  - message id
% .level       - message severity level. One of
%                'info'    - print message
%                'warning' - print message, raise a warning
%                'error'   - print message, throw an error
% .destination - output destination. One of
%                'none'    - silently ignore this message
%                'stdout'  - standard output
%                'stderr'  - standard error output
%                'syslog'  - (UNIX) syslog
%                Warnings and errors will always be logged to the command
%                window and to syslog, if destination == 'syslog'. All
%                other messages will only be logged to the specified location.
% .verbose
% .backtrace   - control verbosity and backtrace, one of 'on' or 'off'
%
% function cfg_message('on'|'off', 'verbose'|'backtrace', msgidregexp)
% Set verbosity and backtrace display for all messages where msgid
% matches msgidregexp. To match a message id exactly, use the regexp
% '^msgid$'.
%
% function cfg_message('none'|'stdout'|'stderr'|'syslog', 'destination', msgidregexp)
% Set destination for all messages matching msgidregexp.
%
% function cfg_message('info'|'warning'|'error', msgidregexp)
% Set severity level for all messages matching msgidregexp.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_message.m 1873 2008-07-01 15:23:40Z volkmar $

rev = '$Rev: 1873 $'; %#ok

if nargin < 1 || isempty(varargin{1})
    return;
end

% Get message settings
msgcfg = cfg_get_defaults('msgcfg');
msgtpl = cfg_get_defaults('msgtpl');
msgdef = cfg_get_defaults('msgdef');
% Helper inline functions
ismsgid = inline('~isempty(regexp(msgidstr,''^([a-zA-Z]\w*:)+[a-zA-Z]\w*$'',''once''))',...
                 'msgidstr');
findcfg = inline('~cellfun(@isempty,regexp({msgcfg.identifier},msgregexp,''once''))', ...
                 'msgcfg','msgregexp');
% Input checks
if ~(ischar(varargin{1}) && size(varargin{1},1) == 1) && ...
        ~(isstruct(varargin{1}) && numel(varargin{1}) == 1 && ...
          all(isfield(varargin{1}, {'message', 'identifier'})))
    cfg_message('matlabbatch:cfg_message', ...
          'First argument must be a one-line character string or an errorstruct.');
    return;
end
if any(strcmpi(varargin{1}, {'on', 'off', 'none', 'stdout', 'stderr', 'syslog'}))
    % Set message properties
    if nargin ~= 3
        cfg_message('matlabbatch:cfg_message', ...
              'Must specify status, property and msgidregexp.');
        return;
    end
    if any(strcmpi(varargin{1}, {'on', 'off'})) && ~any(strcmpi(varargin{2}, ...
                    {'verbose', 'backtrace'}))
        cfg_message('matlabbatch:cfg_message', ...
              ['Message property must be one of ''verbose'' or ' ...
               '''backtrace''.']);
        return;
    elseif any(strcmpi(varargin{1}, {'none', 'stdout', 'stderr', 'syslog'})) && ...
            ~strcmpi(varargin{2}, 'destination')
        cfg_message('matlabbatch:cfg_message', ...
              'Message property must be ''destination''.');
        return;
    end
    if ~ischar(varargin{3}) || size(varargin{3},1) ~= 1
        cfg_message('matlabbatch:cfg_message', ...
              'Third argument must be a one-line character string.');
        return;
    end
    msgstate  = lower(varargin{1});
    msgprop   = lower(varargin{2});
    msgregexp = varargin{3};
    sel = findcfg(msgcfg, msgregexp);
    if any(sel)
        % Set property on all matching messages
        [msgcfg(sel).(msgprop)] = deal(msgstate);
        cfg_get_defaults('msgcfg',msgcfg);
    elseif ismsgid(msgregexp)
        % Add new rule, if msgregexp is a valid id
        msgcfg(end+1)          = msgdef;
        msgcfg(end).identifier = msgregexp;
        msgcfg(end).(msgprop)  = msgstate;
        cfg_get_defaults('msgcfg',msgcfg);
    end
    if ~ismsgid(msgregexp)
        % Update templates
        sel = strcmp(msgregexp, {msgtpl.identifier});
        if any(sel)
            % Update literally matching regexps
            [msgtpl(sel).(msgprop)] = deal(msgstate);
        else
            % Add new template rule
            msgtpl(end+1)          = msgdef;
            msgtpl(end).identifier = msgregexp;
            msgtpl(end).(msgprop)  = msgstate;
        end
        cfg_get_defaults('msgtpl',msgtpl);
    end
elseif any(strcmpi(varargin{1}, {'info', 'warning', 'error'}))
    % Set message level
    if nargin ~= 2
        cfg_message('matlabbatch:cfg_message', ...
              'Must specify level and msgidregexp.');
        return;
    end
    if ~ischar(varargin{2}) || size(varargin{2},1) ~= 1
        cfg_message('matlabbatch:cfg_message', ...
              'Second argument must be a one-line character string.');
        return;
    end
    msglvl    = lower(varargin{1});
    msgregexp = varargin{2};
    sel = findcfg(msgcfg, msgregexp);
    if any(sel)
        % Set message level
        [msgcfg(sel).level] = deal(msglvl);
        cfg_get_defaults('msgcfg',msgcfg);
    elseif ismsgid(msgregexp)
        % Add new rule, if msgregexp is a valid id
        msgcfg(end+1)          = msgdef;
        msgcfg(end).identifier = msgregexp;
        msgcfg(end).level      = msglvl;
        cfg_get_defaults('msgcfg',msgcfg);
    end
    if ~ismsgid(msgregexp)
        % Update templates
        sel = strcmp(msgregexp, {msgtpl.identifier});
        if any(sel)
            % Update literally matching regexps
            [msgtpl(sel).level] = deal(msglvl);
        else
            % Add new template rule
            msgtpl(end+1)          = msgdef;
            msgtpl(end).identifier = msgregexp;
            msgtpl(end).level      = msglvl;
        end
        cfg_get_defaults('msgtpl',msgtpl);
    end
else
    % Issue message
    if isstruct(varargin{1})
        msgid  = varargin{1}.identifier;
        msgfmt = varargin{1}.message;
        extras = {};
        if isfield(varargin{1},'stack')
            stack = varargin{1}.stack;
        else
            stack = dbstack(2);
        end
        % discard other fields
    else
        if nargin < 2
            cfg_message('matlabbatch:cfg_message', ...
                  'Must specify msgid and message.');
        end
        if ~ischar(varargin{2}) || size(varargin{2},1) ~= 1
            cfg_message('matlabbatch:cfg_message', ...
                  'Second argument must be a one-line character string.');
        end
        msgid  = varargin{1};
        msgfmt = varargin{2};
        extras = varargin(3:end);
        stack = dbstack(2);
    end
    % find msgcfg entry for msgid
    sel = strcmp(msgid, {msgcfg.identifier});
    if any(sel)
        cmsgcfg = msgcfg;
    else
        % no identity match found, match against regexp templates
        sel = ~cellfun(@isempty, regexp(msgid, {msgtpl.identifier}));
        cmsgcfg = msgtpl;
        [cmsgcfg(sel).identifier] = deal(msgid);
    end
    switch nnz(sel)
        case 0
            % no matching id (or invalid msgid), use default setting
            cmsgcfg = msgdef;
            cmsgcfg.identifier = msgid;
            issuemsg(cmsgcfg, msgfmt, extras);
        case 1
            % exact match
            issuemsg(cmsgcfg(sel), msgfmt, extras);
        otherwise
            % multiple matches, use last match (i.e. least recently added rule)
            is = find(sel);
            issuemsg(cmsgcfg(is(end)), msgfmt, extras);
    end
end

function issuemsg(msgcfg, msgfmt, extras)
if strcmp(msgcfg.level,'error') || ~strcmp(msgcfg.destination, 'none')
    switch msgcfg.destination
        case 'none'
            fid = -1;
        case 'stdout'
            fid = 1;
        case 'stderr'
            fid = 2;
        case 'syslog'
            if isunix
                fid = '/usr/bin/logger';
            else
                fid = -1;
            end
    end
    msg = sprintf(msgfmt, extras{:});
    switch msgcfg.level
        case 'info'
            if ischar(fid)
                unix(sprintf('%s -t [MATLAB] %s: %s', fid, msgcfg.identifier, ...
                             msg));
            elseif ~ischar(fid) && fid > 0
                % only display, if destination is neither syslog nor none
                fprintf(fid, '%s\n', msg);
                if strcmp(msgcfg.backtrace, 'on')
                    dbstack(2) % This will ignore fid
                end
                if strcmp(msgcfg.verbose, 'on')
                    fprintf(fid, 'To turn off this message, type %s(''none'', ''destination'', ''%s'').\n',...
                            mfilename, msgcfg.identifier);
                end
            end
        case 'warning'
            if ischar(fid)
                unix(sprintf('%s -t [MATLAB] %s: %s', fid, msgcfg.identifier, ...
                             msg));
            elseif ~ischar(fid) && fid > 0
                % only display, if destination is neither syslog nor none
                bsts = warning('query','backtrace'); 
                vsts = warning('query','verbose');
                warning('off', 'backtrace');
                warning('off', 'verbose');
                warning(msgcfg.identifier, msg);
                if strcmp(msgcfg.backtrace, 'on')
                    dbstack(2)
                end
                if strcmp(msgcfg.verbose, 'on')
                    fprintf('To turn off this message, type %s(''none'', ''destination'', ''%s'').\n',...
                            mfilename, msgcfg.identifier);
                end
                warning(bsts);
                warning(vsts);
            end
        case 'error'
            if ischar(fid)
                unix(sprintf('%s -t [MATLAB] %s: %s', fid, msgcfg.identifier, ...
                             msg));
            end
            le.identifier = msgcfg.identifier;
            le.message    = msg;
            le.stack      = dbstack(2);
            error(le);
    end
end