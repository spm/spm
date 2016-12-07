function spm_cli(varargin)
% Command line interface for SPM
%__________________________________________________________________________
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cli.m 6964 2016-12-07 16:50:41Z guillaume $ 


%-Input arguments
%--------------------------------------------------------------------------
if isempty(varargin), return; end
if nargin == 1 && ismember(varargin{1},{'help','--help'})
    cmd = lower(spm('Ver'));
    fprintf([...
        'Usage: %s [NODE] [arg...]\n',...
        '       %s [NODE] help\n',...
        '       %s [ help | --help ]\n',...
        '       %s help --verbose\n',...
        '\n',...
        '[NODE]         Batch node name\n',...
        '--verbose      List all available batch nodes\n',...
        '--help, help   Print usage statement\n'],...
        cmd, cmd, cmd, cmd);
        return
end
if nargin == 2 ...
        && all(ismember(strtok(varargin(1:2),'-'),{'help','verbose'}))
    c0 = cfg_util('getcfg');
    nodes = find_all_modules(c0.values{1})';
    for i=1:numel(nodes)
        [n, c] = find_module(strsplit(nodes{i},'.'), c0);
        h = normalise_help_text(c.help);
        if ~c.hidden, fprintf('  %s    %s\n',nodes{i},h{1}); end
    end
    return
end
node = strsplit(varargin{1},'.');
opts = varargin(2:end);

%-Get the configuration tree
%--------------------------------------------------------------------------
try
    c0 = cfg_util('getcfg');
catch
    spm_jobman('initcfg');
    c0 = cfg_util('getcfg');
end

%-Search module in configuration tree and allow for node shortcut
% (use first matching module if ambiguity)
%--------------------------------------------------------------------------
[n, c0, found] = find_module(node, c0);

if ~found
    error('Cannot find module %s.\n',strjoin(node,'.'));
else
    n = n(2:end); % first is matlabbatch
    if ~isequal(n,node)
        node = n;
        fprintf('Matching module %s.\n',strjoin(node,'.'));
    end
end

%-Help for a specific node
%--------------------------------------------------------------------------
if numel(opts) == 1 && ismember(opts{1},{'help','--help'})
    cmd = lower(spm('Ver'));
    fprintf([...
        'Usage: %s %s [OPTIONS]\n',...
        '\n',...
        'Options:\n'],...
        cmd,strjoin(node,'.'));
    [flags, h1s, ms] = get_items(c0);
    n = cellfun(@numel,flags);
    for i=1:numel(flags)
        if ms(i), m = '(*) '; else m = ''; end
        fprintf('    --%s%s%s%s\n',flags{i},repmat(' ',1,3+max(n)-n(i)),m,h1s{i});
    end
    fprintf([...
        '\n',...
        '    --help%sPrint usage statement\n'],repmat(' ',1,3+max(n)-4));
    return;
end

%-Get options
%--------------------------------------------------------------------------
i = 1;
params = {};
values = {};
while i <= numel(opts)
    flag = strtok(opts{i},'-');
    opt = {};
    i = i + 1;
    for j=i:numel(opts)
        if opts{j}(1) == '-', break; end
        opt{end+1} = opts{j};
        i = i + 1;
    end
    params{end+1} = flag;
    values{end+1} = opt;
end

%-Create batch job structure
%--------------------------------------------------------------------------
job = struct();
for i=1:numel(params)
    prm = params{i};
    if ~any(prm == '.')
        subs = substruct('.',prm);
    else
        clear subs
        p = strsplit(prm,'.');
        for j=1:numel(p)
            subs(j) = substruct('.',p{j});
        end
        prm = p{end};
    end
    job = subsasgn(job, subs, get_val(values{i},prm,c0));
end

%-Execute batch job
%--------------------------------------------------------------------------
subs = reshape([repmat({'.'},1,numel(node));node],1,[]);
subs = substruct(subs{:});
job  = subsasgn(struct(),subs,job);
cfg_message('none','destination','^matlabbatch:run:jobstart$');
cfg_message('none','destination','^matlabbatch:run:jobdone$');
cfg_message('none','destination','^matlabbatch:run:jobfailed$');
cfg_message('none','destination','^matlabbatch:run:modstart$');
cfg_message('none','destination','^matlabbatch:run:moddone$');
%cfg_message('none','destination','^matlabbatch:run:modfailed$');
try
    [out, job] = spm_jobman('run', {job});
catch
    err = lasterror;
    switch err.identifier
        case 'spm:spm_jobman:jobNotFilled'
            error('Incomplete job.');
        case 'matlabbatch:run:jobfailederr'
            error('Job failed.');
        otherwise
            error('Unknown error during job execution.');
    end
end


%==========================================================================
function [node, c0, found] = find_module(node, c0, found)

if nargin < 3, found = false; end
if found, return; end
if strcmp(node{1}, c0.tag)
    if numel(node) == 1
        if isa(c0,'cfg_exbranch'), found = true; end
    else
        if ~isa(c0,'cfg_exbranch')
            for i=1:numel(c0.values)
                [n, c1 , found] = find_module(node(2:end), c0.values{i});
                if found, node = [c0.tag n]; c0 = c1; return; end
            end
        end
    end
else
    if ~isa(c0,'cfg_exbranch')
        for i=1:numel(c0.values)
            [n, c1, found] = find_module(node, c0.values{i});
            if found, node = [c0.tag n]; c0 = c1; return; end
        end
    end
end


%==========================================================================
function [nodes,n] = find_all_modules(c0,nodes,n)

if nargin < 2, nodes = {}; end
if nargin < 3, n = {}; end
if isa(c0,'cfg_exbranch')
    nodes = [nodes strjoin([n c0.tag],'.')];
else
    m = n;
    for i=1:numel(c0.values)
        [nodes, n] = find_all_modules(c0.values{i},nodes,[m c0.tag]);
    end
end


%==========================================================================
function param = get_val(param,flag,c0)

for i=1:numel(c0.val)
    if strcmp(c0.val{i}.tag, flag)
        switch class(c0.val{i})
            case 'cfg_files'
                param = param(:);
            case 'cfg_entry'
                switch c0.val{i}.strtype
                    case {'n','w','r'}
                        param = str2num(char(param));
                    case {'s'}
                        param = char(param);
                end
            case 'cfg_menu'
                if numel(param) > 1
                    error('Only one value accepted for "%s".',flag);
                end
                param = char(param);
                labels = c0.val{i}.labels;
                j = find(ismember(labels,param));
                if nnz(j)
                    param = c0.val{i}.values{j};
                else
                    opts = sprintf('  "%s"\n',labels{:});
                    error(['Cannot translate option "%s": %s.\n',...
                        'Options are:\n%s'],flag, param, opts);
                end
            otherwise
        end
        return;
    end
end
for i=1:numel(c0.val)
    switch class(c0.val{i})
        case 'cfg_branch'
            param = get_val(param,flag,c0.val{i});
        case 'cfg_repeat'
        case 'cfg_choice'
    end
end


%==========================================================================
function [flags, h1s, ms] = get_items(c0)
flags = {}; h1s = {}; ms = [];
for i=1:numel(c0.val) % go recursive instead
    flag = {}; h1 = {}; m = [];
    if isa(c0.val{i},'cfg_repeat')
        % if more than one, they are not stored yet
        for j=1:numel(c0.val{i}.values)
            [flag, h1, m] = get_items(c0.val{i}.values{j});
            flag = [{c0.val{i}.values{j}.tag} flag];
            h1   = [{char(c0.val{i}.values{j}.help{:})} h1];
            m    = [is_mandatory(c0.val{i}.values{j}) m];
        end
    elseif isa(c0.val{i},'cfg_branch')
        [flag, h1, m] = get_items(c0.val{i});
        flag = [{c0.val{i}.tag} flag];
        h1   = [{char(c0.val{i}.help{:})} h1];
        m    = [is_mandatory(c0.val{i}) m];
    else
        flag = c0.val{i}.tag;
        h1   = char(c0.val{i}.help{:});
        m    = is_mandatory(c0.val{i});
    end
    h1    = normalise_help_text(h1);
    flags = [flags flag];
    h1s   = [h1s h1];
    ms    = [ms m];
end


%==========================================================================
function m = is_mandatory(item)
m = isempty(item.val) && isempty(item.def);


%==========================================================================
function h1 = normalise_help_text(h1)
h1 = cellstr(h1);
for i=1:numel(h1)
    if ~isempty(h1{i}) && h1{i}(end) == '.', h1{i} = h1{i}(1:end-1); end
    if numel(h1{i}) > 70 % use "stty size" to get the number of columns
        h1{i} = [h1{i}(1:70-3) '...'];
    end
end
