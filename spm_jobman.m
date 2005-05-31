function varargout = spm_jobman(varargin)
% UI/Batching stuff
%_______________________________________________________________________
% This code is based on an earlier version by Philippe Ciuciu and
% Guillaume Flandin of Orsay, France.
%
% FORMAT spm_jobman
%        spm_jobman('interactive')
%        spm_jobman('interactive',job)
%        spm_jobman('interactive',job,node)
%        spm_jobman('interactive','',node)
% Runs the user interface in interactive mode.
%
% FORMAT spm_jobman('serial')
%        spm_jobman('serial',job)
%        spm_jobman('serial',job,node)
%        spm_jobman('serial','',node)
% Runs the user interface in serial mode.
%
% FORMAT spm_jobman('run',job)
% Runs a job.
%
% FORMAT spm_jobman('help',node)
%        spm_jobman('help',node,width)
% Creates a cell array containing help information.  This is justified
% to be 'width' characters wide. e.g.
%     h = spm_jobman('help','jobs.spatial.coreg.estimate');
%     for i=1:numel(h),fprintf('%s\n',h{i}); end;
%
%     node - indicates which part of the configuration is to be used.
%            For example, it could be 'jobs.spatial.coreg'.
%
%     job  - can be the name of a jobfile (as a .mat or a .xml), or a
%            'jobs' variable loaded from a jobfile.
%
%
% FORMAT spm_jobman('defaults')
% Runs the interactive defaults editor.
%
% FORMAT spm_jobman('pulldown')
% Creates a pulldown 'TASKS' menu in the Graphics window.
%
% FORMAT spm_jobman('chmod')
% Changes the modality for the TASKS pulldown.
% 
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_jobman.m 184 2005-05-31 13:23:32Z john $


if nargin==0
    setup_ui;
else
    switch lower(varargin{1})
    case {'interactive'}
        if nargin>=2
            setup_ui(varargin{2:nargin});
        else
            setup_ui;
        end;

    case {'serial'}
        if nargin>=2,
            serial(varargin{2:nargin});
        else
            serial;
        end;

    case {'run'}
        if nargin<2
            error('Nothing to run');
        end;
        run_job(varargin{2});

    case {'defaults'},
        setup_ui('defaults');

    case {'pulldown'}
       pulldown;

    case {'chmod'}
        if nargin>1,
            chmod(varargin{2});
        end;

    case {'help'}
        if nargin>=2,
            varargout{1} = showdoc(varargin{2:nargin});
        else
            varargout{1} = showdoc;
        end;

    otherwise
        error(['"' varargin{1} '" - unknown option']);
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function defaults_edit(varargin)
setup_ui('defaults');
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function interactive(varargin)
ud = get(varargin{1},'UserData');
if iscell(ud)
    setup_ui(ud{:});
else
    setup_ui;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function setup_ui(varargin)
% Configure the user interface

fg   = clearwin;
fs   = getdef('ui.fs');
if numel(fs)~=1 || ~isnumeric(fs{1}) || numel(fs{1})~=1,
    fs = 12;
else
    fs = fs{1};
end;

col1 = getdef('ui.colour1');
if numel(col1)~=1 || ~isnumeric(col1{1}) || numel(col1{1})~=3,
    col1 = [0.8 0.8 1];
else
    col1 = col1{1};
end;

col2 = getdef('ui.colour2');
if numel(col2)~=1 || ~isnumeric(col2{1}) || numel(col2{1})~=3,
    col2 = [1 1 0.8];
else
    col2 = col2{1};
end;

col3 = getdef('ui.colour3');
if numel(col3)~=1 || ~isnumeric(col3{1}) || numel(col3{1})~=3,
    col3 = [0 0 0];
else
    col3 = col3{1};
end;

t=uicontrol(fg,...
    'Style','listbox',...
    'Units','normalized',...
    'Position',[0.02 0.31 0.48 0.67],...
    'Callback',@click_batch_box,...
    'Tag','batch_box',...
    'ForegroundColor', col3,...
    'BackgroundColor',col1,...
    'FontSize',fs,...
    'Value',1);
c0 = cntxtmnu(t);
c1 = uimenu('Label','Exp/Con All',  'Parent',c0);
uimenu('Label','Expand All',   'Parent',c1,'Callback',@expandall);
uimenu('Label','Contract All', 'Parent',c1,'Callback',@contractall);

t=uicontrol(fg,...
    'Style','listbox',...
    'ListBoxTop',1,...
    'Units','normalized',...
    'Position',[0.02 0.02 0.96 0.28],...
    'Tag','help_box',...
    'FontName','fixedwidth',...
    'FontSize',fs,...
    'ForegroundColor', col3,...
    'BackgroundColor',col2);
set(t,'Value',[], 'Enable', 'inactive', 'Max',2, 'Min',0);
workaround(t);
cntxtmnu(t);

t=uicontrol(fg,...
    'Style','listbox',...
    'Units','normalized',...
    'Position',[0.51 0.74 0.47 0.24],...
    'ForegroundColor', col3,...
    'BackgroundColor',col1,...
    'FontSize',fs,...
    'Tag','opts_box');
cntxtmnu(t);

t=uicontrol(fg,...
    'Style','listbox',...
    'Units','normalized',...
    'Position',[0.51 0.49 0.47 0.24],...
    'Tag','val_box',...
    'ForegroundColor', col3,...
    'BackgroundColor',col2,...
    'Enable', 'inactive',...
    'Value',[],...
    'FontSize',fs,...
    'Max',2, 'Min',0);
set(t,'Value',[], 'Enable', 'on', 'Max',2, 'Min',0,'ListBoxTop',1);
cntxtmnu(t);

t=uicontrol(fg,...
    'Style','listbox',...
    'Units','normalized',...
    'Position',[0.51 0.38 0.47 0.10],...
    'Tag','msg_box',...
    'ForegroundColor', col3,...
    'BackgroundColor',col2,...
    'Enable', 'inactive',...
    'Value',[],...
    'FontSize',fs,...
    'Max',2, 'Min',0);
cntxtmnu(t);

if ~(nargin==1 && ischar(varargin{1}) && strcmp(varargin{1},'defaults')),
    t=uicontrol(fg,...
        'style','pushbutton',...
        'units','normalized',...
        'Position',[0.51 0.31 0.15 0.06],...
        'ForegroundColor', col3,...
        'BackgroundColor',col1,...
        'String','Save',...
        'Callback',@save_job,...
        'FontSize',fs,...
        'Tag','save');
    cntxtmnu(t);

    t=uicontrol(fg,...
        'style','pushbutton',...
        'units','normalized',...
        'Position',[0.67 0.31 0.15 0.06],...
        'ForegroundColor', col3,...
        'BackgroundColor',col1,...
        'String','Load',...
        'Callback',@load_job,...
        'FontSize',fs,...
        'Tag','load');
    cntxtmnu(t);

    t=uicontrol(fg,...
        'style','pushbutton',...
        'units','normalized',...
        'Position',[0.83 0.31 0.15 0.06],...
        'ForegroundColor', col3,...
        'BackgroundColor',col1,...
        'String','Run',...
        'Callback',@run_struct,...
        'Tag','run',...
        'FontSize',fs,...
        'Enable', 'off');
    cntxtmnu(t);
    if nargin>0
        initialise(varargin{:});
    else
        initialise;
    end;
else
    t=uicontrol(fg,...
        'style','pushbutton',...
        'units','normalized',...
        'Position',[0.51 0.31 0.15 0.06],...
        'ForegroundColor', col3,...
        'BackgroundColor',col1,...
        'String','OK',...
        'FontSize',fs,...
        'Callback',@ok_defaults);
    cntxtmnu(t);

    t=uicontrol(fg,...
        'style','pushbutton',...
        'units','normalized',...
        'Position',[0.67 0.31 0.15 0.06],...
        'ForegroundColor', col3,...
        'BackgroundColor',col1,...
        'String','Reset',...
        'FontSize',fs,...
        'Callback',@reset_defaults);
    cntxtmnu(t);

    t=uicontrol(fg,...
        'style','pushbutton',...
        'units','normalized',...
        'Position',[0.83 0.31 0.15 0.06],...
        'ForegroundColor', col3,...
        'BackgroundColor',col1,...
        'String','Cancel',...
        'FontSize',fs,...
        'Callback',@clearwin);
    cntxtmnu(t);
    initialise('defaults');
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function ok_defaults(varargin)
harvest_def(get(batch_box,'UserData'));
clearwin;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function reset_defaults(varargin)
global defaults

modality = [];
if isfield(defaults,'modality'),
    modality = defaults.modality;
end;
spm_defaults;
if ~isfield(defaults,'modality') && ~isempty(modality),
    defaults.modality = modality;
end;
pulldown;
clearwin;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function fg = clearwin(varargin)

fg = spm_figure('findwin','Graphics');
if isempty(fg), fg = spm_figure('Create','Graphics'); end;
delete(findobj(fg,'Parent',fg));
delete(batch_box);
delete(opts_box);
spm('Pointer');
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function expandall(varargin)
c         = expcon(get(gco,'UserData'),1);
str       = get_strings(c);
val       = min(get(gco,'Value'),length(str));
set(gco,'String',str,'Value',val,'UserData',c);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function contractall(varargin)
c         = expcon(get(gco,'UserData'),0);
str       = get_strings(c);
val       = min(get(gco,'Value'),length(str));
set(gco,'String',str,'Value',val,'UserData',c);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = expcon(c,val)
if isfield(c,'expanded'), c.expanded = val; end;
if isfield(c,'val'),
    for i=1:length(c.val),
        c.val{i} = expcon(c.val{i},val);
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function initialise(job,node)
% load the config file, possibly adding a job
% to it, and generally tidy it up.  The batch box
% is updated to contain the current structure.
if nargin<1, job = '';      end;
if nargin<2, node = 'jobs'; end;
c = initialise_struct(job);
set(batch_box,'UserData',start_node(c,node));
update_ui;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function run_job(job)
% Run a job. This function is not called via the UI.
c = initialise_struct(job);
if all_set(c),
    run_struct1(c);
else
    error('This job is not ready yet');
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function click_batch_box(varargin)
% Called when the batch box is clicked

remove_string_box;
if strcmp(get(get(varargin{1},'Parent'),'SelectionType'),'open')
    run_in_current_node(@expand_contract,false);
    str       = get_strings(get(batch_box,'UserData'));
    val       = min(get(batch_box,'Value'),length(str));
    set(batch_box,'String',str,'Value',val);
else
    run_in_current_node(@click_batch_box_fun,false);
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function update_ui
[str,sts] = get_strings(get(batch_box,'UserData'));
val       = min(get(batch_box,'Value'),length(str));
set(batch_box,'String',str,'Value',val);
run_in_current_node(@click_batch_box_fun,false);

run_but   = findobj(0,'Tag','run');
if sts
    set(run_but,'Enable', 'on');
else
    set(run_but,'Enable', 'off');
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [str,sts] = get_strings(c)
[unused,str,sts] = start_node(c,@get_strings1,0);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,str,sts] = get_strings1(c,l)
% Return the cell vector of strings displayed in the batch
% box, and also whether the job is runable.  This is a recursive
% function that a branch of the data structure is passed to.

if isfield(c,'hidden'),
    sts = 1;
    str = '';
    return;
end;

sts = 1;
mrk = '  ';
nam = '';
if isfield(c,'datacheck') && ~isempty(c.datacheck),
    mrk = ' <-@';
    sts = 0;
end;
if isfield(c,'name'), nam = c.name; end;

if isfield(c,'expanded') %  && isfield(c,'val') && ~isempty(c.val)

    if strcmp(c.type,'repeat') && isfield(c,'num')
        sts = (length(c.val) >= c.num(1)) &&...
              (length(c.val) <= c.num(2));
        if ~sts
            mrk = ' <-X';
        end;
    end;

    if c.expanded
        str = {[repmat('.   ',1,l) '-' nam]};
        for i=1:length(c.val)
            [c.val{i},str1,sts1] = get_strings1(c.val{i},l+1);
            sts = sts && sts1;
            if ~isempty(str1),
                str = {str{:},str1{:}};
            end;
        end;
    else
        str = {[repmat('.   ',1,l) '+' nam]};
        for i=1:length(c.val)
            [c.val{i},unused,sts1] = get_strings1(c.val{i},l+1);
            sts = sts && sts1;
        end;
        if ~sts,
            mrk    = ' <-X';
        end;
    end;
    str{1} = [str{1} mrk];
else
    switch c.type
    case {'files','menu','entry','const','choice'}
        if isempty(c.val)
            mrk = ' <-X';
            sts = 0;
        end;
    case 'repeat' % don't know, whether this ever gets executed
        if isfield(c,'num')
            sts = (length(c.val) >= c.num(1)) &&...
                  (length(c.val) <= c.num(2));
            if ~sts
                mrk = ' <-X';
            end;
        end;
    end;
    str = {[repmat('.   ',1,l) ' ' nam mrk]};
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function ok = all_set(c)
% Return whether all the necessary fields have been filled in.
% This is a recursive function that a branch of the data structure
% is passed to.

if isfield(c,'hidden'),
    ok = true;
    return;
end;

ok = true;

switch c.type
case {'files','menu','entry','const','choice'}
    if isempty(c.val)
        ok = false;
    end;

case {'branch','repeat'}
    if strcmp(c.type,'repeat') && isfield(c,'num')
        ok = ok && (length(c.val) >= c.num(1)) && ...
                   (length(c.val) <= c.num(2));
    end;
    if ok,
        for i=1:length(c.val),
            ok1 = all_set(c.val{i});
            ok  = ok && ok1;
            if ~ok, break; end;
        end;
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function click_options_box(varargin)
% The function called when the options box is clicked

dat = get(opts_box,'UserData');
if ~isempty(dat)
    fun = dat{get(opts_box,'Value')};
    run_in_current_node(fun.fun,true,fun.args{:});
    if fun.redraw, update_ui; end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function varargout = run_in_current_node(varargin)
% Get the current node, and run a job in it
varargout = {};
c   = get(batch_box,'UserData');
n   = get(batch_box,'Value');
show_msg('');
va = {c,@run_in_current_node1,n,varargin{:}};
[c,unused,varargout{:}] = start_node(va{:});
set(batch_box,'UserData',c);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,varargout] = start_node(c,fun,varargin)
persistent node;
if isempty(node)
    node = {'jobs'};
end;
varargout = {};
if nargin==2
    tmp = [0 find([fun '.']=='.')];
    node = {};
    for i=1:length(tmp)-1
        node = {node{:},fun((tmp(i)+1):(tmp(i+1)-1))};
    end;
    c = make_nodes(c,node);
    return;
end;
va = {c,node,fun,varargin{:}};
[c,unused,varargout{1:nargout-1}] = start_node1(va{:});
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = make_nodes(c,node)
if isfield(c,'tag') && ~strcmp(gettag(c),node{1})
    return;
end;
if isfield(c,'tag')
    node = {node{2:end}};
    if isempty(node), return; end;
end;
if isfield(c,'values') && (~isfield(c,'val') ||...
   isempty(c.val) || ~iscell(c.val) ||...
  ~strcmp(gettag(c.val{1}),node{1}))
    for i=1:length(c.values)
        if strcmp(gettag(c.values{i}),node{1})
            c.val = {c.values{i}};
            c.val{1} = make_nodes(c.val{1},node);
        end;
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,sts,varargout] = start_node1(c,node,fun,varargin)
varargout = {};
if nargin<4, varargin = {}; end;
if isfield(c,'tag') && ~strcmp(gettag(c),node{1})
    varargout = cell(1,nargout-2);
    sts       = 0;
    return;
end;
if isfield(c,'tag'),
    node = {node{2:end}};
end;
if isempty(node)
    [c,varargout{1:nargout-2}] = feval(fun,c,varargin{:});
    sts = 1;
    return;
end;
if isfield(c,'val'),
    for i=1:length(c.val)
        [c.val{i},sts,varargout{1:nargout-2}] = start_node1(c.val{i},node,fun,varargin{:});
        if sts, return; end;
    end;
end;
error('No such node'); % Should never execute
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,n,varargout] = run_in_current_node1(c,n,fun,modify,varargin)
% Satellite for run_in_current_node

varargout = {};
if nargin<5, varargin = {}; end;
if isfield(c,'hidden'), return; end;

n         = n-1;
if n<0, return; end;
if n==0,
    if ~isfield(c,'datacheck'),
        show_msg('');
    else
        show_msg(c.datacheck);
    end;
    [c,varargout{:}] = feval(fun,c,varargin{:});
end;

if isfield(c,'expanded') && ~isempty(c.val)
    if c.expanded,
        val   = c.val;
        c.val = {};
        for i=1:length(val)
            [tmp,n,varargout{:}] = run_in_current_node1(val{i},n,fun,modify,varargin{:});
            if ~isempty(tmp)
                if iscell(tmp)
                    c.val = {c.val{:}, tmp{:}};
                else
                    c.val = {c.val{:}, tmp};
                end;
            end;
        end;
        if modify && isfield(c,'check'),
            if all_set(c),
                [unused,val] = harvest(c);
                c.datacheck  = feval(c.check,val);
            end;
        end;
        if isfield(c,'datacheck') && ~isempty(c.datacheck),
            show_msg(c.datacheck);
        end;
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = click_batch_box_fun(c,unused)
% Mostly set up the options box, but also update the help and value boxes

if ~isempty(opts_box)
    dat  = {};
    str  = {};

    switch c.type
    case {'const'}
        % do nothing

    case {'files'}
        str = {str{:}, 'Specify Files'};
        dat = {dat{:}, struct('fun',@file_select,'args',{{}},'redraw',1)};

    case {'menu'}
        str = {str{:}, 'Specify Menu Item'};
        dat = {dat{:}, struct('fun',@menu_entry,'args',{{}},'redraw',0)};

    case {'entry'}
        str = {str{:}, 'Specify Text'};
        dat = {dat{:}, struct('fun',@text_entry,'args',{{}},'redraw',1)};

    case {'branch','choice','repeat'}
        if strcmp(c.type,'repeat')
            for i=1:length(c.values)
                if ~isfield(c.values{i},'hidden'),
                    str = {str{:},['New "' c.values{i}.name '"']};
                    dat = {dat{:}, struct('fun',@series,'args',{{c.values{i}}},'redraw',1)};
                end;
            end;

        elseif strcmp(c.type,'choice')
            for i=1:length(c.values)
                if ~isfield(c.values{i},'hidden'),
                    str = {str{:},['Choose "' c.values{i}.name '"']};
                    dat = {dat{:}, struct('fun',@choose,'args',{{c.values{i}}},'redraw',1)};
                end;
            end;
        end;
    end;

    if isfield(c,'removable')
        str = {str{:},['Remove Item "' c.name '"'],['Replicate Item "' c.name '"']};
        dat = {dat{:}, struct('fun',@remove,'args',{{}},'redraw',1),...
                       struct('fun',@replicate,'args',{{}},'redraw',1)};
    end;

    if ~same(get(opts_box,'String')',str)
        val = 1;
    else
        val = get(opts_box,'Value');
    end;
    val = max(min(length(str),val),min(1,length(str)));
    set(opts_box,'String',str,'UserData',dat,'Callback',@click_options_box,'Value',val);
end;

show_value(c);

% Update help
txt = '';
if isfield(c,'help'), txt = c.help; end;
help_box = findobj(0,'tag','help_box');
if ~isempty(help_box)
    set(help_box,'String','                    ');
    workaround(help_box);
    ext = get(help_box,'Extent');
    pos = get(help_box,'position');
    pw   = floor(pos(3)/ext(3)*21-4);
    set(help_box,'String',spm_justify(pw,txt));
    workaround(help_box);
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = expand_contract(c,unused)
% Expand/contract a node (for visualisation)

if isfield(c,'expanded') && ~isempty(c.val)
    if c.expanded
        c.expanded = false;
    else
        c.expanded = true;
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = menu_entry(c,unused)
% Setup opts box for menu entry

if isfield(c,'labels') && isfield(c,'values')
    str = {};
    dat = {};
    if ~isempty(c.val)
        val = c.val{1};
    else
        val = '<UNDEFINED>';
    end;
    dv = 1;
    for i=1:length(c.values)
        if ~(ischar(val) && strcmp(val,'<UNDEFINED>')) && same(c.values{i},val)
            str{i} = ['* ' c.labels{i}];
            dv = i;
        else
            str{i} = ['  ' c.labels{i}];
        end;
        dat{i} = struct('fun',@menuval,'args',{{c.values{i}}},'redraw',1);
    end;
    set(opts_box,'String',str,'UserData',dat,'Callback',@click_options_box,'Value',dv);
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c] = file_select(c)
% Select files
try
    set(batch_box,'Enable', 'inactive');
    set(opts_box, 'Enable', 'inactive');

    addvfiles(c.id);
    if ~isempty(c.val),
        sel = c.val{1};
    else
        sel = '';
    end;
    [s,ok] = spm_select(c.num,c.filter,c.name,sel);
    if ok, c.val{1} = cellstr(s); end;
catch
end;
spm_select('clearvfiles');
set(batch_box,'Enable', 'on');
set(opts_box, 'Enable', 'on');
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = text_entry(c)
% Create a string box and prompt, while hiding the opts box

opts_box = findobj(gcf,'tag','opts_box');
pos = get(opts_box,'Position');
set(opts_box,'Visible','off');
fs  = get(opts_box,'FontSize');
fa  = get(opts_box,'FontAngle');
fw  = get(opts_box,'FontWeight');
clf = get(opts_box,'ForegroundColor');
clb = get(opts_box,'BackgroundColor');

n = [];
m = [];
if isfield(c,'num'),    n = c.num;    end;
if isfield(c,'extras'), m = c.extras; end;

newstring = uicontrol(gcf,...
    'Style','edit',...
    'Units','normalized',...
    'Position',[pos(1:3) pos(4)/2],...
    'Tag','string_box',...
    'HorizontalAlignment','left',...
    'FontSize',fs,...
    'FontAngle',fa,...
    'FontWeight',fw,...
    'ForegroundColor',clf,...
    'BackgroundColor',clb,...
    'Callback',@accept_string);
cntxtmnu(newstring);

strM='';
switch lower(c.strtype)
case 's', TTstr='enter string';
case 'e', TTstr='enter expression - evaluated';
case 'n', TTstr='enter expression - natural number(s)';
        if ~isempty(m), strM=sprintf(' (in [1,%d])',m); TTstr=[TTstr,strM]; end
case 'w', TTstr='enter expression - whole number(s)';
        if ~isempty(m), strM=sprintf(' (in [0,%d])',m); TTstr=[TTstr,strM]; end
case 'i', TTstr='enter expression - integer(s)';
case 'r', TTstr='enter expression - real number(s)';
        if ~isempty(m), TTstr=[TTstr,sprintf(' in [%g,%g]',min(m),max(m))]; end
case 'c', TTstr='enter indicator vector e.g. 0101...  or abab...';
        if ~isempty(m) && isfinite(m), strM=sprintf(' (%d)',m); end
otherwise, TTstr='enter expression';
end

if isempty(n), n=NaN; end
n=n(:); if length(n)==1, n=[n,1]; end; dn=length(n);
if any(isnan(n)) || (prod(n)==1 && dn<=2) || (dn==2 && min(n)==1 && isinf(max(n)))
    str = '';
   lstr = '';
elseif dn==2 && min(n)==1
    str = sprintf('[%d]',max(n));
   lstr = [str,'-vector: '];
elseif dn==2 && sum(isinf(n))==1
    str = sprintf('[%d]',min(n));
   lstr = [str,'-vector(s): '];
else
    str='';
    for i = 1:dn,
        if isfinite(n(i)),
            str = sprintf('%s,%d',str,n(i));
        else
            str = sprintf('%s,*',str);
        end
    end
     str = ['[',str(2:end),']'];
    lstr = [str,'-matrix: '];
end
strN = sprintf('%s',lstr);
col  = get(gcf,'Color');
uicontrol(gcf,'Style','text',...
    'Units','normalized',...
    'BackgroundColor',col,...
    'String',[strN,strM,TTstr],...
    'Tag','string_prompt',...
    'HorizontalAlignment','Left',...
    'FontSize',fs,...
    'FontAngle',fa,...
    'FontWeight',fw,...
    'Position',[pos(1:3)+[0 pos(4)/2 0] pos(4)/2]);

if isfield(c,'val') && ~isempty(c.val)
    val = c.val{1};
    if ischar(val)
        tmp = val';
        tmp = tmp(:)';
        set(newstring,'String',tmp);
    elseif isnumeric(val)
        if ndims(val)>2,
            set(newstring,'String',['reshape([', num2str(val(:)'), '],[' num2str(size(val)) '])']);
        else
            if size(val,1)==1,
                set(newstring,'String',num2str(val(:)'));
            elseif size(val,2)==1,
                set(newstring,'String',['[' num2str(val(:)') ']''']);
            else
                str = '';
                for i=1:size(val,1),
                    str = [str ' ' num2str(val(i,:)) ';'];
                end;
                set(newstring,'String',str);
            end;
        end;
    end;
end;

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function accept_string(varargin)
% Accept a typed in string?
run_in_current_node(@stringval,true,get(varargin{1},'String'));
update_ui;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function remove_string_box
% Delete the string box and prompt, making the opts box
% visible again
delete(findobj(0,'Tag','string_box'));
delete(findobj(0,'Tag','string_prompt'));
set(opts_box,'Visible','on');
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function show_value(c)
% Update the value box

valtxt = {'<UNDEFINED>'};
% col    = [1 0 0];

switch c.type
case {'menu'}
    if isfield(c,'val') && ~isempty(c.val)
        if isfield(c,'labels') && isfield(c,'values')
            val = c.val{1};
            for i=1:length(c.values)
                if same(c.values{i},val)
                    valtxt = c.labels{i};
                    % col    = [0 0 0];
                    break;
                end;
            end;
        end;
    end;

case {'const','entry'}
    if isfield(c,'val') && ~isempty(c.val)
        val = c.val{1};
        if isempty(val)
            valtxt = '<Empty>';
        else
            if isnumeric(val)
                sz = size(val);
                if length(sz)>2
                    valtxt = sprintf('%dx',sz);
                    valtxt = [valtxt(1:(end-1)) ' Numeric Array'];
                else
                    valtxt = cell(size(val,1),1);
                    for i=1:size(val,1)
                        valtxt{i} = sprintf('%14.8g ',val(i,:));
                    end;
                end;
            elseif(ischar(val))
                valtxt = val;
            else
                valtxt = 'Can not display';
            end;
        end;
        % col    = [0 0 0];
    end;

case {'files'}
    if isfield(c,'val') && ~isempty(c.val)
        if isempty(c.val{1}) || isempty(c.val{1}{1})
            valtxt = '<None>';
        else
            valtxt = c.val{1};
        end;
        % col    = [0 0 0];
    end;

case {'choice'}
    if isfield(c,'val') && length(c.val)==1
        valtxt = {'A choice, where',['"' c.val{1}.name '"'], 'is selected.'};
        % col    = [0.5 0.5 0.5];
    else
        valtxt = {'A choice, with', 'nothing selected.'};
        % col    = [0.5 0.5 0.5];
    end;

case {'repeat'}
    ln     = length(c.val); s = 's'; if ln==1, s = ''; end;
    valtxt = ['A series containing ' num2str(length(c.val)) ' item' s '.'];
    % col    = [0.5 0.5 0.5];

case {'branch'}
    ln     = length(c.val); s = 's'; if ln==1, s = ''; end;
    valtxt = {['A branch holding ' num2str(length(c.val)) ' item' s '.']};
    if isfield(c,'prog'),
        valtxt = {valtxt{:}, '', '   User specified values',...
                                 '   from this branch will be',...
                                 '   collected and passed to',...
                                 '   an executable function',...
                                 '   when the job is run.'};
    end;
    % col    = [0.5 0.5 0.5];

otherwise
    valtxt = 'What on earth is this???';
end;
val_box = findobj(0,'tag','val_box');
if ~isempty(val_box)
    % set(val_box,'String',valtxt,'ForegroundColor',col);
    set(val_box,'String',valtxt);
    workaround(val_box);
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = remove(c,varargin)
% Remove node c
if strcmp(questdlg(['Remove "' c.name '"?'],'Confirm','Yes','No','Yes'),'Yes'),
    c     = {};
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = replicate(c,varargin)
% Replicate node c
c     = {c,uniq_id(c)};
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = choose(c,val)
% Specify the value of c to be val
c.val{1}   = uniq_id(val);
c.expanded = true;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = menuval(c,val)
% Specify the value of c to be val
c.val{1} = val;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = stringval(c,val)
% Accept (or not) a typed in string

msg = 'SUCCESS: Values accepted';

switch c.strtype
case {'s'}
    c.val{1} = val;
    show_value(c);
    remove_string_box;

case {'s+'}
    msg = 'FAILURE: Cant do s+ yet';
    beep;
    remove_string_box;

otherwise
    n = Inf;
    if isfield(c,'num')
        n      = c.num;
    end;
    if isfield(c,'extras')
        [val,msg] = spm_eeval(val,c.strtype,n,c.extras);
    else
        [val,msg] = spm_eeval(val,c.strtype,n);
    end;

    if ischar(val)
        beep;
        msg = ['FAILURE: ' msg];
    else
        c.val{1} = val;
        show_value(c);
        remove_string_box;
        msg = ['SUCCESS: ' msg];
    end;
end;
show_msg(msg);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = series(c,new)
% Add a new repeat
c.expanded    = true;
new.removable = true;
if isfield(c,'val')
    c.val = {c.val{:},uniq_id(new)};
else
    c.val = {uniq_id(new)};
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function setdef(strin,val)
global defaults

if isempty(defaults), spm_defaults; end;
if ischar(val) && strcmp(val,'<UNDEFINED>'), return; end;
o  = find(strin=='.');
df = cell(1,length(o)+1);
prev = 1;
for i=1:length(o),
    df{i} = strin(prev:(o(i)-1));
    prev  = o(i)+1;
end;
df{end}  = strin(prev:end);
defaults = setdef_sub(defaults,df,val);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function d = setdef_sub(d,df,val)
if isempty(df),
    d = val;
else
    if ~isfield(d,df{1}),d.(df{1}) = []; end;
    d.(df{1}) = setdef_sub(d.(df{1}),df(2:end),val);
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function val = getdef(strin)
% Get value of one of the SPM defaults

global defaults

val = getdef_sub(defaults,strin);
if ischar(val) && strcmp(val,'<UNDEFINED>')
    val = {};
else
    val = {val};
end;
%fprintf('%s\n',strin);
%disp(val)
%disp('----');
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = getdef_sub(defs,field)
% Satellite function for getdef

o = find(field=='.');
if isempty(o)
    if isfield(defs,field)
        c = defs.(field);
    else
        c = '<UNDEFINED>';
    end;
    return;
end;
if isfield(defs,field(1:(o-1)))
    c = getdef_sub(defs.(field(1:(o-1))),field((o+1):end));
else
    c = '<UNDEFINED>';
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function t = same(a,b)
% Are two data structures identical?

% Innocent until proven guilty
t = true;

% Check the dimensions
if isempty(a) && isempty(b)
    t = true; return;
end;
sa = size(a);
sb = size(b);
if length(sa) ~= length(sb)
    t = false; return;
end;
if ~all(sa==sb)
    t = false; return;
end;

% Check the classes
ca = class(a);
if ~strcmp(ca,class(b)), t = false; return; end;

% Recurse through data structure
switch ca
case {'double','single','sparse','char','int8','uint8',...
      'int16','uint16','int32','uint32','logical'}
    msk = ((a==b) | (isnan(a)&isnan(b)));
    if ~all(msk(:)), t = false; return; end;

case {'struct'}
    fa = fieldnames(a);
    fb = fieldnames(b);
    if length(fa) ~= length(fb), t = false; return; end;
    for i=1:length(fa)
        if ~strcmp(fa{i},fb{i}), t = false; return; end;
        for j=1:length(a)
            if ~same(a(j).(fa{i}),b(j).(fb{i}))
                t = false; return;
            end;
        end;
    end;

case {'cell'}
    for j=1:length(a(:))
        if ~same(a{j},b{j}), t = false; return; end;
    end;

case {'function_handle'}
    if ~strcmp(func2str(a),func2str(b))
        t = false; return;
    end;

otherwise
    warning(['Unknown class "' ca '"']);
    t = false;
end;

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function t = batch_box
t = findobj(0,'tag','batch_box');
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function t = opts_box
t = findobj(0,'tag','opts_box');
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function save_job(varargin)
% Save a batch job
c = get(batch_box,'UserData'); spm('Pointer','Watch');
[unused,jobs] = harvest(c);
spm('Pointer');
%eval([tag '=val;']);
cll = {'*.mat','Matlab .mat file';'*.xml','XML file'};
[filename, pathname, FilterIndex] = uiputfile(cll,'Save job as');
if ischar(filename)
    [unused,unused,ext] = fileparts(filename);
    if isempty(ext)
        ext = cll{FilterIndex}(2:end);
        filename = [filename ext];
    end
    if strcmp(ext,'.xml')
        spm('Pointer','Watch');
        savexml(fullfile(pathname,filename),'jobs');
        spm('Pointer');
    elseif strcmp(ext,'.mat')
        if str2double(version('-release'))>=14,
            save(fullfile(pathname,filename),'-V6','jobs');
        else
            save(fullfile(pathname,filename),'jobs');
        end;
    else
        questdlg(['Unknown extension (' ext ')'],'Nothing saved','OK','OK');
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function load_job(varargin)
% Load a batch job

cll = {'*.mat','Matlab .mat file';'*.xml','XML file'};
[filename, pathname] = uigetfile(cll,'Load job file');
if ischar(filename)
    [unused,nam,ext] = fileparts(filename);
    if strcmp(ext,'.xml')
        spm('Pointer','Watch');
        try
            loadxml(fullfile(pathname,filename),'jobs');
        catch
            questdlg('LoadXML failed',filename,'OK','OK');
            return;
        end;
        spm('Pointer');
    elseif strcmp(ext,'.mat')
        try
            load(fullfile(pathname,filename),'jobs');
        catch
            questdlg('Load failed',filename,'OK','OK');
        end;
    else
        questdlg(['Unknown extension (' ext ')'],'Nothing loaded','OK','OK');
    end;
    if exist('jobs','var')
        spm('Pointer','Watch');
        initialise(jobs);
        spm('Pointer');
    else
        questdlg(['No jobs (' nam ext ')'],'No jobs','OK','OK');
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function run_struct(varargin)
% Get data structure from handle, and run it

c = get(batch_box,'UserData');
run_struct1(c);
disp('--------------------------');
disp('Done.');
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function run_struct1(c)
% Run a batch job from a data structure

if isfield(c,'prog')
    prog = c.prog;
    [unused,val] = harvest(c);
    disp('--------------------------');
    disp(['Running "' c.name '"']);
    [Finter,unused,CmdLine] = spm('FnUIsetup',c.name);
    spm('Pointer','Watch');
    spm('FigName',[c.name ': running'],Finter,CmdLine);

    if 0
        try
            feval(prog,val);
            spm('FigName',[c.name ': done'],Finter,CmdLine);
        catch
            disp(['An error occurred when running "' c.name '"']);
            disp( '--------------------------------');
            disp(lasterr);
            disp( '--------------------------------');
            spm('FigName',[c.name ': failed'],Finter,CmdLine);
            errordlg({['An error occurred when running "' c.name '"'],lasterr},'SPM Jobs');
        end;
        spm('Pointer');
    else
        feval(prog,val);
        spm('FigName',[c.name ': done'],Finter,CmdLine);
        spm('Pointer');
    end;

else
    if isfield(c,'val')
        for i=1:length(c.val)
            run_struct1(c.val{i});
        end;
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function harvest_def(c)
switch c.type,
case{'const','menu','files','entry'},
    if isfield(c,'def') && numel(c.val)==1,
        setdef(c.def,c.val{1});
    end;

otherwise
    for i=1:length(c.val),
        harvest_def(c.val{i});
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [tag,val,typ] = harvest(c)
% Take a data structure, and extract what is needed to save it
% as a batch job

tag = 'unknown';
val = [];
typ = c.type;
switch(typ)
case {'const','menu','files','entry'}
    tag = gettag(c);
    if ~isempty(c.val)
        val = c.val{1};
    else
        val = '<UNDEFINED>';
    end;

case {'branch'}
    tag = gettag(c);
    if isfield(c,'val')
        val = [];
        for i=1:length(c.val)
            [tag1,val1] = harvest(c.val{i});
            val.(tag1)  = val1;
        end;
    end;

case {'repeat'}
    tag = gettag(c);
    val = {};
    if isfield(c,'val')
        for i=1:length(c.val),
            [tag1,val1,typ1] = harvest(c.val{i});
            if length(c.values)==1 && strcmp(typ1,'branch'),
                if i==1
                    val    = val1;
                else
                    val(i) = val1;
                end;
            else
                if length(c.values)>1,
                    if iscell(val1)
                        val1 = struct(tag1,{val1});
                    else
                        val1 = struct(tag1,val1);
                    end;
                end;
                if i==1
                    val = {val1};
                else
                    val = {val{:}, val1};
                end;
            end;
        end;
    end;

case {'choice'}
    if isfield(c,'tag'), tag = gettag(c); end;
    [tag1,val1] = harvest(c.val{1});
    if iscell(val1)
        val = struct(tag1,{val1});
    else
        val = struct(tag1,val1);
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function tag = gettag(c)
% Get a tag field - possibly from one of the kids

if (strcmp(c.type,'repeat') || strcmp(c.type,'choice')) && numel(c.values)>0
    tag = gettag(c.values{1});
    for i=2:length(c.values)
        if ~strcmp(tag,gettag(c.values{i}))
            tag = c.tag;
            return;
        end;
    end;
else
    tag = c.tag;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c0 = cntxtmnu(ob)
c0 = uicontextmenu('Parent',get(ob,'Parent'));
set(ob,'uicontextmenu',c0);
c1 = uimenu('Label','Font', 'Parent',c0);
     uimenu('Label','Plain',      'Parent',c1,'Callback','set(gco,''FontWeight'',''normal'',''FontAngle'',''normal'');');
     uimenu('Label','Bold',       'Parent',c1,'Callback','set(gco,''FontWeight'',''bold'',  ''FontAngle'',''normal'');');
     uimenu('Label','Italic',     'Parent',c1,'Callback','set(gco,''FontWeight'',''normal'',''FontAngle'',''italic'');');
     uimenu('Label','Bold-Italic','Parent',c1,'Callback','set(gco,''FontWeight'',''bold'',  ''FontAngle'',''italic'');');
c1 = uimenu('Label','Fontsize','Parent',c0);
fs = [8 9 10 12 14 16 18]; % [20 24 28 32 36 44 48 54 60 66 72 80 88 96];
for i=fs,
    uimenu('Label',sprintf('%-3d',i),'Parent',c1,'Callback',@fszoom,'UserData',i);
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function workaround(t)
set(t,'Value',[], 'Enable', 'on', 'Max',2, 'Min',0,'ListBoxTop',1);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function addvfiles(id,c)
if nargin<2,
    c  = get(batch_box,'UserData');
end;
vf =addvfiles1(c,id,{});
spm_select('clearvfiles');
spm_select('addvfiles',vf);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [vf,sts]=addvfiles1(c,id,vf)
sts = 0;
if ~isstruct(c) || ~isfield(c,'type'),
    return;
end;

if isfield(c,'vfiles'),
    if ~find_id(c,id),
        [c,unused,ok] = get_strings1(c,0);
        if ok,
            [unused,job] = harvest(c);
            vf1          = feval(c.vfiles,job);
            vf           = {vf{:}, vf1{:}};
        end;
    else
        sts = 1;
    end;
    return;
end;

switch c.type,
case {'repeat','choice','branch'},
    for i=1:length(c.val),
        [vf,sts]=addvfiles1(c.val{i},id,vf);
        if sts, return; end;
    end;
end;

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function ok = find_id(c,id)
ok = 0;
if ~isstruct(c) || ~isfield(c,'type'),
    return;
end;
if c.id==id,
    ok = 1;
    return;
end;
if isfield(c,'val'),
    for i=1:length(c.val),
        ok = find_id(c.val{i},id);
        if ok, return; end;
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = uniq_id(c)
if ~isstruct(c) || ~isfield(c,'type'),
    return;
end;
c.id = rand(1);
if isfield(c,'val'),
    for i=1:length(c.val),
        c.val{i} = uniq_id(c.val{i});
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function show_msg(txt)
lb  = findobj('tag','msg_box');
if isempty(txt),
    set(lb,'String',{});
else
    msg = get(lb,'String');
    if iscell(txt),
        msg = {msg{:} txt{:}};
    else
        msg = {msg{:} txt};
    end;
    set(lb,'String',msg);
end;
drawnow;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%function beep
%fprintf('%c',7);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function pulldown
% Create a pulldown for individual jobs
c  = initialise_struct;
fg = spm_figure('findwin','Graphics');
if isempty(fg), return; end;
set(0,'ShowHiddenHandles','on');
delete(findobj(fg,'tag','jobs'));
set(0,'ShowHiddenHandles','off');
f0 = uimenu(fg,'Label','TASKS','HandleVisibility','off','tag','jobs');
pulldown1(f0,c,c.tag);
uimenu(f0,'Label','Batch','CallBack',@interactive,'Separator','on');
uimenu(f0,'Label','Defaults','CallBack',@defaults_edit,'Separator','off');
f1 = uimenu(f0,'Label','Sequential');
pulldown2(f1,c,c.tag);

if 0, % Currently unused
    f1 = uimenu(f0,'Label','Modality');
    modalities = {'FMRI','PET','EEG'};
    for i=1:length(modalities)
        tmp = modalities{i};
        if strcmp(tmp,getdef('modality')),
            tmp = ['*' tmp];
        else
            tmp = [' ' tmp];
        end;
        uimenu(f1,'Label',tmp,'CallBack',@chmod);
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function pulldown1(f0,c0,tag0)
if ~isfield(c0,'values'), return; end;
for i=1:length(c0.values),
    c1 = c0.values{i};
    if isstruct(c1) && ~isfield(c1,'hidden'),
        tag1 = tag0;
        if isfield(c1,'tag'),
            tag1 = [tag1 '.' c1.tag];
        end;
        if isfield(c1,'prog'),
            uimenu(f0,'Label',c1.name,'CallBack',@interactive,'UserData',{'',tag1});
        else
            f1 = uimenu(f0,'Label',c1.name);
            pulldown1(f1,c1,tag1);
        end;
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function pulldown2(f0,c0,tag0)
if ~isfield(c0,'values'), return; end;
for i=1:length(c0.values),
    c1 = c0.values{i};
    if isstruct(c1) && ~isfield(c1,'hidden'),
        tag1 = tag0;
        if isfield(c1,'tag'),
            tag1 = [tag1 '.' c1.tag];
        end;
        if isfield(c1,'prog'),
             if findcheck(c1),
                 uimenu(f0,'Label',c1.name,'Enable','off');
             else
                 uimenu(f0,'Label',c1.name,'CallBack',@run_serial,'UserData',{'',tag1});
             end;
        else
            f1 = uimenu(f0,'Label',c1.name);
            pulldown2(f1,c1,tag1);
        end;
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function hascheck = findcheck(c)
hascheck = false;
if ~isstruct(c) || ~isfield(c,'type'), return; end;
if isfield(c,'check'),
    hascheck = true;
    return;
end;
if isfield(c,'values'),
    for i=1:numel(c.values),
        hascheck = findcheck(c.values{i});
        if hascheck, return; end;
    end;
end;
if isfield(c,'val'),
    for i=1:numel(c.val),
        hascheck = findcheck(c.val{i});
        if hascheck, return; end;
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function chmod(mod,varargin)
global defaults
if isempty(defaults), spm_defaults; end;
if ischar(mod),
    %if strcmpi(defaults.modality,mod),
    %    spm('ChMod',mod);
    %end;
    defaults.modality = mod;
else
    tmp = get(mod,'Label');
    defaults.modality = tmp(2:end);
end;
pulldown;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function fszoom(varargin)
fs = sscanf(get(varargin{1},'Label'),'%d');
set(gco,'FontSize',fs);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function run_serial(varargin)
ud = get(varargin{1},'UserData');
if iscell(ud)
    serial(ud{:});
else
    serial;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function serial(job,node)

if nargin<2, node = 'jobs'; end;

fg = spm_figure('FindWin','Interactive');
if isempty(fg), fg = spm('CreateIntWin'); end;
delete(findobj(fg,'Parent',fg));
t=uicontrol(fg,...
    'Style','listbox',...
    'Units','normalized',...
    'Position',[0.02 0.02 0.96 0.62],...
    'Tag','help_box2',...
    'FontName','fixedwidth',...
    'FontSize',12,...
    'BackgroundColor',[1 1 1]);
set(t,'Value',[], 'Enable', 'inactive', 'Max',2, 'Min',0);
workaround(t);
cntxtmnu(t);
spm('Pointer');
drawnow;

if nargin>0,
    c = initialise_struct(job);
else
    c = initialise_struct;
end;

c       = start_node(c,node);
c       = start_node(c,@run_ui,{});
spm_input('!DeleteInputObj');
delete(findobj(fg,'Parent',fg));
[unused,jobs] = harvest(c);
%savexml('job_last.xml','jobs');
run_job(jobs);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = run_ui(c,varargin)

nnod=1;
while(1),
    nod = nnod;
    [ci,unused,hlp] = get_node(c,nod);
    if isempty(ci), break; end;

    help_box = findobj(0,'tag','help_box2');
    if ~isempty(help_box),
        set(help_box,'String','                    ');
        workaround(help_box);
        ext = get(help_box,'Extent');
        pos = get(help_box,'position');
        pw   = floor(pos(3)/ext(3)*21-4);
        set(help_box,'String',spm_justify(pw,hlp));
        workaround(help_box);

        if isfield(c,'prog'),
            try
                set(help_box,'HandleVisibility','off');
                [Finter,unused,CmdLine] = spm('FnUIsetup',c.name);
                spm('FigName',[c.name ': setup'],Finter,CmdLine);
            catch
            end;
            set(help_box,'HandleVisibility','on');
        end;
    end;

    pos = 1;

    switch ci.type,
    case {'const','files','menu','entry'}
        nnod = nod + 1;

        vl = {'<UNDEFINED>'};
        if isfield(ci,'def'), vl = getdef(ci.def); end;
        if numel(vl)~=0 && (~ischar(vl{1}) || ~strcmp(vl{1},'<UNDEFINED>')),
            getit = 0;
            if ~isfield(ci,'val') || ~iscell(ci.val) || isempty(ci.val),
                ci.val = vl;
            end;
        else
            getit = 1;
        end;

        switch ci.type,
        case {'const'}
        case {'files'}
            num      = ci.num;

            if getit,
                if ~isempty(ci.val),
                    sel = ci.val{1};
                else
                    sel = '';
                end;
                addvfiles(ci.id,c);
                [ci.val{1},ok] = spm_select(num,ci.filter,ci.name,sel);
                if ~ok,
                    error('File Selector was deleted.');
                end;
                spm_select('clearvfiles');
                ci.val{1} = cellstr(ci.val{1});
            end;

        case {'menu'}
            dv = [];
            if getit,
                if isfield(ci,'val') && ~isempty(ci.val),
                    for i=1:length(ci.values)
                        if same(ci.values{i},ci.val{1})
                            dv = i;
                        end;
                    end;
                end;
                lab = ci.labels{1};
                for i=2:length(ci.values),
                    lab = [lab '|' ci.labels{i}];
                end;
                if isempty(dv),
                    ind   = spm_input(ci.name,pos,'m',lab,1:length(ci.values));
                else
                    ind   = spm_input(ci.name,pos,'m',lab,1:length(ci.values),dv);
                end;
                ci.val = {ci.values{ind}};
            end;

        case {'entry'}
            n1  = Inf;
            if isfield(ci,'num'), n1 = ci.num; end;
            if getit,
                val = '';
                if isfield(ci,'val') && ~isempty(ci.val) && ~strcmp(ci.val{1},'<UNDEFINED>'),
                    val = ci.val{1};
                end;
                if isfield(ci,'extras')
                    val = spm_input(ci.name,pos,ci.strtype,val,n1,ci.extras);
                else
                    val = spm_input(ci.name,pos,ci.strtype,val,n1);
                end;
                ci.val{1} = val;
            end;
        end;

    case {'repeat'},
        lab = 'Done';
        for i=1:length(ci.values)
            lab = [lab '|New "' ci.values{i}.name '"'];
        end;
        tmp = spm_input(ci.name,pos,'m',lab,0:length(ci.values));
        if tmp,
            ci.val = {ci.val{:}, uniq_id(ci.values{tmp})};
            nnod = nod;
        else
            nnod = nod + 1;
        end;

    case {'choice'}
        nnod = nod + 1;
        lab  = ci.values{1}.name;
        for i=2:length(ci.values)
            lab = [lab '|' ci.values{i}.name];
        end;
        tmp = spm_input(ci.name,pos,'m',lab,1:length(ci.values));
        ci.val = {uniq_id(ci.values{tmp})};

    end;
    c = set_node(c,nod,ci);
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [ci,n,hlp] = get_node(c,n)
ci  = [];
hlp = '';
switch c.type,
case {'const','files','menu','entry','choice'}
    n = n-1;
    if n==0,
        ci  = c;
        hlp = {['* ' upper(c.name)],c.help{:},'',''};
        return;
    end;

case 'branch',
    for i=1:numel(c.val),
        [ci,n,hlp] = get_node(c.val{i},n);
        if ~isempty(ci),
            hlp = {hlp{:},repmat('=',1,20),'',['* ' upper(c.name)],c.help{:},'',''};
            return;
        end;
    end;

case 'repeat',
    for i=1:numel(c.val),
        [ci,n,hlp] = get_node(c.val{i},n);
        if ~isempty(ci),
            hlp = {hlp{:},repmat('=',1,20),'',['* ' upper(c.name)],c.help{:},'',''};
            return;
        end;
    end;
    n = n-1;
    if n==0,
        ci  = c;
        hlp = {['* ' upper(c.name)],c.help{:},'',''};
        return;
    end;

end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,n] = set_node(c,n,ci)
switch c.type,
case {'const','files','menu','entry','choice'}
    n = n-1;
    if n==0,
        c = ci;
        return;
    end;

case 'branch',
    for i=1:numel(c.val),
        [c.val{i},n] = set_node(c.val{i},n,ci);
        if n<0,return; end;
    end;

case 'repeat',
    for i=1:numel(c.val),
        [c.val{i},n] = set_node(c.val{i},n,ci);
        if n<0,return; end;
    end;
    n = n-1;
    if n==0,
        c = ci;
        return;
    end;

end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = initialise_struct(job)
% load the config file, possibly adding a job
% to it, and generally tidy it up.  The batch box
% is updated to contain the current structure.

persistent c0
if isempty(c0),
    c0 = spm_config;
    c0 = tidy_struct(c0);
end;
c = insert_defs(c0);
if nargin==1 && ischar(job) && strcmp(job,'defaults'),
    c      = defsub(c,{});
    c.name = 'SPM Defaults';
else
    c      = hide_null_jobs(c);
    if nargin>0 && ~isempty(job),
        job = fromfile(job);
        c   = job_to_struct(c,job,'jobs');
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,defused] = defsub(c,defused)
if nargin<2, defused = {}; end;
if isfield(c,'prog'), c = rmfield(c,'prog'); end;
switch c.type,
case {'const'}
    c = [];

case {'menu','entry','files'}
    if ~isfield(c,'def') || any(strcmp(c.def,defused)),
        c = [];
    else
        defused = {defused{:},c.def};
    end;

case {'branch'}
    msk = true(1,length(c.val));
    for i=1:length(c.val),
        [c.val{i},defused] = defsub(c.val{i},defused);
        msk(i)   = ~isempty(c.val{i});
    end;
    c.val = c.val(msk);
    if isempty(c.val), c = []; end;

case {'choice','repeat'}
    c.type = 'branch';
    c.val  = c.values;
    c      = rmfield(c,'values');
    [c,defused] = defsub(c,defused);

end;
if isfield(c,'vfiles'), c = rmfield(c,'vfiles'); end;
if isfield(c,'check'),  c = rmfield(c,'check');  end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = insert_defs(c)
% Recursively descend through the tree structure, 
% and assigning default values.
if ~isstruct(c) || ~isfield(c,'type'),
    return;
end;
switch c.type
case {'menu','entry','files'}
    if isfield(c,'def')
        c.val = getdef(c.def);
        if strcmp(c.type,'files') && ~isempty(c.val)
            c.val = {cellstr(c.val{1})};
        end;
    end;

case {'repeat','choice'},
    if isfield(c,'values')
        for i=1:numel(c.values)
            c.values{i} = insert_defs(c.values{i});
        end;
    end;
end;
if isfield(c,'val')
    for i=1:numel(c.val)
        c.val{i} = insert_defs(c.val{i});
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = tidy_struct(c)
% Recursively descend through the tree structure, cleaning up
% fields that may be missing and adding an 'expanded' field
% where necessary.

if ~isstruct(c) || ~isfield(c,'name') || ~isfield(c,'type')
    return;
end;
c.id = rand(1);

if ~isfield(c,'help'), c.help = {}; end;
if ischar(c.help), c.help = {c.help}; end;

switch c.type
case {'const'}
    if ~isfield(c,'tag')
        disp(c); warning(['No tag field for "' c.name '"']);
        c.tag = 'unknown';
    end;
    if ~isfield(c,'val')
        disp(c); warning(['No val field for "' c.name '"']);
        c.val = {'<UNDEFINED>'};
    end;

case {'menu'}
    if ~isfield(c,'tag')
        disp(c); warning(['No tag field for "' c.name '"']);
        c.tag = 'unknown';
    end;

    if ~isfield(c,'labels') || ~isfield(c,'values')
        disp(c); warning(['No labels and values field for "' c.name '"']);
        c.labels = {};
        c.values = {};
    end;
    if length(c.labels) ~= length(c.values)
        disp(c); warning(['Labels and values fields incompatible for "' c.name '"']);
        c.labels = {};
        c.values = {};
    end;
    if ~isfield(c,'help'), c.help = {'Option selection by menu'}; end;

case {'entry'}
    if ~isfield(c,'tag')
        disp(c); warning(['No tag field for "' c.name '"']);
        c.tag = 'unknown';
    end;
    if ~isfield(c,'strtype')
        disp(c); warning(['No strtype field for "' c.name '"']);
        c.strtype = 'e';
    end;
    if ~isfield(c,'num')
        disp(c); warning(['No num field for "' c.name '"']);
        c.num = [1 1];
    end;
    if length(c.num)~=2
        disp(c); warning(['Num field for "' c.name '" is wrong length']);
        c.num = [Inf 1];
    end;
    if ~isfield(c,'help'), c.help = {'Option selection by text entry'}; end;

case {'files'}
    if ~isfield(c,'tag')
        disp(c); warning(['No tag field for "' c.name '"']);
        c.tag = 'unknown';
    end;
    if ~isfield(c,'filter')
        disp(c); warning(['No filter field for "' c.name '"']);
        c.filter = '*';
    end;
    if ~isfield(c,'num')
        disp(c); warning(['No num field for "' c.name '"']);
        c.num = Inf;
    end;
    if length(c.num)~=1 && length(c.num)~=2
        disp(c); warning(['Num field for "' c.name '" is wrong length']);
        c.num = Inf;
    end;
    if ~isfield(c,'help'), c.help = {'File selection'}; end;

case {'branch'}
    if ~isfield(c,'tag')
        disp(c); warning(['No tag field for "' c.name '"']);
        c.tag = 'unknown';
    end;
    if ~isfield(c,'val')
        disp(c); warning(['No val field for "' c.name '"']);
        c.val = {};
    end;

    c.expanded = false;
    if ~isfield(c,'help'), c.help = {'Branch structure'}; end;

case {'choice'}
    if ~isfield(c,'tag')
        disp(c); warning(['No tag field for "' c.name '"']);
        c.tag = 'unknown';
    end;
    if ~isfield(c,'values') || ~iscell(c.values)
        disp(c); error(['Bad values for "' c.name '"']);
    end;
    for i=1:length(c.values)
        c.values{i} = tidy_struct(c.values{i});
    end;
    if ~isfield(c,'val') || ~iscell(c.val) || length(c.val) ~= 1
        c.val = {c.values{1}};
    end;
    c.expanded = true;
    if ~isfield(c,'help'), c.help = {'Choice structure'}; end;

case {'repeat'}
    if ~isfield(c,'values') || ~iscell(c.values)
        disp(c); error(['Bad values for "' c.name '"']);
    end;
    for i=1:length(c.values)
        c.values{i} = tidy_struct(c.values{i});
    end;
    if length(c.values)>1 && ~isfield(c,'tag')
        disp(c); warning(['No tag field for "' c.name '"']);
        c.tag = 'unknown';
    end;
    if length(c.values)==1 && isfield(c,'tag')
    %    disp(c); warning(['"' c.name '" has unused tag']);
        c = rmfield(c,'tag');
    end;
    c.expanded = true;
    if ~isfield(c,'help'), c.help = {'Repeated structure'}; end;

    if isfield(c,'num') && numel(c.num)==1,
        if finite(c.num),
            c.num = [c.num c.num];
        else
            c.num = [0 c.num];
        end;
    end;

end;

if ~isfield(c,'val'), c.val = {}; end;

%switch c.type
%case {'menu','entry','files'}
%    %if isempty(c.val)
%        if isfield(c,'def')
%            c.val = getdef(c.def);
%            if strcmp(c.type,'files') && ~isempty(c.val)
%                c.val = {cellstr(c.val{1})};
%            end;
%        end;
%    %end;
%end;

if isfield(c,'val')
    for i=1:length(c.val)
        c.val{i} = tidy_struct(c.val{i});
    end;
end;
if isfield(c,'values') && strcmp(c.type,'repeat')
    for i=1:length(c.values)
        c.values{i} = tidy_struct(c.values{i});
    end;
end;

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function job = fromfile(job)
if ischar(job)
    [pth,nam,ext] = fileparts(job);
    if strcmp(ext,'.xml')
        spm('Pointer','Watch');
        try
            loadxml(job,'jobs');
        catch
            questdlg('LoadXML failed',job,'OK','OK');
            return;
        end;
        spm('Pointer');
    elseif strcmp(ext,'.mat')
        try
            load(job,'jobs');
        catch
            questdlg('Load failed',job,'OK','OK');
            return;
        end;
    else
        questdlg(['Unknown extension (' ext ')'],'Nothing loaded','OK','OK');
        return;
    end;
    job = jobs;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = hide_null_jobs(c)
c = hide_null_jobs1(c);
c = hide_null_jobs2(c);
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,flg] = hide_null_jobs1(c)
if ~isstruct(c) || ~isfield(c,'type'),
    flg = true;
    return;
end;

if ~include(c),
    flg      = false;
    c.hidden = true;
    return;
end;

switch c.type,
case {'repeat','branch','choice'},
    msk1 = true;
    msk2 = true;
    if isfield(c,'val') && ~isempty(c.val),
        msk1 = true(1,numel(c.val));
        for i=1:length(c.val)
            [c.val{i},msk1(i)] = hide_null_jobs1(c.val{i});
        end;
    end;
    if isfield(c,'values') && ~isempty(c.values),
        msk2 = true(1,numel(c.values));
        for i=1:length(c.values)
            [c.values{i},msk2(i)] = hide_null_jobs1(c.values{i});
        end;
    end;
    flg = any(msk1) || any(msk2);
    if ~flg, c.hidden = true; end;
otherwise
    flg = true;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function ok = include(c)
% Check that the modality is OK
ok = true;
if isfield(c,'modality'),
    mod = getdef('modality');
    if ~isempty(mod),
        mod = mod{1};
        ok  = false;
        for i=1:length(c.modality),
            if strcmpi(c.modality{i},mod),
                ok = true;
                return;
            end;
        end;
   end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [c,flg] = hide_null_jobs2(c)
if ~isstruct(c) || ~isfield(c,'type')
    flg = 0;
    return;
end;
if isfield(c,'prog'),
    flg = 1;
    return;
end;

switch c.type,
case {'repeat','branch','choice'},
    flg = 0;
    msk1 = [];
    msk2 = [];
    if isfield(c,'val'),
        msk1 = ones(1,numel(c.val));
        for i=1:length(c.val)
            [c.val{i},msk1(i)] = hide_null_jobs2(c.val{i});
        end;
    end;
    if isfield(c,'values'),
        msk2 = ones(1,numel(c.values));
        for i=1:length(c.values)
            [c.values{i},msk2(i)] = hide_null_jobs2(c.values{i});
        end;
    end;
    if (sum(msk1) + sum(msk2))>0, flg = 1; end;
    if ~flg, c.hidden = 1; end;
otherwise
    flg = 0;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function c = job_to_struct(c,job,tag)
% Modify a structure based on a batch job
if isstruct(c) && isfield(c,'hidden'),
    return;
end;
switch c.type
case {'const','menu','files','entry'}
    if ~strcmp(gettag(c),tag), return; end;
    if ischar(job) && strcmp(job,'<UNDEFINED>')
        c.val = {};
    else
        c.val{1} = job;
    end;

case {'branch'}
    if ~strcmp(gettag(c),tag), return; end;
    if ~isstruct(job), return; end;

    tag = fieldnames(job);
    for i=1:length(tag)
        for j=1:length(c.val)
            if strcmp(gettag(c.val{j}),tag{i})
                c.val{j} = job_to_struct(c.val{j},job.(tag{i}),tag{i});
                break;
            end;
        end;
    end;

case {'choice'}
    if ~strcmp(gettag(c),tag), return; end;
    if ~isstruct(job), return; end;
    tag = fieldnames(job);
    if length(tag)>1, return; end;
    tag = tag{1};
    for j=1:length(c.values)
        if strcmp(gettag(c.values{j}),tag)
            c.val = {job_to_struct(c.values{j},job.(tag),tag)};
        end;
    end;

case {'repeat'}
    if ~strcmp(gettag(c),tag), return; end;
    if length(c.values)==1 && strcmp(c.values{1}.type,'branch')
        if ~isstruct(job), return; end;
        for i=1:length(job)
            if strcmp(gettag(c.values{1}),tag)
                c.val{i} = job_to_struct(c.values{1},job(i),tag);
                c.val{i}.removable = true;
            end;
        end;
    elseif length(c.values)>1
        if ~iscell(job), return; end;
        for i=1:length(job)
            tag = fieldnames(job{i});
            if length(tag)>1, return; end;
            tag = tag{1};
            for j=1:length(c.values)
                if strcmp(gettag(c.values{j}),tag)
                    c.val{i} = job_to_struct(c.values{j},job{i}.(tag),tag);
                    c.val{i}.removable = true;
                    break;
                end;
            end;
        end;
    else
        if ~iscell(job), return; end;
        for i=1:length(job)
            c.val{i} = job_to_struct(c.values{1},job{i},tag);
            c.val{i}.removable = true;
        end;
    end;

end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function doc = showdoc(str,wid)
if nargin<1, str = ''; end;
if nargin<2, wid = 60; end;

tmp = [0 find([str '.']=='.')];
node = {};
for i=1:length(tmp)-1,
    tmp1 = str((tmp(i)+1):(tmp(i+1)-1));
    if ~isempty(tmp1),
        node = {node{:},tmp1};
    end;
end;
if numel(node)>1 && strcmp(node{1},'jobs'),
    node = node(2:end);
end;

c = initialise_struct;
doc = showdoc1(node,c,wid);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function doc = showdoc1(node,c,wid)
doc   = {};
if isempty(node),
    doc = showdoc2(c,'',wid);
    return;
end;

if isfield(c,'values'),
    for i=1:numel(c.values),
        if isfield(c.values{i},'tag') && strcmp(node(1),c.values{i}.tag),
            doc = showdoc1(node(2:end),c.values{i},wid);
            return;
        end;
    end;
end;

if isfield(c,'val'),
    for i=1:numel(c.val),
        if isfield(c.val{i},'tag') && strcmp(node(1),c.val{i}.tag),
            doc = showdoc1(node(2:end),c.val{i},wid);
            return;
        end;
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function doc = showdoc2(c,lev,wid)
doc = {''};
if ~isempty(lev) && sum(lev=='.')==1,
        % doc = {doc{:},repmat('_',1,80),''};
end;

if isfield(c,'name'),
    str   = sprintf('%s %s', lev, c.name);
    %under = repmat('-',1,length(str));
    doc = {doc{:},str};
    % if isfield(c,'modality'),
    %     txt = 'Only for ';
    %     for i=1:numel(c.modality),
    %         txt = [txt ' ' c.modality{i}];
    %     end;
    %     doc = {doc{:},'',txt, ''};
    %end;
    if isfield(c,'help');
        hlp = spm_justify(wid,c.help);
        doc = {doc{:},hlp{:}};
    end;

    switch (c.type),
    case {'repeat'},
        if length(c.values)==1,
            doc = {doc{:}, '', sprintf('Repeat "%s", any number of times.',c.values{1}.name)};
        else
            doc = {doc{:}, '', 'Any of the following options can be chosen, any number of times'};
            i = 0;
            for ii=1:length(c.values),
                if isstruct(c.values{ii}) && isfield(c.values{ii},'name'),
                    i    = i+1;
                    doc = {doc{:}, sprintf('    %2d) %s', i,c.values{ii}.name)};
                end;
            end;
        end;
        doc = {doc{:},''};

    case {'choice'},
        doc = {doc{:}, '', 'Any one of these options can be selected:'};
        i = 0;
        for ii=1:length(c.values),
            if isstruct(c.values{ii}) && isfield(c.values{ii},'name'),
                i    = i+1;
                doc = {doc{:}, sprintf('    %2d) %s', i,c.values{ii}.name)};
            end;
        end;
        doc = {doc{:},''};

    case {'branch'},
        doc = {doc{:}, '', sprintf('This item contains %d fields:', length(c.val))};
        i = 0;
        for ii=1:length(c.val),
            if isstruct(c.val{ii}) && isfield(c.val{ii},'name'),
                i    = i+1;
                doc = {doc{:}, sprintf('    %2d) %s', i,c.val{ii}.name)};
            end;
        end;
        doc = {doc{:},''};

    case {'menu'},
        doc = {doc{:}, '', 'One of these values is chosen:'};
        for k=1:length(c.labels),
            doc = {doc{:}, sprintf('    %2d) %s', k, c.labels{k})};
        end;
        doc = {doc{:},''};

    case {'files'},
        if length(c.num)==1 && isfinite(c.num(1)) && c.num(1)>=0,
            tmp = spm_justify(wid,sprintf('A "%s" file is selected by the user.',c.filter));
        else
            tmp = spm_justify(wid,sprintf('"%s" files are selected by the user.\n',c.filter));
        end;
        doc = {doc{:}, '', tmp{:}, ''};

    case {'entry'},
        switch c.strtype,
        case {'e'},
            d = 'Evaluated statements';
        case {'n'},
            d = 'Natural numbers';
        case {'r'},
            d = 'Real numbers';
        case {'w'},
            d = 'Whole numbers';
        otherwise,
            d = 'Values';
        end;
        tmp = spm_justify(wid,sprintf('%s are entered.',d));
        doc = {doc{:}, '', tmp{:}, ''};
    end;

    i = 0;
    doc = {doc{:},''};
    if isfield(c,'values'),
        for ii=1:length(c.values),
            if isstruct(c.values{ii}) && isfield(c.values{ii},'name'),
                i    = i+1;
                lev1 = sprintf('%s%d.', lev, i);
                doc1 = showdoc2(c.values{ii},lev1,wid);
                doc  = {doc{:}, doc1{:}};
            end;
        end;
    end;
    if isfield(c,'val'),
        for ii=1:length(c.val),
            if isstruct(c.val{ii}) && isfield(c.val{ii},'name'),
                i    = i+1;
                lev1 = sprintf('%s%d.', lev, i);
                doc1 = showdoc2(c.val{ii},lev1,wid);
                doc  = {doc{:}, doc1{:}};
            end;
        end;
    end;
    doc = {doc{:}, ''};
end;


