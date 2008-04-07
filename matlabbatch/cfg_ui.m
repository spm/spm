function varargout = cfg_ui(varargin)
% CFG_UI M-File for cfg_ui.fig
%      CFG_UI, by itself, creates a new CFG_UI or raises the existing
%      singleton*.
%
%      H = CFG_UI returns the handle to a new CFG_UI or the handle to
%      the existing singleton*.
%
%      CFG_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CFG_UI.M with the given input arguments.
%
%      CFG_UI('Property','Value',...) creates a new CFG_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cfg_ui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cfg_ui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_ui.m 1312 2008-04-07 07:47:40Z volkmar $

rev = '$Rev: 1312 $';

% edit the above text to modify the response to help cfg_ui

% Last Modified by GUIDE v2.5 25-Mar-2008 19:32:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cfg_ui_OpeningFcn, ...
                   'gui_OutputFcn',  @cfg_ui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%% Local functions
% Most callbacks just refer to one of these local_ functions.

% --------------------------------------------------------------------
function local_DelMod(hObject)
handles = guidata(hObject);
udmodlist = get(handles.modlist,'userdata');
val = get(handles.modlist,'value');
if ~isempty(udmodlist.cmod)
    cfg_util('delfromjob',udmodlist.cjob, udmodlist.id{val});
    local_showjob(hObject);
end;

% --------------------------------------------------------------------
function local_ReplMod(hObject)
handles = guidata(hObject);
udmodlist = get(handles.modlist,'userdata');
val = get(handles.modlist,'value');
if ~isempty(udmodlist.cmod)
    cfg_util('replicate',udmodlist.cjob, udmodlist.id{val});
    local_showjob(hObject);
end;

%% Input evaluation
% --------------------------------------------------------------------
function [val sts] = local_eval_valedit(varargin)
% for security reasons, separate string evaluation from valedit_Callback
% (evaluation might overwrite variables). Uses evalc to suppress any
% console output.
val = [];
sts = false;
try
    % try to evaluate str as rvalue
    val = evalin('base', varargin{1});
    sts = true;
catch
    everr = lasterror;
    if strcmp(everr.identifier, 'MATLAB:m_invalid_lhs_of_assignment')
        try
            evalin('base', varargin{1});
            val = evalin('base','val');
            % test if val variable exists
            if ~exist('val','var')
                error('cfg_ui:local_eval_valedit:noval','No variable ''val'' assigned.');
            end;
            sts = true;
        catch
            sts = false;
            val = [];
            everr = lasterror;
            msgbox(everr.message,'Evaluation error','modal');
        end;
    end;
end;

%% Dynamic Menu Creation
% --------------------------------------------------------------------
function local_setmenu(hObject)
% delete previous added menu, if any
handles = guidata(hObject);
prevmenu = findobj(handles.cfg_ui, 'Tag', 'AddedAppMenu');
if ~isempty(prevmenu)
    delete(prevmenu);
end;
% get strings and ids
[id,stop,val]=cfg_util('listcfgall',[],cfg_findspec({{'hidden',false}}),{'name','level'});
str = val{1};
lvl = cat(1,val{2}{:});
% strings and ids are in preorder - if stop is true, then we are at leaf
% level, if lvl(k) <= lvl(k-1) we are at siblings/parent level
% remember last menu at lvl for parents, start at lvl > 1 (top parent is
% figure menu)
lastmenulvl = zeros(1, max(lvl));
lastmenulvl(1) = handles.cfg_ui;
toplevelmenus = [];
toplevelids   = {};
% 1st entry is top of matlabbatch config tree, applications start at 2nd entry
for k = 2:numel(lvl)
    if stop(k)
        label = sprintf('New: %s', str{k});
        udata = id{k};
        cback = @local_addtojob;
    else
        label = str{k};
        udata = [];
        cback = '';
    end;
    cm = uimenu('parent',lastmenulvl(lvl(k)-1), 'Label',label, 'Userdata',udata, ...
                'Callback',cback, 'tag','AddedAppMenu');
    lastmenulvl(lvl(k)) = cm;
    if lvl(k) == 2
        toplevelmenus(end+1) = cm;
        toplevelids{end+1}   = id{k};
    end;
end;
% add defaults manipulation entries
for k =1:numel(toplevelmenus)
    cm = uimenu('Parent',toplevelmenus(k), 'Label','Load Defaults', ...
                'Callback',@local_loaddefs, 'Userdata',toplevelids{k}, ...
                'tag','AddedAppMenu', 'Separator','on');
    cm = uimenu('Parent',toplevelmenus(k), 'Label','Save Defaults', ...
                'Callback',@local_savedefs, 'Userdata',toplevelids{k}, ...
                'tag','AddedAppMenu');
    cm = uimenu('parent',toplevelmenus(k), 'Label','Edit Defaults', ...
                'Callback',@local_editdefs, 'Userdata',toplevelids{k}, ...
                'tag','AddedAppMenu', 'Separator','on');
end;
% --------------------------------------------------------------------
function local_addtojob(varargin)
id  = get(gcbo, 'userdata');
handles = guidata(gcbo);
udmodlist = get(handles.modlist, 'userdata');
cfg_util('addtojob', udmodlist.cjob, id);
local_showjob(gcbo);
% --------------------------------------------------------------------
function local_loaddefs(varargin)
appid = get(gcbo, 'Userdata');
[file sts] = cfg_getfile(1, '.*\.m$','Load Defaults from');
if sts
    cfg_util('initdef', appid, file);
end;
% --------------------------------------------------------------------
function local_savedefs(varargin)
appid = get(gcbo, 'Userdata');
[tag, def] = cfg_util('harvestdef', appid);
[file path] = uiputfile({'*.m','MATLAB .m file'}, 'Save Defaults as', ...
                        sprintf('%s_defaults.m', tag));
if ~ischar(file)
    return;
end;
fid = fopen(fullfile(path, file),'w');
if fid < 1
    warning('matlabbatch:cfg_ui:savedefs', ...
            'Save failed: no defaults written to %s.', ...
            fullfile(path, file));
    return;
end;
[defstr tagstr] = gencode(def, tag);
[u1 funcname] = fileparts(file);
fprintf(fid, 'function %s = %s\n', tagstr, funcname);
for k = 1:numel(defstr)
    fprintf(fid, '%s\n', defstr{k});
end;
fclose(fid);
% --------------------------------------------------------------------
function local_editdefs(varargin)
% Defaults edit mode bypasses local_showjob, but uses all other GUI
% callbacks. Where behaviour/GUI visibility is different for
% normal/defaults mode, a check for the presence of udmodlist(1).defid is
% performed.
handles = guidata(gcbo);
% Disable application menus & edit menus
set(findobj(handles.cfg_ui, 'Tag', 'AddedAppMenu'), 'Enable','off');
set(findobj(handles.cfg_ui,'-regexp', 'Tag','.*(Del)|(Repl)Mod$'),'Enable','off');
set(findobj(handles.cfg_ui,'-regexp','Tag','^MenuEditVal.*'), 'Enable', 'off');
% Change current menu to 'Quit'
set(gcbo, 'Enable','on', 'Callback',@local_editdefsquit, ...
          'Label','Quit Defaults');
set(get(gcbo, 'Parent'), 'Enable','on');
% Get module list for application
appid = get(gcbo, 'Userdata');
[id,stop,val]=cfg_util('listcfg', appid, cfg_findspec({{'hidden',false}}), ...
                       {'name'});
udmodlist = get(handles.modlist, 'userdata');
udmodlist(1).defid = id;
set(handles.modlist, 'String',val{1}, 'Value',1, 'Userdata',udmodlist);
local_showmod(gcbo);
% --------------------------------------------------------------------
function local_editdefsquit(varargin)
handles = guidata(gcbo);
set(findobj(handles.cfg_ui, 'Tag', 'AddedAppMenu'), 'Enable','on');
set(gcbo, 'Enable','on', 'Callback',@local_editdefs, ...
          'Label','Edit Defaults');
% remove defs field from udmodlist
udmodlist = rmfield(get(handles.modlist, 'userdata'), 'defid');
if numel(fieldnames(udmodlist)) == 0
    udmodlist = local_init_modlist;
end;
set(handles.modlist, 'userdata',udmodlist);
local_showjob(gcbo);
%% Show job contents
% --------------------------------------------------------------------
function local_showjob(obj,cjob)
handles = guidata(obj);
if nargin == 1
    % udmodlist should be initialised here
    udmodlist = get(handles.modlist,'userdata');
    cjob = udmodlist.cjob;
else
    % set cjob, if supplied
    udmodlist = local_init_modlist;
    udmodlist(1).cjob = cjob;
end;
[id str sts dep sout] = cfg_util('showjob',cjob);
if isempty(str)
    str = {'No Modules in Batch'};
    mrk = {' '};
    cmod = 1;
    udmodlist.cmod = [];
    set(findobj(handles.cfg_ui,'-regexp', 'Tag','.*(Del)|(Repl)Mod$'),'Enable','off');
else
    if isempty(udmodlist.cmod)
        cmod = 1;
    else
        cmod = min(get(handles.modlist,'value'), numel(str));
        if udmodlist.cmod ~= cmod
            set(handles.module, 'Userdata',[]);
        end;
    end
    udmodlist.id = id;
    udmodlist.sout = sout;
    udmodlist.cmod = cmod;
    set(findobj(handles.cfg_ui,'-regexp', 'Tag','.*(Del)|(Repl)Mod$'),'Enable','on');
end;
for k = 1:numel(sts)
    if ~sts(k)
        mrk{k} = '  <-X';
    elseif dep(k)
        mrk{k} = '  DEP';
    else
        mrk{k} = ' ';
    end;
end;
str = cellstr(cat(2,strvcat(str), strvcat(mrk)));
set(handles.modlist, 'string', str, 'userdata',udmodlist, 'value', cmod);
if ~isempty(sts) && all(sts)
    set(findobj(handles.cfg_ui,'-regexp', 'Tag','.*File(Run)|(RunSerial)$'),'Enable','on');
else
    set(findobj(handles.cfg_ui,'-regexp', 'Tag','.*File(Run)|(RunSerial)$'),'Enable','off');
end    
local_showmod(obj);

%% Show Module Contents
% --------------------------------------------------------------------
function local_showmod(obj)
handles = guidata(obj);
udmodlist = get(handles.modlist, 'userdata');
if ~isempty(udmodlist.cmod)
    cmod = get(handles.modlist, 'value');
    % fill module box with module contents
    if isfield(udmodlist, 'defid')
        % list defaults
        dflag = true;
        cid = {udmodlist.defid{cmod} ,[]};
    else
        dflag = false;
        cid = {udmodlist.cjob, udmodlist.id{cmod}, []};
    end;
    [id stop contents] = ...
        cfg_util('listmod', cid{:}, ...
                 cfg_findspec({{'hidden',false}}), ...
                 cfg_tropts({{'hidden', true}},1,Inf,1,Inf,dflag), ...
                 {'name','val','labels','values','class','level', ...
                  'all_set','all_set_item'});
    if isempty(id) || ~cfg_util('isitem_mod_id', id{1})
        % Module not found without hidden flag
        % Try to list top level entry of module anyway, but not module items.
        [id stop contents] = ...
            cfg_util('listmod', cid{:}, ...
                     cfg_findspec({}), ...
                     cfg_tropts({{'hidden', true}},1,1,1,1,dflag), ...
                     {'name','val','labels','values','class','level', ...
                      'all_set','all_set_item'});
    end;
    for k = 1:numel(contents{1})
        indent = repmat('  ', 1, contents{6}{k}-1);
        if contents{8}{k}
            if any(strcmp(contents{5}{k}, {'cfg_menu','cfg_files','cfg_entry'})) && ...
                    isa(contents{2}{k}{1}, 'cfg_dep')
                if numel(contents{2}{k}{1}) == 1
                    datastr{k} = sprintf('DEP %s', contents{2}{k}{1}.sname);
                else
                    datastr{k} = sprintf('DEP (%d outputs)', numel(contents{2}{k}{1}));
                end;
            else
                switch contents{5}{k}
                    case 'cfg_menu',
                        for l = 1:numel(contents{4}{k})
                            if isequal(contents{2}{k}{1}, contents{4}{k}{l})
                                datastr{k} = contents{3}{k}{l};
                                break;
                            end;
                        end;
                    case 'cfg_files',
                        if numel(contents{2}{k}{1}) == 1
                            if length(contents{2}{k}{1}{1}) < 15 && ~isempty(contents{2}{k}{1}{1})
                                datastr{k} = contents{2}{k}{1}{1};
                            elseif isempty(contents{2}{k}{1}{1})
                                datastr{k} = ' ';
                            else
                                datastr{k} = sprintf('...%s', contents{2}{k}{1}{1}(end-10:end));
                            end;
                        else
                            datastr{k} = sprintf('%d files', numel(contents{2}{k}{1}));
                        end;
                    case 'cfg_entry'
                        csz = size(contents{2}{k}{1});
                        % TODO use gencode like string formatting
                        if ischar(contents{2}{k}{1}) && any(csz(1:2) == 1)
                            if numel(contents{2}{k}{1}) < 15
                                datastr{k} = contents{2}{k}{1};
                            else
                                datastr{k} = sprintf('%s...', contents{2}{k}{1}(1:10));
                            end;
                        elseif isnumeric(contents{2}{k}{1}) && any(csz(1:2) == 1)
                            % always display line vector as summary
                            datastr{k} = num2str(contents{2}{k}{1}(:)');
                        else
                            szstr = sprintf('%dx', csz);
                            datastr{k} = sprintf('%s %s', szstr(1:end-1), class(contents{2}{k}{1}));
                        end;
                    otherwise
                        datastr{k} = ' ';
                end;
            end;
        else
            datastr{k} = '<-X';
        end;
        namestr{k} = sprintf('%s%s  ', indent, contents{1}{k});
    end;
    str = cellstr(cat(2, strvcat(namestr), strjust(strvcat(datastr),'right')));
    udmodule = get(handles.module, 'userdata');
    if isempty(udmodule)
        citem = 1;
    else
        % try to find old item in new module struct - this may change due
        % to repeat/choice changes
        oldid = udmodule.id{udmodule.oldvalue};
        citem = 1;
        for k = 1:numel(id)
            if isequal(id{k}, oldid)
                citem = k;
                break;
            end;
        end;
    end;
    udmodule.contents = contents;
    udmodule.id = id;
    udmodule.oldvalue = citem;
    set(handles.module, 'String', str, 'Value', citem, 'userdata', udmodule);
    % set help box to module help
    [id stop help] = cfg_util('listmod', cid{:}, cfg_findspec, ...
                              cfg_tropts(cfg_findspec,1,1,1,1,false), {'help'});
    set(handles.helpbox, 'String', spm_justify(handles.helpbox, help{1}{1}), 'Value', 1);
    udmodlist(1).cmod = cmod;
    set(handles.modlist, 'userdata', udmodlist);
    local_showvaledit(obj);
else
    set(handles.module, 'String','No Module selected', 'Value',1,'Userdata',[]);
    set(findobj(handles.cfg_ui,'-regexp','Tag','^Btn.*'), 'Visible', 'off');
    set(findobj(handles.cfg_ui,'-regexp','Tag','^MenuEditVal.*'), 'Enable', 'off');
    set(handles.valshow, 'String','', 'Visible','off');
    set(handles.valshowLabel, 'Visible','off');
    % set help box to matlabbatch top node help
    [id stop help] = cfg_util('listcfgall', [], cfg_findspec({{'tag','matlabbatch'}}), {'help'});
    set(handles.helpbox, 'String', spm_justify(handles.helpbox, help{1}{1}), 'Value', 1);
end;

%% Show Item
% --------------------------------------------------------------------
function local_showvaledit(obj)
handles = guidata(obj);
udmodlist = get(handles.modlist, 'userdata');
cmod = get(handles.modlist, 'value');
udmodule = get(handles.module, 'userdata');
value = get(handles.module, 'value');
set(findobj(handles.cfg_ui,'-regexp', 'Tag','^BtnVal.*'), 'Visible','off');
set(findobj(handles.cfg_ui,'-regexp', 'Tag','^MenuEditVal.*'), 'Enable','off');
set(handles.valshow,'String', '','Visible','off');
set(handles.valshowLabel, 'Visible','off');
switch(udmodule.contents{5}{value})
    case {'cfg_entry','cfg_files'}
        if ~isempty(udmodule.contents{2}{value}) && isa(udmodule.contents{2}{value}{1}, 'cfg_dep')
            str = {'Reference from'};
            for k = 1:numel(udmodule.contents{2}{value}{1}) % we may have multiple dependencies
                str{k+1} = udmodule.contents{2}{value}{1}(k).sname; % return something to be printed
            end;
        elseif ~isempty(udmodule.contents{2}{value})
            str = gencode(udmodule.contents{2}{value}{1},'val');
        else
            str = '';
        end;
        set(handles.valshow,'String', str, 'Visible','on', 'Value', 1);
        set(handles.valshowLabel, 'Visible','on');
        if ~isfield(udmodlist, 'defid')
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*AddDep$'), ...
                'Visible','on', 'Enable','on');
        end;
        if strcmp(udmodule.contents{5}{value},'cfg_files')
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*SelectFiles$'), ...
                'Visible','on', 'Enable','on');
        else
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*EditVal$'), ...
                'Visible','on', 'Enable','on');
        end
        set(findobj(handles.cfg_ui,'-regexp','Tag','.*ClearVal$'), ...
            'Visible','on', 'Enable','on');
    case 'cfg_menu'
        set(findobj(handles.cfg_ui,'-regexp','Tag','.*EditVal$'), 'Visible','on', 'Enable','on');
        set(findobj(handles.cfg_ui,'-regexp','Tag','.*ClearVal$'), ...
            'Visible','on', 'Enable','on');
    case 'cfg_choice'
        if ~isfield(udmodlist, 'defid')
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*EditVal$'), ...
                'Visible','on', 'Enable','on');
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*ClearVal$'), ...
                'Visible','on', 'Enable','on');
        end;
    case {'cfg_repeat'}
        if ~isfield(udmodlist, 'defid')
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*AddItem$'), ...
                'Visible','on', 'Enable','on');
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*ReplItem$'), ...
                'Visible','on', 'Enable','on');
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*DelItem$'), ...
                'Visible','on', 'Enable','on');
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*ClearVal$'), ...
                'Visible','on', 'Enable','on');
        end;
end;
if isfield(udmodlist, 'defid')
    cmid = {udmodlist.defid{cmod}};
else
    cmid = {udmodlist.cjob udmodlist.id{cmod}};
end;
[id stop help] = cfg_util('listmod', cmid{:}, udmodule.id{value}, cfg_findspec, ...
                          cfg_tropts(cfg_findspec,1,1,1,1,false), {'help'});
set(handles.helpbox, 'string', spm_justify(handles.helpbox, help{1}{1}));

%% Value edit dialogues
% --------------------------------------------------------------------

function local_valedit_edit(hObject)
% Normal mode. Depending on strtype, put '' or [] around entered
% input. If input has ndims > 2, isn't numeric or char, proceed with
% expert dialogue.
handles = guidata(hObject);
if strcmp(get(findobj(handles.cfg_ui,'Tag','MenuEditExpertEdit'), ...
              'checked'),'on')
    local_valedit_expert_edit(hObject);
    return;
end;
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
cmod = get(handles.modlist, 'Value');
udmodlist = get(handles.modlist, 'Userdata');
val = udmodule.contents{2}{value};
if isfield(udmodlist, 'defid')
    cmid = {udmodlist.defid{cmod}};
else
    cmid = {udmodlist.cjob, udmodlist.id{cmod}};
end;
[id stop strtype] = cfg_util('listmod', cmid{:}, udmodule.id{value}, cfg_findspec, ...
                             cfg_tropts(cfg_findspec,1,1,1,1,false), {'strtype'});
if isempty(val)
    if strtype{1}{1} == 's'
        val = {''};
    else
        val = {[]};
    end;
end;
% If we can't handle this, use expert mode
if ndims(val{1}) > 2 || ~(ischar(val{1}) || isnumeric(val{1}))
    local_valedit_expert_edit(hObject);
    return;
end
if strtype{1}{1} == 's'
    instr = val;
    encl  = '''''';
else
    try
        instr = {num2str(val{1})};
        encl  = '[]';
    catch
        local_valedit_expert_edit(hObject);
        return;
    end;
end;
sts = false;
while ~sts
    % estimate size of input field based on instr
    % Maximum widthxheight 140x20, minium 60x2
    szi = size(instr{1});
    mxwidth = 140;
    rdup = ceil(szi(2)/mxwidth)+3;
    szi = max(min([szi(1)*rdup szi(2)],[20 140]),[2,60]);
    str = inputdlg(strvcat('Enter a value.', ...
                           ' ', ...
                           'To clear a value, clear the input field and accept.', ...
                           ' ', ...
                           ['Accept input with CTRL-RETURN, cancel with ' ...
                        'ESC.']), ...
                   udmodule.contents{1}{value}, ...
                   szi,instr);
    if iscell(str) && isempty(str)
        % User has hit cancel button
        return;
    end;
    % save instr in case of evaluation error
    instr = str;
    % str{1} is a multiline char array
    % 1) cellify it
    % 2) add newline to each string
    % 3) concatenate into one string
    cstr = cellstr(str{1});
    str = strcat(cstr, {char(10)});
    str = cat(2, str{:});
    % Evaluation is encapsulated to avoid users compromising this function
    % context
    [val sts] = local_eval_valedit(str);
    if ~sts
        % try with matching value enclosure
        if strtype{1}{1} == 's'
            str = strcat({encl(1)}, cstr, {encl(2)}, {char(10)});
        else
            cestr = {encl(1) cstr{:} encl(2)};
            str = strcat(cestr, {char(10)});
        end;
        str = cat(2, str{:});
        % Evaluation is encapsulated to avoid users compromising this function
        % context
        [val sts] = local_eval_valedit(str);
        if ~sts
            uiwait(msgbox(sprintf('Input could not be evaluated.'),'Evaluation error','modal'));
        end;
    end;
end;
% This code will only be reached if a new value has been set
local_setvaledit(hObject, val);

function local_valedit_expert_edit(hObject)
% Expert mode, use full flexibility of evaluated expressions.
handles = guidata(hObject);
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
val = udmodule.contents{2}{value};
sts = false;
if ~isempty(val) && isa(val{1}, 'cfg_dep')
    local_valedit_AddDep(hObject);
    return;
else
    instr = {strvcat(get(handles.valshow, 'String'))};
end;
while ~sts
    % estimate size of input field based on instr
    % Maximum widthxheight 140x20, minium 60x2
    szi = size(instr{1});
    mxwidth = 140;
    rdup = ceil(szi(2)/mxwidth)+3;
    szi = max(min([szi(1)*rdup szi(2)],[20 140]),[2,60]);
    str = inputdlg(strvcat('Enter a valid MATLAB expression.', ...
                           ' ', ...
                           ['Strings must be enclosed in single quotes ' ...
                        '(''A''), multiline arrays in brackets ([ ]).'], ...
                           ' ', ...
                           'To clear a value, enter an empty cell ''{}''.', ...
                           ' ', ...
                           ['Accept input with CTRL-RETURN, cancel with ' ...
                        'ESC.']), ...
                   udmodule.contents{1}{value}, ...
                   szi,instr);
    if iscell(str) && isempty(str)
        % User has hit cancel button
        return;
    end;
    % save instr in case of evaluation error
    instr = str;
    % str{1} is a multiline char array
    % 1) cellify it
    % 2) add newline to each string
    % 3) concatenate into one string
    str = strcat(cellstr(str{1}), {char(10)});
    str = cat(2, str{:});
    % Evaluation is encapsulated to avoid users compromising this function
    % context
    [val sts] = local_eval_valedit(str);
    if ~sts
        uiwait(msgbox(sprintf('Input could not be evaluated. Possible reasons are:\n1) Input should be a vector or matrix, but is not enclosed in ''['' and '']'' brackets.\n2) Input should be a character or string, but is not enclosed in '' single quotes.\n3) Input should be a MATLAB variable, but is misspelled.\n4) Input should be a MATLAB expression, but has syntax errors.'),'Evaluation error','modal'));
    end;
end;
% This code will only be reached if a new value has been set
local_setvaledit(hObject, val);

% --------------------------------------------------------------------
function local_valedit_list(hObject)

handles = guidata(hObject);
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
cval = -1;
if strcmp(udmodule.contents{5}{value},'cfg_choice')
    % compare tag, not filled entries
    cmpsubs = substruct('.','tag');
else
    cmpsubs =struct('type',{},'subs',{});
end;
valsubs = substruct('{}',{1});
for l = 1:numel(udmodule.contents{4}{value})
    valuesubs = substruct('{}',{l});
    if ~isempty(udmodule.contents{2}{value}) && isequal(subsref(udmodule.contents{2}{value},[valsubs cmpsubs]), subsref(udmodule.contents{4}{value},[valuesubs cmpsubs]))
        mrk{l} = '*';
        cval = l;
    else
        mrk{l} = ' ';
    end;
end;
if strcmp(udmodule.contents{5}{value},'cfg_choice')
    for k = 1:numel(udmodule.contents{4}{value})
        str{k} = udmodule.contents{4}{value}{k}.name;
    end;
else
    str = udmodule.contents{3}{value};
end;
if cval > 0
    inival = cval;
else
    inival = 1;
end;
% for some reason, MATLAB wants to have the list size to be specified in
% pixels, not chars. Try to work this out based on the number of
% characters and a font size of 12.
szi = min(size(strvcat(str)')+1,[140 60])*13;
[val sts] = listdlg('Name',udmodule.contents{1}{value}, ...
                    'ListString',strcat(mrk, str), 'SelectionMode','single', ...
                    'InitialValue',inival, 'ListSize', szi);
if ~sts || val == cval
    % Selection cancelled or nothing changed
    return;
end;
local_setvaledit(hObject, val);

%% Set changed values
% --------------------------------------------------------------------
function local_setvaledit(obj, val)
handles = guidata(obj);
udmodlist = get(handles.modlist, 'Userdata');
cmod = get(handles.modlist, 'Value');
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
if isfield(udmodlist, 'defid')
    cfg_util('setdef', udmodlist(1).defid{cmod}, udmodule.id{citem}, val);
    local_showmod(obj);
else
    cfg_util('setval', udmodlist.cjob, udmodlist.id{cmod}, udmodule.id{citem}, val);
    cfg_util('harvest', udmodlist.cjob, udmodlist.id{cmod});
    local_showjob(obj);
end;

%% Button/Menu callbacks
% These callbacks are called from both the Edit menu and the GUI buttons
% --------------------------------------------------------------------
function local_valedit_SelectFiles(hObject)
handles = guidata(hObject);
udmodlist = get(handles.modlist, 'Userdata');
cmod = get(handles.modlist, 'Value');
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
if isfield(udmodlist,'defid')
    cmid = {udmodlist.defid{cmod}};
else
    cmid = {udmodlist.cjob, udmodlist.id{cmod}};
end;
[unused1 unused2 contents] = cfg_util('listmod', cmid{:}, udmodule.id{citem},{},cfg_tropts({},1,1,1,1,false),{'val','num','filter','name','dir','ufilter'});
if isempty(contents{1}{1}) || isa(contents{1}{1}{1}, 'cfg_dep')
    inifile = '';
else
    inifile = contents{1}{1}{1};
end;
[val sts] = cfg_getfile(contents{2}{1}, contents{3}{1}, contents{4}{1}, inifile, contents{5}{1}, contents{6}{1});
if sts
    local_setvaledit(hObject, cellstr(val));
end;

% --------------------------------------------------------------------
function local_valedit_EditValue(hObject)
handles = guidata(hObject);
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
if any(strcmp(udmodule.contents{5}{value},{'cfg_choice','cfg_menu'}))
    local_valedit_list(hObject);
else
    local_valedit_edit(hObject);
end;

% --------------------------------------------------------------------
function local_valedit_ClearVal(hObject)
local_setvaledit(hObject, {});

% --------------------------------------------------------------------
function local_valedit_AddDep(hObject)
handles = guidata(hObject);
udmodlist = get(handles.modlist, 'Userdata');
cmod = get(handles.modlist, 'Value');
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
sout = cat(2,udmodlist(1).sout{1:cmod-1});
smatch = false(size(sout));
% loop over sout to find all souts that match the current item
for k = 1:numel(sout)
    smatch(k) = cfg_util('match', udmodlist.cjob, udmodlist.id{cmod}, udmodule.id{citem}, sout(k).tgt_spec);
end;
sout = sout(smatch);
if isempty(sout)
    set(findobj(handles.cfg_ui,'-regexp','Tag','.*AddDep$'), 'Enable','off');
    uiwait(msgbox('No matching dependencies found.','Add Dependency','modal'));
else
    % for some reason, MATLAB wants to have the list size to be specified in
    % pixels, not chars. Try to work this out based on the number of
    % characters and a font size of 12.
    str = strvcat(sout.sname);
    szi = min(size(strvcat(str)')+1,[140 20])*13;
    [val sts] = listdlg('Name',udmodule.contents{1}{citem}, 'ListString',str, ...
                        'ListSize',szi);
    if sts
        local_setvaledit(hObject, sout(val));
    end;
end;

% --------------------------------------------------------------------
function local_valedit_AddItem(hObject)
handles = guidata(hObject);
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
if numel(udmodule.contents{4}{citem}) > 1
    % If there is a choice, prepare a selection list
    for k = 1:numel(udmodule.contents{4}{citem})
        str{k} = udmodule.contents{4}{citem}{k}.name;
    end;
    [val sts] = listdlg('Name', udmodule.contents{1}{citem}, 'ListString',str, 'SelectionMode','single');
else
    % There is no choice, just add the only item
    val = 1;
    sts = true;
end
if sts
    local_setvaledit(hObject,[val inf]);
end;

% --------------------------------------------------------------------
function local_valedit_ReplItem(hObject)
handles = guidata(hObject);
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
if numel(udmodule.contents{2}{citem}) > 1
    % If there is a choice, prepare a selection list
    for k = 1:numel(udmodule.contents{2}{citem})
        str{k} = sprintf('%s (%d)', udmodule.contents{2}{citem}{k}.name, k);
    end;
    [val sts] = listdlg('Name', udmodule.contents{1}{citem}, 'ListString',str, 'SelectionMode','single');
else
    % There is no choice, just replicate the only item
    val = 1;
    sts = true;
end
if sts
    local_setvaledit(hObject,[-1 val]);
end;

% --------------------------------------------------------------------
function local_valedit_DelItem(hObject)
handles = guidata(hObject);
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
str = {};
for k = 1:numel(udmodule.contents{2}{citem})
    str{k} = udmodule.contents{2}{citem}{k}.name;
end;
if ~isempty(str)
    [val sts] = listdlg('Name', udmodule.contents{1}{citem}, 'ListString', str);
    if sts
        % delete from last to first item, otherwise order would be not
        % preserved
        val = sort(val, 'descend');
        for k = 1:numel(val)
            local_setvaledit(hObject, [Inf val(k)]);
        end;
    end;
end;

%% Automatic Callbacks
% --- Executes just before cfg_ui is made visible.
function cfg_ui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cfg_ui (see VARARGIN)

% Add configuration specific menu items
local_setmenu(hObject);

% Check udmodlist
udmodlist = get(handles.modlist, 'userdata');
if isempty(udmodlist)
    udmodlist = local_init_modlist;
    udmodlist.cjob = cfg_util('initjob');
    set(handles.modlist, 'userdata', udmodlist);
end;

% show job
local_showjob(hObject);

% Choose default command line output for cfg_ui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cfg_ui wait for user response (see UIRESUME)
% uiwait(handles.cfg_ui);


% --- Outputs from this function are returned to the command line.
function varargout = cfg_ui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% File Menu Callbacks
% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MenuFileNew_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist, 'userdata');
if ~isempty(udmodlist.cmod)
    cfg_util('deljob',udmodlist(1).cjob);
end;
udmodlist = local_init_modlist;
udmodlist.cjob = cfg_util('initjob');
set(handles.modlist, 'userdata', udmodlist);
local_showjob(hObject);

% --------------------------------------------------------------------
function MenuFileLoad_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[files sts] = cfg_getfile([1 Inf], '.*\.m$', 'Load Job File(s)');
if sts
    udmodlist = get(handles.modlist, 'userdata');
    cfg_util('deljob',udmodlist(1).cjob);
    udmodlist = local_init_modlist;
    udmodlist.cjob = cfg_util('initjob', cellstr(files));
    set(handles.modlist, 'userdata', udmodlist);
    local_showjob(hObject);
end;

% --------------------------------------------------------------------
function MenuFileSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file path idx] = uiputfile({'*.m','Matlab m file'}, 'Save Job');
if isnumeric(file) && file == 0
    return;
end;
udmodlist = get(handles.modlist, 'userdata');
cfg_util('savejob', udmodlist.cjob, fullfile(path, file));

% --------------------------------------------------------------------
function MenuFileRun_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(get(0,'Children'),'Pointer','watch');
udmodlist = get(handles.modlist, 'userdata');
cfg_util('run',udmodlist(1).cjob);
set(get(0,'Children'),'Pointer','arrow');

% --------------------------------------------------------------------
function MenuFileRunSerial_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileRunSerial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(get(0,'Children'),'Pointer','watch');
udmodlist = get(handles.modlist, 'userdata');
cfg_util('runserial',udmodlist(1).cjob);
set(get(0,'Children'),'Pointer','arrow');

% --------------------------------------------------------------------
function MenuFileAddApp_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileAddApp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file path idx] = uigetfile({'*.m','Matlab m file'}, 'Load Application Configuration');
if isnumeric(file) && file == 0
    return;
end;
addpath(path);
[p fun e v] = fileparts(file);
cfg_util('addapp', fun);
local_setmenu(hObject);

% --------------------------------------------------------------------
function MenuFileClose_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.cfg_ui);

%% Edit Menu Callbacks
% --------------------------------------------------------------------
function MenuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MenuEditUpdateView_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditUpdateView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_showjob(hObject);

% --------------------------------------------------------------------
function MenuEditReplMod_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditReplMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_ReplMod(hObject);

% --------------------------------------------------------------------
function MenuEditDelMod_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditDelMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_DelMod(hObject);

% --------------------------------------------------------------------
function MenuEditValAddItem_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValAddItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_AddItem(hObject);

% --------------------------------------------------------------------
function MenuEditValReplItem_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValReplItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_ReplItem(hObject);

% --------------------------------------------------------------------
function MenuEditValDelItem_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValDelItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_DelItem(hObject);

% --------------------------------------------------------------------
function MenuEditValEditVal_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValEditVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_EditValue(hObject);

% --------------------------------------------------------------------
function MenuEditValClearVal_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValClearVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_ClearVal(hObject);

% --------------------------------------------------------------------
function MenuEditValSelectFiles_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValSelectFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_SelectFiles(hObject);

% --------------------------------------------------------------------
function MenuEditValAddDep_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValAddDep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_AddDep(hObject);

% --------------------------------------------------------------------
function MenuEditExpertEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditExpertEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(gcbo,'checked'),'on')
    set(gcbo, 'checked', 'off');
else
    set(gcbo, 'checked', 'on');
end;

%% Module List Callbacks
% --- Executes on selection change in modlist.
function modlist_Callback(hObject, eventdata, handles)
% hObject    handle to modlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns modlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modlist
udmodlist = get(handles.modlist, 'userdata');
if ~isempty(udmodlist.cmod)
    if ~isfield(udmodlist, 'defid')
        cfg_util('harvest', udmodlist(1).id{udmodlist(1).cmod});
        local_showjob(hObject);
    else
        local_showmod(hObject);
    end;
end;

% --- Executes during object creation, after setting all properties.
function modlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Module List Context Menu Callbacks

% --------------------------------------------------------------------
function CmModlistReplMod_Callback(hObject, eventdata, handles)
% hObject    handle to CmModlistReplMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_ReplMod(hObject);

% --------------------------------------------------------------------
function CmModlistDelMod_Callback(hObject, eventdata, handles)
% hObject    handle to CmModlistDelMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_DelMod(hObject);

% --------------------------------------------------------------------
function CmModlist_Callback(hObject, eventdata, handles)
% hObject    handle to CmModlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Module Callbacks
% --- Executes on selection change in module.
function module_Callback(hObject, eventdata, handles)
% hObject    handle to module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns module contents as cell array
%        contents{get(hObject,'Value')} returns selected item from module

% Selection change is called both when there is a real selection change,
% but also if return is hit or there is a double click

value = get(hObject,'Value');
udmodule = get(hObject,'Userdata');
if isempty(udmodule)
    return;
end;
if udmodule.oldvalue ~= value
    udmodule.oldvalue = value;
    set(hObject, 'Userdata', udmodule);
    local_showvaledit(hObject);
end;
if strcmp(get(handles.cfg_ui,'SelectionType'),'open')
    % open modal MenuEdit window, do editing
    % Unfortunately, MATLAB focus behaviour makes it impossible to do this
    % in the existing valshow object - if this object looses focus, it will
    % call its callback without checking why it lost focus.
    % Call appropriate input handler for editable and selection types
    switch udmodule.contents{5}{value}
        case {'cfg_entry'},
            local_valedit_edit(hObject);
        case { 'cfg_files'},
            local_valedit_SelectFiles(hObject);
        case {'cfg_choice', 'cfg_menu'},
            local_valedit_list(hObject);
    end;
end;

% --- Executes during object creation, after setting all properties.
function module_CreateFcn(hObject, eventdata, handles)
% hObject    handle to module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Value Display Callbacks
% --- Executes during object creation, after setting all properties.
function valshow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: MenuEdit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function valshow_Callback(hObject, eventdata, handles)
% hObject    handle to valshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valshow as text
%        str2double(get(hObject,'String')) returns contents of valshow as a double

%% GUI Buttons
% --- Executes on button press in BtnValSelectFiles.
function BtnValSelectFiles_Callback(hObject, eventdata, handles)
% hObject    handle to BtnValSelectFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_SelectFiles(hObject);

% --- Executes on button press in BtnValEditVal.
function BtnValEditVal_Callback(hObject, eventdata, handles)
% hObject    handle to BtnValEditVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_EditValue(hObject);

% --- Executes on button press in BtnValAddDep.
function BtnValAddDep_Callback(hObject, eventdata, handles)
% hObject    handle to BtnValAddDep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_AddDep(hObject);

% --- Executes on button press in BtnValAddItem.
function BtnValAddItem_Callback(hObject, eventdata, handles)
% hObject    handle to BtnValAddItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_AddItem(hObject);

% --- Executes on button press in BtnValDelItem.
function BtnValDelItem_Callback(hObject, eventdata, handles)
% hObject    handle to BtnValDelItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_DelItem(hObject);


% --- Executes on button press in BtnValReplItem.
function BtnValReplItem_Callback(hObject, eventdata, handles)
% hObject    handle to BtnValReplItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_ReplItem(hObject);


% --------------------------------------------------------------------
function modlist = local_init_modlist
% Initialise modlist to empty struct
% Don't initialise defid field - this will be added by defaults editor
modlist = struct('cjob',[],'cmod',[],'id',[],'sout',[]);
