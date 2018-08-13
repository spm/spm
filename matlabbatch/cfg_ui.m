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
% $Id: cfg_ui.m 7394 2018-08-13 16:24:53Z spm $

rev = '$Rev: 7394 $'; %#ok

% edit the above text to modify the response to help cfg_ui

% Last Modified by GUIDE v2.5 18-Jul-2014 13:38:10

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

%%% Data structure initialisation
% --------------------------------------------------------------------
function udmodlist = local_init_udmodlist
% Initialise udmodlist to empty struct
% Don't initialise defid field - this will be added by defaults editor
udmodlist = struct('cjob',[],'cmod',[],'id',[],'sout',[],'modified',false,'wd','');

%%% Module management
% --------------------------------------------------------------------
function local_DelMod(hObject)
handles = guidata(hObject);
udmodlist = get(handles.modlist,'userdata');
val = get(handles.modlist,'value');
if ~isempty(udmodlist.cmod)
    cfg_util('delfromjob',udmodlist.cjob, udmodlist.id{val});
    udmodlist.modified = true;
    set(handles.modlist,'userdata',udmodlist);
    local_showjob(hObject);
end;

% --------------------------------------------------------------------
function local_ReplMod(hObject)
handles = guidata(hObject);
udmodlist = get(handles.modlist,'userdata');
val = get(handles.modlist,'value');
if ~isempty(udmodlist.cmod)
    cfg_util('replicate',udmodlist.cjob, udmodlist.id{val});
    udmodlist.modified = true;
    set(handles.modlist,'userdata',udmodlist);
    local_showjob(hObject);
end;

% --------------------------------------------------------------------
function local_AddMod(varargin)
id  = get(gcbo, 'userdata');
handles = guidata(gcbo);
udmodlist = get(handles.modlist, 'userdata');
% add module to job, harvest to initialise its virtual outputs
mod_job_id = cfg_util('addtojob', udmodlist.cjob, id);
cfg_util('harvest', udmodlist.cjob, mod_job_id);
udmodlist.modified = true;
set(handles.modlist,'userdata',udmodlist);
local_showjob(gcbo);

%%% GUI management
% --------------------------------------------------------------------
function local_setmenu(parent, id, cb, dflag)
% parent: menu parent
% id: id to start listcfgall
% cb: callback to add/run new nodes
% dflag: add defaults edit (true/false)

% delete previous added menu, if any
prevmenu = findobj(parent, 'Tag', 'AddedAppMenu');
if ~isempty(prevmenu)
    delete(prevmenu);
end;
% get strings and ids
[id,stop,val]=cfg_util('listcfgall',id,cfg_findspec({{'hidden',false}}),{'name','level'});
str = val{1};
lvl = cat(1,val{2}{:});
% strings and ids are in preorder - if stop is true, then we are at leaf
% level, if lvl(k) <= lvl(k-1) we are at siblings/parent level
% remember last menu at lvl for parents, start at lvl > 1 (top parent is
% figure menu)
lastmenulvl = zeros(1, max(lvl));
lastmenulvl(1) = parent;
toplevelmenus = [];
toplevelids   = {};
% 1st entry is top of matlabbatch config tree, applications start at 2nd entry
for k = 2:numel(lvl)
    label = str{k};
    if stop(k)
        udata = id{k};
        cback = cb;
    else
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
hkeys = cell(1,numel(toplevelmenus)+2);
hkeys{1} = 'f';
hkeys{2} = 'e';
for k =1:numel(toplevelmenus)
    % add hot keys
    clabel = get(toplevelmenus(k),'Label');
    for l = 1:numel(clabel)
        if ~isspace(clabel(l)) && ~any(strcmpi(clabel(l),hkeys))
            hkeys{k+2} = lower(clabel(l));
            clabel = [clabel(1:l-1) '&' clabel(l:end)];
            break;
        end;
    end;
    set(toplevelmenus(k),'Label',clabel);
    if dflag
        % add defaults manipulation entries
        % disable Load/Save
        %cm = uimenu('Parent',toplevelmenus(k), 'Label','Load Defaults', ...
        %           'Callback',@local_loaddefs, 'Userdata',toplevelids{k}, ...
        %           'tag','AddedAppMenu', 'Separator','on');
        %cm = uimenu('Parent',toplevelmenus(k), 'Label','Save Defaults', ...
        %           'Callback',@local_savedefs, 'Userdata',toplevelids{k}, ...
        %           'tag','AddedAppMenu');
        cm = uimenu('parent',toplevelmenus(k), 'Label','Edit Defaults', ...
                    'Callback',@local_editdefs, 'Userdata',toplevelids{k}, ...
                    'tag','AddedAppMenu', 'Separator','on');
    end;
end;

% --------------------------------------------------------------------
function local_setfont(obj,fs)
handles = guidata(obj);
cfg_get_defaults([mfilename '.lfont'], fs);
set(handles.modlist, fs{:});
set(handles.module, fs{:});
set(handles.valshow, fs{:});
set(handles.helpbox, fs{:});

% --------------------------------------------------------------------
function local_pointer(ptr)
shh = get(0,'showhiddenhandles');
set(0,'showhiddenhandles','on');
C = get(0,'Children');
for i = 1:numel(C)
    try, set(C(i),'Pointer',ptr); end
end
drawnow;
set(0,'showhiddenhandles',shh);

%%% Defaults editing
% --------------------------------------------------------------------
function local_loaddefs(varargin)
appid = get(gcbo, 'Userdata');
[file, sts] = cfg_getfile(1, '.*\.m$','Load Defaults from');
if sts
    cfg_util('initdef', appid, file{1});
end;

% --------------------------------------------------------------------
function local_savedefs(varargin)
appid = get(gcbo, 'Userdata');
[tag, def] = cfg_util('harvestdef', appid);
[file, path] = uiputfile({'*.m','MATLAB .m file'}, 'Save Defaults as', ...
                        sprintf('%s_defaults.m', tag));
if ~ischar(file)
    return;
end;
fname      = fullfile(path, file);
[fid, msg] = fopen(fname, 'wt');
if fid == -1
    cfg_message('matlabbatch:fopen', 'Failed to open ''%s'' for writing:\n%s', fname, msg);
end
[defstr, tagstr] = gencode(def, tag);
[u1, funcname] = fileparts(file);
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
% Disable application menus & file, edit menus
set(findobj(handles.cfg_ui, 'Tag', 'AddedAppMenu'), 'Enable','off');
set(findobj(handles.cfg_ui, 'Tag', 'MenuFile'), 'Enable', 'off');
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
udmodlist(1).cmod  = 1;
set(handles.modlist, 'Value',1, 'ListboxTop',1, 'Userdata',udmodlist, 'String',val{1});
local_showmod(gcbo);

% --------------------------------------------------------------------
function local_editdefsquit(varargin)
handles = guidata(gcbo);
set(findobj(handles.cfg_ui, 'Tag', 'AddedAppMenu'), 'Enable','on');
set(findobj(handles.cfg_ui, 'Tag', 'MenuFile'), 'Enable', 'on');
set(gcbo, 'Enable','on', 'Callback',@local_editdefs, ...
          'Label','Edit Defaults');
% remove defs field from udmodlist
udmodlist = rmfield(get(handles.modlist, 'userdata'), 'defid');
if numel(fieldnames(udmodlist)) == 0
    udmodlist = local_init_udmodlist;
end;
set(handles.modlist, 'userdata',udmodlist);
local_showjob(gcbo);

%%% Contents display
% --------------------------------------------------------------------
function local_showjob(obj,cjob)
handles = guidata(obj);
if nargin == 1
    % udmodlist should be initialised here
    udmodlist = get(handles.modlist,'userdata');
    cjob = udmodlist.cjob;
else
    % set cjob, if supplied
    udmodlist = local_init_udmodlist;
    udmodlist(1).cjob = cjob;
    % move figure onscreen
    cfg_onscreen(obj);
    set(obj,'Visible','on');
end;
[id, str, sts, dep, sout] = cfg_util('showjob',cjob);
if isempty(str)
    str = {'No Modules in Batch'};
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
    mrk = cell(size(sts));
    [mrk{dep}] = deal('DEP');
    [mrk{~sts}] = deal('<-X');
    [mrk{~dep & sts}] = deal('');
    str = cfg_textfill(handles.modlist, str, mrk, false);
end;
ltop = cfg_ui_getListboxTop(handles.modlist, cmod, numel(str));
set(handles.modlist, 'userdata',udmodlist, 'value', cmod, 'ListboxTop', ltop, 'string', str);
if ~isempty(sts) && all(sts)
    set(findobj(handles.cfg_ui,'-regexp', 'Tag','.*File(Run)|(RunSerial)$'),'Enable','on');
else
    set(findobj(handles.cfg_ui,'-regexp', 'Tag','.*File(Run)|(RunSerial)$'),'Enable','off');
end    
local_showmod(obj);

% --------------------------------------------------------------------
function local_showmod(obj)
handles = guidata(obj);
udmodlist = get(handles.modlist, 'userdata');
if ~isempty(udmodlist.cmod)
    cmod = get(handles.modlist, 'value');
    % fill module box with module contents
    dflag = isfield(udmodlist, 'defid');
    if isfield(udmodlist, 'defid')
        % list defaults
        cmid = udmodlist.defid{cmod};
    else
        cmid = {udmodlist.cjob, udmodlist.id{cmod}};
    end;
    [id, namestr, datastr, contents] = cfg_ui_util('showmod', cmid, dflag);
    str = cfg_textfill(handles.module,namestr,datastr,true);
    udmodule = get(handles.module, 'userdata');
    if isempty(udmodule)
        citem = min(2,numel(id));
    else
        % try to find old item in new module struct - this may change due
        % to repeat/choice changes
        % Try to make first input item active, not module head
        oldid = udmodule.id{udmodule.oldvalue};
        citem = find(cellfun(@(cid)isequal(cid, oldid),id));
        if isempty(citem)
            citem = min(2,numel(id));
        end
    end;
    udmodule.contents = contents;
    udmodule.id = id;
    udmodule.oldvalue = citem;
    ltop = cfg_ui_getListboxTop(handles.module, citem, numel(str));
    set(handles.moduleHead,'String',sprintf('Current Module: %s', contents{1}{1}));
    set(handles.module, 'Value', citem, 'ListboxTop', ltop, 'userdata', udmodule, 'String', str);
    udmodlist(1).cmod = cmod;
    set(handles.modlist, 'userdata', udmodlist);
    local_showvaledit(obj);
    uicontrol(handles.module);
else
    set(handles.module, 'Value',1,'ListboxTop',1,'Userdata',[], 'String',{'No Module selected'});
    set(handles.moduleHead,'String','No Current Module');
    set(findobj(handles.cfg_ui,'-regexp','Tag','^Btn.*'), 'Visible', 'off');
    set(findobj(handles.cfg_ui,'-regexp','Tag','^MenuEditVal.*'), 'Enable', 'off');
    set(handles.valshow, 'String','', 'Visible','off');
    set(handles.valshowLabel, 'Visible','off');
    % set help box to matlabbatch top node help
    [id, stop, help] = cfg_util('listcfgall', [], cfg_findspec({{'tag','matlabbatch'}}), {'showdoc'});
    set(handles.helpbox, 'Value',1, 'ListboxTop',1, 'String',cfg_justify(handles.helpbox, help{1}{1}));
end;

% --------------------------------------------------------------------
function local_showvaledit(obj)
handles = guidata(obj);
fig = handles.cfg_ui;
udmodlist = get(handles.modlist, 'userdata');
udmodule = get(handles.module, 'userdata');
cmod = get(handles.modlist, 'value');
dflag = isfield(udmodlist, 'defid');
citem = get(handles.module, 'value');
if dflag
    ciid = [udmodlist.defid(cmod) udmodule.id(citem)];
else
    ciid = {udmodlist.cjob udmodlist.id{cmod} udmodule.id{citem}};
end;
contents = cellfun(@(c)subsref(c, substruct('{}',{citem})), udmodule.contents, 'UniformOutput', false);
sout = cat(2, udmodlist.sout{1:cmod-1});
cfg_ui_util('showvaledit', fig, ciid, contents, sout, dflag, [], @()local_valedit_update(obj));
drawnow;

%%% Value edit dialogues
% --------------------------------------------------------------------
function local_valedit_EditValue(hObject)
handles = guidata(hObject);
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
cmod = get(handles.modlist, 'Value');
udmodlist = get(handles.modlist, 'Userdata');
val = udmodule.contents{2}{value};
dflag = isfield(udmodlist, 'defid');
if dflag
    ciid = {udmodlist.defid(cmod) udmodule.id{value}};
else
    ciid = {udmodlist.cjob, udmodlist.id{cmod} udmodule.id{value}};
end;
itemname = udmodule.contents{1}{value};
cfg_ui_util('valedit_EditValue', ciid, itemname, val);

% --------------------------------------------------------------------
function local_valedit_update(hObject)
% update GUI after some value has changed
handles = guidata(hObject);
udmodlist = get(handles.modlist, 'Userdata');
dflag = isfield(udmodlist, 'defid');
if dflag
    local_showmod(hObject);
else
    udmodlist.modified = true;
    set(handles.modlist,'userdata',udmodlist);
    local_showjob(hObject);
end;

% --------------------------------------------------------------------
function local_valedit_dep(hObject)
handles = guidata(hObject);
udmodlist = get(handles.modlist, 'userdata');
cmod = get(handles.modlist, 'value');
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
sout = cfg_ui_util('showvaldeps', udmodlist.cjob, udmodlist.id{cmod}, udmodule.id{citem}, cat(2,udmodlist(1).sout{1:cmod-1}));
str = {sout.sname};
[val, sts] = listdlg('Name',udmodule.contents{1}{citem}, 'ListString',str);
if sts
    cfg_ui_util('setvaledit', {udmodlist.cjob, udmodlist.id{cmod}, udmodule.id{citem}}, sout(val), false);
end;
local_valedit_update(hObject);

% --------------------------------------------------------------------
function local_clearvaledit(hObject)
[ciid, dflag] = local_get_ciid(hObject);
cfg_ui_util('setvaledit', ciid, {}, dflag);
local_valedit_update(hObject);

% --------------------------------------------------------------------
function local_preview(hObject)
[ciid, dflag] = local_get_ciid(hObject);
cfg_ui_util('preview', ciid, dflag);

% --------------------------------------------------------------------
function [ciid, dflag] = local_get_ciid(hObject)
handles = guidata(hObject);
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
cmod = get(handles.modlist, 'Value');
udmodlist = get(handles.modlist, 'Userdata');
dflag = isfield(udmodlist, 'defid');
if dflag
    ciid = {udmodlist.defid(cmod) udmodule.id{value}};
else
    ciid = {udmodlist.cjob, udmodlist.id{cmod} udmodule.id{value}};
end;

% --------------------------------------------------------------------
function cmd = local_check_job_modified(udmodlist, action)
switch lower(action)
    case 'replace'
        queststr = ['The current batch contains unsaved changes. '...
            'Do you want to replace it with another batch?'];
        answers  = {'Continue','Cancel'};
        defans   = 'Continue';
    case 'quit'
        queststr = ['The current batch contains unsaved changes. Do you want to quit ' ...
            'anyway or do you want to hide the batch window ' ...
            'instead?'];
        answers  = {'Quit','Cancel','Hide'};
        defans   = 'Quit';
end
if udmodlist.modified
    cmd = questdlg(queststr, 'Unsaved Changes', answers{:}, defans);
else
    cmd = defans;
end;

%%% Code display callbacks
% --------------------------------------------------------------------
function local_ShowCode_Copy(ob, ev, ctxt)
str = get(ctxt,'String');
sel = get(ctxt,'Value');
str = str(sel);
clipboard('copy',sprintf('%s\n',str{:}));

% --------------------------------------------------------------------
function local_ShowCode_SelAll(ob, ev, ctxt)
set(ctxt,'Value', 1:numel(get(ctxt, 'String')));

% --------------------------------------------------------------------
function local_ShowCode_UnSelAll(ob, ev, ctxt)
set(ctxt,'Value', []);

%% Automatic Callbacks
% --- Executes just before cfg_ui is made visible.
function cfg_ui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cfg_ui (see VARARGIN)

% move figure onscreen
cfg_onscreen(hObject);

% Add configuration specific menu items
local_setmenu(handles.cfg_ui, [], @local_AddMod, true);

% Check udmodlist
udmodlist = get(handles.modlist, 'userdata');
if isempty(udmodlist) || ~(~isempty(udmodlist.cjob) && cfg_util('isjob_id', udmodlist.cjob))
    udmodlist = local_init_udmodlist;
    udmodlist.cjob = cfg_util('initjob');
    set(handles.modlist, 'userdata', udmodlist);
end;

% set initial font
lf = cfg_get_defaults([mfilename '.lfont']);
local_setfont(hObject, lf);

% set ExpertEdit checkbox
set(handles.MenuViewExpertEdit, ...
    'checked',cfg_get_defaults([mfilename '.ExpertEdit']));

% show job
local_showjob(hObject);

% Choose default command line output for cfg_ui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Executes when cfg_ui is resized.
function cfg_ui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to cfg_ui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% this is just "Update View"
MenuViewUpdateView_Callback(hObject, eventdata, handles);

% --- Executes when user attempts to close cfg_ui.
function cfg_ui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to cfg_ui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
udmodlist = get(handles.modlist,'userdata');
cmd = local_check_job_modified(udmodlist, 'quit');
switch lower(cmd)
    case 'quit'
        if ~isempty(udmodlist.cjob)
            cfg_util('deljob', udmodlist.cjob);
        end;
        delete(hObject);
%         set(hObject,'Visible','off');
%         udmodlist = local_init_udmodlist;
%         udmodlist.cjob = cfg_util('initjob');
%         set(handles.modlist,'userdata',udmodlist);
    case 'hide'
        set(hObject,'Visible','off');
end;

% --- Outputs from this function are returned to the command line.
function varargout = cfg_ui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%% File Menu Callbacks
% --------------------------------------------------------------------
function MenuFileNew_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist, 'userdata');
cmd = local_check_job_modified(udmodlist, 'replace');
if strcmpi(cmd, 'continue')
    if ~isempty(udmodlist.cmod)
        cfg_util('deljob',udmodlist(1).cjob);
    end;
    udmodlist = local_init_udmodlist;
    udmodlist.cjob = cfg_util('initjob');
    set(handles.modlist, 'userdata', udmodlist);
    local_showjob(hObject);
end;

% --------------------------------------------------------------------
function MenuFileLoad_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist, 'userdata');
cmd = local_check_job_modified(udmodlist, 'replace');
if strcmpi(cmd, 'continue')
    [files, sts] = cfg_getfile([1 Inf], 'batch', 'Load Job File(s)', {}, udmodlist.wd);
    if sts
        local_pointer('watch');
        cfg_util('deljob',udmodlist(1).cjob);
        udmodlist = local_init_udmodlist;
        try
            udmodlist.wd = fileparts(files{1});
            udmodlist.cjob = cfg_util('initjob', files);
        catch
            l = lasterror;
            errordlg(l.message,'Error loading job', 'modal');
        end    
        set(handles.modlist, 'userdata', udmodlist);
        set(handles.module, 'userdata', []);
        local_showjob(hObject);
        local_pointer('arrow');
    end;
end;

% --------------------------------------------------------------------
function MenuFileSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist, 'userdata');
opwd = pwd;
if ~isempty(udmodlist.wd)
    try
        cd(udmodlist.wd);
    end
end;
[file, pth, idx] = uiputfile({'*.mat','Matlab .mat File';...
                    '*.m','Matlab .m Script File'}, 'Save Job');
cd(opwd);
if isnumeric(file) && file == 0
    return;
end;
local_pointer('watch');
[p, n, e] = fileparts(file);
if isempty(e) || ~any(strcmp(e,{'.mat','.m'}))
    e1 = {'.mat','.m'};
    e2 = e1{idx};
    file = sprintf('%s%s', n, e);
else
    file = n;
    e2 = e;
end
try
    cfg_util('savejob', udmodlist.cjob, fullfile(pth, [file e2]));
    udmodlist.modified = false;
    udmodlist.wd = pth;
    set(handles.modlist,'userdata',udmodlist);
catch
    l = lasterror;
    errordlg(l.message,'Error saving job', 'modal');
end
local_pointer('arrow');

% --------------------------------------------------------------------
function MenuFileScript_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileScript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------

udmodlist = get(handles.modlist, 'userdata');
opwd = pwd;
if ~isempty(udmodlist.wd)
    try
        cd(udmodlist.wd);
    end
end;
[file, pth, idx] = uiputfile({'*.m','Matlab .m Script File'},...
    'Script File name');
cd(opwd);
if isnumeric(file) && file == 0
    return;
end;
local_pointer('watch');
[p, n, e] = fileparts(file);
try
    cfg_util('genscript', udmodlist.cjob, pth, [n '.m']);
    udmodlist.modified = false;
    udmodlist.wd = pth;
    set(handles.modlist,'userdata',udmodlist);
    if ~isdeployed
        edit(fullfile(pth, [n '.m']));
    end
catch
    l = lasterror;
    errordlg(l.message,'Error generating job script', 'modal');
end
local_pointer('arrow');

% --------------------------------------------------------------------
function MenuFileRun_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_pointer('watch');
udmodlist = get(handles.modlist, 'userdata');
try
    cfg_util('run',udmodlist(1).cjob);
catch
    le = lasterror;
    if strcmpi(questdlg(sprintf('An error occured during job execution. Please see the MATLAB command window for details.\n\nSave error information?'),'Error in job execution', 'Yes','No','Yes'), 'yes')
        opwd = pwd;
        if ~isempty(udmodlist.wd)
            try
                cd(udmodlist.wd);
            end
        end;
        [file, pth, idx] = uiputfile({'*.mat','Matlab .mat File'},...
            'Error .mat File name');
        cd(opwd);
        if ~(isnumeric(file) && file == 0)
            [u1, ojob] = cfg_util('harvest',udmodlist(1).cjob);
            [u1, rjob] = cfg_util('harvestrun',udmodlist(1).cjob);
            outputs   = cfg_util('getalloutputs', udmodlist(1).cjob);
            diarystr  = cfg_util('getdiary',udmodlist(1).cjob);
            [p, n, e] = fileparts(file);
            save(fullfile(pth, [n '.mat']), 'ojob', 'rjob', 'outputs', 'diarystr');
        end;
    end
end;
local_pointer('arrow');

% --------------------------------------------------------------------
function MenuFileRunSerial_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileRunSerial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_pointer('watch');
udmodlist = get(handles.modlist, 'userdata');
cfg_util('runserial',udmodlist(1).cjob);
local_pointer('arrow');

% --------------------------------------------------------------------
function MenuFileAddApp_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileAddApp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist, 'userdata');
cont = local_check_job_modified(udmodlist, 'replace');
if cont
    [file, sts] = cfg_getfile([1 1], '.*\.m$', 'Load Application Configuration');
    if sts
        udmodlist = get(handles.modlist, 'userdata');
        if ~isempty(udmodlist.cmod)
            cfg_util('deljob',udmodlist(1).cjob);
        end;
        [p, fun, e] = fileparts(file{1});
        addpath(p);
        cfg_util('addapp', fun);
        local_setmenu(handles.cfg_ui, [], @local_AddMod, true);
        udmodlist = local_init_udmodlist;
        udmodlist.cjob = cfg_util('initjob');
        set(handles.modlist, 'userdata', udmodlist);
        local_showjob(hObject);
    end;
end;

% --------------------------------------------------------------------
function MenuFileClose_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.cfg_ui);

%%% Edit Menu Callbacks
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

local_clearvaledit(hObject);

% --------------------------------------------------------------------
function MenuEditValAddDep_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValAddDep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_dep(hObject);

%%% View Menu Callbacks
% --------------------------------------------------------------------
function MenuViewUpdateView_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewUpdateView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function seems to be called on startup without guidata - do nothing
% there
if ~isempty(handles) 
    local_setmenu(handles.cfg_ui, [], @local_AddMod, true);
    udmodlist = get(handles.modlist,'Userdata');
    if isstruct(udmodlist)
        if isfield(udmodlist,'defid')
            local_showmod(hObject);
        else
            local_showjob(hObject);
        end;
    end
end;

% --------------------------------------------------------------------
function MenuViewFontSize_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewFontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lf = cfg_get_defaults([mfilename '.lfont']);
fs = uisetfont(cell2struct(lf(2:2:end),lf(1:2:end),2));
if isstruct(fs)
    % construct argument list for set
    fn = fieldnames(fs);
    fs = struct2cell(fs);
    local_setfont(hObject,[fn'; fs']);
    MenuViewUpdateView_Callback(hObject, eventdata, handles);
end;

% --------------------------------------------------------------------
function MenuViewExpertEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewExpertEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(gcbo,'checked'),'on')
    newstate = 'off';
else
    newstate = 'on';
end;
set(gcbo, 'checked', newstate);
cfg_get_defaults([mfilename '.ExpertEdit'], newstate);

% --------------------------------------------------------------------
function MenuViewShowCode_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewShowCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist, 'userdata');
[un, matlabbatch] = cfg_util('harvest', udmodlist.cjob);
str = gencode(matlabbatch);
fg  = findobj(0,'Type','figure','Tag',[mfilename 'ShowCode']);
if isempty(fg)
    fg   = figure('Menubar','none', 'Toolbar','none', 'Tag',[mfilename 'ShowCode'], 'Units','normalized', 'Name','Batch Code Browser', 'NumberTitle','off');
    ctxt = uicontrol('Parent',fg, 'Style','listbox', 'Units','normalized', 'Position',[0 0 1 1], 'FontName','FixedWidth','Tag',[mfilename 'ShowCodeList']);
else
    figure(fg);
    ctxt = findobj(fg,'Tag',[mfilename 'ShowCodeList']); 
end
um = uicontextmenu;
um1 = uimenu('Label','Copy', 'Callback',@(ob,ev)local_ShowCode_Copy(ob,ev,ctxt), 'Parent',um);
um1 = uimenu('Label','Select all', 'Callback',@(ob,ev)local_ShowCode_SelAll(ob,ev,ctxt), 'Parent',um);
um1 = uimenu('Label','Unselect all', 'Callback',@(ob,ev)local_ShowCode_UnSelAll(ob,ev,ctxt), 'Parent',um);
set(ctxt, 'Min',-1, 'Max',numel(str), 'UIContextMenu',um, 'Value',[], 'ListboxTop',1);
set(ctxt, 'String',str);

%%% Module List Callbacks
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
        local_showjob(hObject);
    else
        local_showmod(hObject);
    end;
end;
% Return focus to modlist - otherwise it would be on current module
uicontrol(handles.modlist);

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

%%% Module List Context Menu Callbacks
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

%%% Module Callbacks
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
    doopen = false;
else
    doopen = strcmp(get(handles.cfg_ui,'SelectionType'),'open');
    % When the user hits SPACE, the listbox callback will be executed twice
    % out of order. Once, it will be called to open an item for editing.
    % Meanwhile, a second call will advance to the next line in the module
    % list. To make sure we show the right data, we try to restore the
    % original value setting from the 'open' call.
    set(hObject, 'Value', value);
end
local_showvaledit(hObject);
if doopen
    % open modal MenuEdit window, do editing
    % Unfortunately, MATLAB focus behaviour makes it impossible to do this
    % in the existing valshow object - if this object looses focus, it will
    % call its callback without checking why it lost focus.
    % Call appropriate input handler for editable and selection types
    local_valedit_EditValue(hObject);
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

%%% Value Display Callbacks
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

%%% GUI Buttons
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

local_valedit_dep(hObject);

% --------------------------------------------------------------------
function CmValEditVal_Callback(hObject, eventdata, handles)
% hObject    handle to CmValEditVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_EditValue(hObject);

% --------------------------------------------------------------------
function CmValAddDep_Callback(hObject, eventdata, handles)
% hObject    handle to CmValAddDep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_dep(hObject);

% --------------------------------------------------------------------
function CmValClearVal_Callback(hObject, eventdata, handles)
% hObject    handle to CmValClearVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_clearvaledit(hObject);


% --------------------------------------------------------------------
function MenuFileMultiBatch_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileMultiBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist,'userdata');
cjob = udmodlist.cjob;
if cfg_util('isjob_id', cjob) 
    njob = cfg_util('clonejob', cjob);
    cfg_ui_multibatch(njob);
end


% --------------------------------------------------------------------
function MenuViewPreview_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_preview(hObject);

% --------------------------------------------------------------------
function CmValPreview_Callback(hObject, eventdata, handles)
% hObject    handle to CmValPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_preview(hObject);
