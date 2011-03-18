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
% $Id: cfg_ui.m 4255 2011-03-18 13:11:03Z volkmar $

rev = '$Rev: 4255 $'; %#ok

% edit the above text to modify the response to help cfg_ui

% Last Modified by GUIDE v2.5 04-Mar-2010 16:09:52

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

%% Input evaluation
% --------------------------------------------------------------------
function [val, sts] = local_eval_valedit(varargin)
% for security reasons, separate string evaluation from valedit_Callback
% (evaluation might overwrite variables). Uses evalc to suppress any
% console output.
val = [];
sts = false;
try
    % 1st, try to convert into numeric matrix without evaluation
    % This converts expressions like '1 -1' into [1 -1] instead of
    % evaluating them
    [val sts] = str2num(varargin{1}); %#ok<ST2NM>
    if ~sts
        % try to evaluate str as rvalue
        val = evalin('base', varargin{1});
        sts = true;
    end
catch
    everr = lasterror;
    if strcmp(everr.identifier, 'MATLAB:m_invalid_lhs_of_assignment')
        try
            evalin('base', varargin{1});
            val = evalin('base','val');
            % test if val variable exists
            if ~exist('val','var')
                cfg_message('cfg_ui:local_eval_valedit:noval','No variable ''val'' assigned.');
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
function local_addtojob(varargin)
id  = get(gcbo, 'userdata');
handles = guidata(gcbo);
udmodlist = get(handles.modlist, 'userdata');
% add module to job, harvest to initialise its virtual outputs
mod_job_id = cfg_util('addtojob', udmodlist.cjob, id);
cfg_util('harvest', udmodlist.cjob, mod_job_id);
udmodlist.modified = true;
set(handles.modlist,'userdata',udmodlist);
local_showjob(gcbo);
% --------------------------------------------------------------------
function local_loaddefs(varargin)
appid = get(gcbo, 'Userdata');
[file sts] = cfg_getfile(1, '.*\.m$','Load Defaults from');
if sts
    cfg_util('initdef', appid, file{1});
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
fid = fopen(fullfile(path, file), 'wt');
if fid < 1
    cfg_message('matlabbatch:savefailed', ...
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
    udmodlist = local_init_udmodlist;
    udmodlist(1).cjob = cjob;
    % move figure onscreen
    cfg_onscreen(obj);
    set(obj,'Visible','on');
end;
[id str sts dep sout] = cfg_util('showjob',cjob);
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
ltop = local_getListboxTop(handles.modlist, cmod, numel(str));
set(handles.modlist, 'userdata',udmodlist, 'value', cmod, 'ListboxTop', ltop, 'string', str);
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
                  'all_set','all_set_item','num'});
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
    set(handles.moduleHead,'String',sprintf('Current Module: %s', contents{1}{1}));
    namestr = cell(1,numel(contents{1}));
    datastr = cell(1,numel(contents{1}));
    namestr{1} = sprintf('Help on: %s',contents{1}{1});
    datastr{1} = '';
    for k = 2:numel(contents{1})
        if contents{6}{k}-2 > 0
            indent = [' ' repmat('. ', 1, contents{6}{k}-2)];
        else
            indent = '';
        end
        if contents{8}{k} || (isfield(udmodlist,'defid') && ~isempty(contents{2}{k}))
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
                        datastr{k} = 'Unknown selection';
                        for l = 1:numel(contents{4}{k})
                            if isequalwithequalnans(contents{2}{k}{1}, contents{4}{k}{l})
                                datastr{k} = contents{3}{k}{l};
                                break;
                            end;
                        end;
                    case 'cfg_files',
                        if numel(contents{2}{k}{1}) == 1
                            if isempty(contents{2}{k}{1}{1})
                                datastr{k} = ' ';
                            else
                                datastr{k} = contents{2}{k}{1}{1};
                            end;
                        else
                            datastr{k} = sprintf('%d files', numel(contents{2}{k}{1}));
                        end;
                    case 'cfg_entry'
                        csz = size(contents{2}{k}{1});
                        % TODO use gencode like string formatting
                        if ischar(contents{2}{k}{1}) && ...
                                numel(csz) == 2 && any(csz(1:2) == 1)
                            datastr{k} = contents{2}{k}{1};
                        elseif (isnumeric(contents{2}{k}{1}) || ...
                                islogical(contents{2}{k}{1})) && ...
                                numel(csz) == 2 && any(csz(1:2) == 1) &&...
                                numel(contents{2}{k}{1}) <= 4
                            % always display line vector as summary
                            datastr{k} = mat2str(contents{2}{k}{1}(:)');
                        elseif any(csz == 0)
                            switch class(contents{2}{k}{1})
                                case 'char',
                                    datastr{k} = '''''';
                                case 'double',
                                    datastr{k} = '[]';
                                otherwise
                                    datastr{k} = sprintf('%s([])', ...
                                                         class(contents{2}{k}{1}));
                            end;
                        else
                            szstr = sprintf('%dx', csz);
                            datastr{k} = sprintf('%s %s', ...
                                                 szstr(1:end-1), class(contents{2}{k}{1}));
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
    ltop = local_getListboxTop(handles.module, citem, numel(str));
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
    [id stop help] = cfg_util('listcfgall', [], cfg_findspec({{'tag','matlabbatch'}}), {'showdoc'});
    set(handles.helpbox, 'Value',1, 'ListboxTop',1, 'String',cfg_justify(handles.helpbox, help{1}{1}));
end;

%% Show Item
% --------------------------------------------------------------------
function local_showvaledit(obj)
handles = guidata(obj);
udmodlist = get(handles.modlist, 'userdata');
cmod = get(handles.modlist, 'value');
udmodule = get(handles.module, 'userdata');
citem = get(handles.module, 'value');
set(findobj(handles.cfg_ui,'-regexp', 'Tag','^BtnVal.*'), 'Visible','off');
set(findobj(handles.cfg_ui,'-regexp', 'Tag','^MenuEditVal.*'), 'Enable','off');
set(findobj(handles.cfg_ui,'-regexp', 'Tag','^CmVal.*'), 'Visible','off');
delete(findobj(handles.cfg_ui,'-regexp', 'Tag','^MenuEditVal.*Dyn$'))
set(findobj(handles.cfg_ui,'-regexp', 'Tag','^valshow.*'), 'Visible','off');
set(handles.valshow,'String', '','Min',0,'Max',0,'Callback',[]);
set(handles.valshowLabel, 'String',sprintf('Current Item: %s',udmodule.contents{1}{citem}));
switch(udmodule.contents{5}{citem})
    case {'cfg_entry','cfg_files'}
        if ~isempty(udmodule.contents{2}{citem}) && isa(udmodule.contents{2}{citem}{1}, 'cfg_dep')
            str = {'Reference from'};
            for k = 1:numel(udmodule.contents{2}{citem}{1}) % we may have multiple dependencies
                str{k+1} = udmodule.contents{2}{citem}{1}(k).sname; % return something to be printed
            end;
        elseif ~isempty(udmodule.contents{2}{citem})
            if ndims(udmodule.contents{2}{citem}{1}) <= 2 
                if ischar(udmodule.contents{2}{citem}{1})
                    str = cellstr(udmodule.contents{2}{citem}{1});
                elseif iscellstr(udmodule.contents{2}{citem}{1})
                    str = udmodule.contents{2}{citem}{1};
                elseif isnumeric(udmodule.contents{2}{citem}{1}) || ...
                        islogical(udmodule.contents{2}{citem}{1})
                    str = cellstr(num2str(udmodule.contents{2}{citem}{1}));
                else
                    str = gencode(udmodule.contents{2}{citem}{1},'val');
                end;
            else
               str = gencode(udmodule.contents{2}{citem}{1},'val');
            end;
        else
            str = '';
        end;
        set(handles.valshow, 'Visible','on', 'Value',1, 'ListboxTop',1,'String', str);
        set(handles.valshowLabel, 'Visible','on');
        if ~isfield(udmodlist, 'defid')
            sout = local_showvaledit_deps(obj);
            if ~isempty(sout)
                set(findobj(handles.cfg_ui,'-regexp','Tag','.*AddDep$'), ...
                    'Visible','on', 'Enable','on');
            end;
        end;
        if strcmp(udmodule.contents{5}{citem},'cfg_files')
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*SelectFiles$'), ...
                'Visible','on', 'Enable','on');
        else
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*EditVal$'), ...
                'Visible','on', 'Enable','on');
        end
        set(findobj(handles.cfg_ui,'-regexp','Tag','.*ClearVal$'), ...
            'Visible','on', 'Enable','on');
    case {'cfg_menu','cfg_choice'}
        if strcmp(udmodule.contents{5}{citem},'cfg_menu') || ~isfield(udmodlist, 'defid')
            cval = -1;
            if strcmp(udmodule.contents{5}{citem},'cfg_choice')
                % compare tag, not filled entries
                cmpsubs = substruct('.','tag');
            else
                cmpsubs = struct('type',{},'subs',{});
            end;
            valsubs = substruct('{}',{1});
            nitem = numel(udmodule.contents{4}{citem});
            mrk = cell(1,nitem);
            for l = 1:nitem
                valuesubs = substruct('{}',{l});
                if ~isempty(udmodule.contents{2}{citem}) && isequal(subsref(udmodule.contents{2}{citem},[valsubs cmpsubs]), subsref(udmodule.contents{4}{citem},[valuesubs cmpsubs]))
                    mrk{l} = '*';
                    cval = l;
                else
                    mrk{l} = ' ';
                end;
            end;
            if strcmp(udmodule.contents{5}{citem},'cfg_choice')
                str = cell(1,nitem);
                for k = 1:nitem
                    str{k} = udmodule.contents{4}{citem}{k}.name;
                end;
            else
                str = udmodule.contents{3}{citem};
            end;
            str = strcat(mrk(:), str(:));
            udvalshow = local_init_udvalshow;
            udvalshow.cval = cval;
            if cval == -1
                cval = 1;
            end;
            ltop = local_getListboxTop(handles.valshow, cval, numel(str));
            set(handles.valshow, 'Visible','on', 'Value',cval, 'ListboxTop',ltop, 'String',str, ...
                              'Callback',@local_valedit_list, ...
                              'Keypressfcn',@local_valedit_key, ...
                              'Userdata',udvalshow);
            set(handles.valshowLabel, 'Visible','on');
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*EditVal$'), ...
                'Visible','on', 'Enable','on');
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*ClearVal$'), ...
                'Visible','on', 'Enable','on');
        end;
    case {'cfg_repeat'}
        if ~isfield(udmodlist, 'defid')
            udvalshow = local_init_udvalshow;
            udvalshow.cval = 1;
            % Already selected items
            ncitems = numel(udmodule.contents{2}{citem});
            str3 = cell(ncitems,1);
            cmd3 = cell(ncitems,1);
            for k = 1:ncitems
                str3{k} = sprintf('Delete: %s (%d)',...
                                          udmodule.contents{2}{citem}{k}.name, k);
                cmd3{k} = [Inf k];
                uimenu(handles.MenuEditValDelItem, ...
                       'Label',sprintf('%s (%d)', ...
                                       udmodule.contents{2}{citem}{k}.name, k), ...
                       'Callback',@(ob,ev)local_setvaledit(ob,cmd3{k},ev), ...
                       'Tag','MenuEditValDelItemDyn');
                uimenu(handles.CmValDelItem, ...
                       'Label',sprintf('%s (%d)', ...
                                       udmodule.contents{2}{citem}{k}.name, k), ...
                       'Callback',@(ob,ev)local_setvaledit(ob,cmd3{k},ev), ...
                       'Tag','CmValDelItemDyn');
            end
            % Add/Replicate callbacks will be shown only if max number of
            % items not yet reached
            if ncitems < udmodule.contents{9}{citem}(2)
                % Available items
                naitems = numel(udmodule.contents{4}{citem});
                str1 = cell(naitems,1);
                cmd1 = cell(naitems,1);
                for k = 1:naitems
                    str1{k} = sprintf('New: %s', udmodule.contents{4}{citem}{k}.name);
                    cmd1{k} = [k Inf];
                    uimenu(handles.MenuEditValAddItem, ...
                           'Label',udmodule.contents{4}{citem}{k}.name, ...
                           'Callback',@(ob,ev)local_setvaledit(ob,cmd1{k},ev), ...
                           'Tag','MenuEditValAddItemDyn');
                    uimenu(handles.CmValAddItem, ...
                           'Label',udmodule.contents{4}{citem}{k}.name, ...
                           'Callback',@(ob,ev)local_setvaledit(ob,cmd1{k},ev), ...
                           'Tag','CmValAddItemDyn');
                end;
                str2 = cell(ncitems,1);
                cmd2 = cell(ncitems,1);
                for k = 1:ncitems
                    str2{k} = sprintf('Replicate: %s (%d)',...
                                      udmodule.contents{2}{citem}{k}.name, k);
                    cmd2{k} = [-1 k];
                    uimenu(handles.MenuEditValReplItem, ...
                           'Label',sprintf('%s (%d)', ...
                                           udmodule.contents{2}{citem}{k}.name, k), ...
                           'Callback',@(ob,ev)local_setvaledit(ob,cmd2{k},ev), ...
                           'Tag','MenuEditValReplItemDyn');
                    uimenu(handles.CmValReplItem, ...
                           'Label',sprintf('%s (%d)', ...
                                           udmodule.contents{2}{citem}{k}.name, k), ...
                           'Callback',@(ob,ev)local_setvaledit(ob,cmd2{k},ev), ...
                           'Tag','CmValReplItemDyn');
                end
                set(findobj(handles.cfg_ui,'-regexp','Tag','.*AddItem$'), ...
                    'Visible','on', 'Enable','on');
                if ncitems > 0
                    set(findobj(handles.cfg_ui,'-regexp','Tag','.*ReplItem$'), ...
                        'Visible','on', 'Enable','on');
                    set(findobj(handles.cfg_ui,'-regexp','Tag','.*ReplItem$'), ...
                        'Visible','on', 'Enable','on');
                end
            else
                str1 = {};
                str2 = {};
                cmd1 = {};
                cmd2 = {};
            end
            str = [str1(:); str2(:); str3(:)];
            udvalshow.cmd = [cmd1(:); cmd2(:); cmd3(:)];
            set(handles.valshow, 'Visible','on', 'Value',1, 'ListboxTop',1, 'String', str, ...
                'Callback',@local_valedit_repeat, ...
                'KeyPressFcn', @local_valedit_key, ...
                'Userdata',udvalshow);
            set(handles.valshowLabel, 'Visible','on');
            set(findobj(handles.cfg_ui,'-regexp','Tag','^Btn.*EditVal$'), ...
                'Visible','on', 'Enable','on');
            if ncitems > 0
                set(findobj(handles.cfg_ui,'-regexp','Tag','.*DelItem$'), ...
                    'Visible','on', 'Enable','on');
            end;
            set(findobj(handles.cfg_ui,'-regexp','Tag','.*ClearVal$'), ...
                'Visible','on', 'Enable','on');
        end;
end;
if isfield(udmodlist, 'defid')
    cmid = udmodlist.defid(cmod);
else
    cmid = {udmodlist.cjob udmodlist.id{cmod}};
end;
[id stop help] = cfg_util('listmod', cmid{:}, udmodule.id{citem}, cfg_findspec, ...
                          cfg_tropts(cfg_findspec,1,1,1,1,false), {'showdoc'});
set(handles.helpbox, 'Value',1, 'ListboxTop',1, 'string',cfg_justify(handles.helpbox, help{1}{1}));
drawnow;

%% List matching dependencies
% --------------------------------------------------------------------
function sout = local_showvaledit_deps(obj)
handles = guidata(obj);
udmodlist = get(handles.modlist, 'userdata');
cmod = get(handles.modlist, 'value');
udmodule = get(handles.module, 'userdata');
citem = get(handles.module, 'value');
sout = cat(2,udmodlist(1).sout{1:cmod-1});
smatch = false(size(sout));
% loop over sout to find whether there are dependencies that match the current item
for k = 1:numel(sout)
    smatch(k) = cfg_util('match', udmodlist.cjob, udmodlist.id{cmod}, udmodule.id{citem}, sout(k).tgt_spec);
end;
sout = sout(smatch);

%% Value edit dialogue
% --------------------------------------------------------------------

function local_valedit_edit(hObject)
% Normal mode. Depending on strtype, put '' or [] around entered
% input. If input has ndims > 2, isn't numeric or char, proceed with
% expert dialogue.
handles = guidata(hObject);
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
cmod = get(handles.modlist, 'Value');
udmodlist = get(handles.modlist, 'Userdata');
val = udmodule.contents{2}{value};
if isfield(udmodlist, 'defid')
    cmid = udmodlist.defid(cmod);
else
    cmid = {udmodlist.cjob, udmodlist.id{cmod}};
end;
[id stop strtype] = cfg_util('listmod', cmid{:}, udmodule.id{value}, cfg_findspec, ...
                             cfg_tropts(cfg_findspec,1,1,1,1,false), {'strtype'});
if isempty(val) || isa(val{1}, 'cfg_dep')
    % silently clear cfg_deps
    if strtype{1}{1} == 's'
        val = {''};
    else
        val = {[]};
    end;
end;
% If requested or we can't handle this, use expert mode
expmode = strcmp(cfg_get_defaults([mfilename '.ExpertEdit']), 'on') ||...
    ndims(val{1}) > 2 || ~(ischar(val{1}) || isnumeric(val{1}) || islogical(val{1}));
% Generate code for current value, if not empty
% Set dialog texts
if expmode
    if ~isequal(val, {''})
        instr = gencode(val{1},'val');
        % remove comments and put in 1-cell multiline char array
        nc = cellfun(@isempty,regexp(instr,'^\s*%'));
        instr = {char(instr(nc))};
    else
        instr = {''};
    end
    hlptxt = char({'Enter a valid MATLAB expression.', ...
        ' ', ...
        ['Strings must be enclosed in single quotes ' ...
        '(''A''), multiline arrays in brackets ([ ]).'], ...
        ' ', ...
        'To clear a value, enter an empty cell ''{}''.', ...
        ' ', ...
        ['Accept input with CTRL-RETURN, cancel with ' ...
        'ESC.']});
    failtxt = {'Input could not be evaluated. Possible reasons are:',...
        '1) Input should be a vector or matrix, but is not enclosed in ''['' and '']'' brackets.',...
        '2) Input should be a character or string, but is not enclosed in '' single quotes.',...
        '3) Input should be a MATLAB variable, but is misspelled.',...
        '4) Input should be a MATLAB expression, but has syntax errors.'};
else
    if strtype{1}{1} == 's'
        instr = val;
        encl  = {'''' ''''};
    else
        try
            instr = {num2str(val{1})};
        catch
            instr = {''};
        end;
        encl  = {'[' ']'};
    end;
    hlptxt = char({'Enter a value.', ...
        ' ', ...
        'To clear a value, clear the input field and accept.', ...
        ' ', ...
        ['Accept input with CTRL-RETURN, cancel with ' ...
        'ESC.']});
    failtxt = {'Input could not be evaluated.'};
end
sts = false;
while ~sts
    % estimate size of input field based on instr
    % Maximum widthxheight 140x20, minimum 60x2
    szi = size(instr{1});
    mxwidth = 140;
    rdup = ceil(szi(2)/mxwidth)+3;
    szi = max(min([szi(1)*rdup szi(2)],[20 140]),[2,60]);
    str = inputdlg(hlptxt, ...
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
    % context - graphics handles are made invisible to avoid accidental
    % damage
    hv = local_disable(handles.cfg_ui,'HandleVisibility');
    [val sts] = local_eval_valedit(str);
    local_enable(handles.cfg_ui,'HandleVisibility',hv);
    % for strtype 's', val must be a string
    sts = sts && (~strcmp(strtype{1}{1},'s') || ischar(val));
    if ~sts
        if ~expmode
            % try with matching value enclosure
            if strtype{1}{1} == 's'
                if ishandle(val) % delete accidentally created objects
                    delete(val);
                end
                % escape single quotes and place the whole string in single quotes
                str = strcat(encl(1), strrep(cstr,'''',''''''), encl(2), {char(10)});
            else
                cestr = [encl(1); cstr(:); encl(2)]';
                str = strcat(cestr, {char(10)});
            end;
            str = cat(2, str{:});
            % Evaluation is encapsulated to avoid users compromising this function
            % context - graphics handles are made invisible to avoid accidental
            % damage
            hv = local_disable(handles.cfg_ui,'HandleVisibility');
            [val sts] = local_eval_valedit(str);
            local_enable(handles.cfg_ui,'HandleVisibility',hv);
        end;
        if ~sts % (Still) no valid input
            uiwait(msgbox(failtxt,'Evaluation error','modal'));
        end;
    end;
end
% This code will only be reached if a new value has been set
local_setvaledit(hObject, val);

% --------------------------------------------------------------------
function en = local_disable(hObject, property)
% disable property in all ui objects, returning their previous state in
% cell list en
handles = guidata(hObject);
c = findall(handles.cfg_ui);
en = cell(size(c));
sel = isprop(c,property);
en(sel) = get(c(sel),property);
set(c(sel),property,'off');

% --------------------------------------------------------------------
function local_enable(hObject, property, en)
% reset original property status. if en is empty, return without changing
% anything. 
if ~isempty(en)
    handles = guidata(hObject);
    c = findall(handles.cfg_ui);
    sel = isprop(c,property);
    set(c(sel),{property},en(sel));
end

%% Value choice dialogue
% --------------------------------------------------------------------
function local_valedit_accept(hObject,val)
handles = guidata(hObject);
udvalshow = get(handles.valshow, 'Userdata');
local_enable(hObject, 'Enable', udvalshow.en);
uiresume(handles.cfg_ui);
set(findobj(handles.cfg_ui, '-regexp', 'Tag','^valshowBtn.*'), ...
    'Visible','off', 'Callback',[]);
if nargin == 1
    val = get(handles.valshow, 'Value');
end
local_setvaledit(hObject,val);

% --------------------------------------------------------------------
function local_valedit_cancel(hObject)
handles = guidata(hObject);
udvalshow = get(handles.valshow, 'Userdata');
local_enable(hObject, 'Enable', udvalshow.en);
uiresume(handles.cfg_ui);
set(findobj(handles.cfg_ui, '-regexp', 'Tag','^valshowBtn.*'), ...
    'Visible','off', 'Callback',[]);
uicontrol(handles.module);

% --------------------------------------------------------------------
function local_valedit_key(hObject, data, varargin)
if strcmpi(data.Key,'escape')
    % ESC must be checked here
    local_valedit_cancel(hObject);
else
    % collect key info for evaluation in local_valedit_list
    handles = guidata(hObject);
    udvalshow = get(handles.valshow, 'Userdata');
    udvalshow.key = data;
    set(handles.valshow, 'Userdata',udvalshow);
end

% --------------------------------------------------------------------
function local_valedit_list(hObject,varargin)
handles = guidata(hObject);
udvalshow = get(handles.valshow, 'Userdata');
val = get(handles.valshow, 'Value');
if ((isempty(udvalshow.key) || ...
        strcmpi(udvalshow.key.Key,'return')) && ...
        isequal(hObject, handles.valshow)) || ...
        isequal(hObject, handles.valshowBtnAccept)
    % callback called from handles.valshow, finish editing and set value
    if val ~= udvalshow.cval
        local_valedit_accept(hObject);
    else
        local_valedit_cancel(hObject);
    end;
elseif  (~isempty(udvalshow.key) && strcmpi(udvalshow.key.Key,'escape') && ...
        isequal(hObject, handles.valshow)) || ...
        isequal(hObject, handles.valshowBtnCancel)
    local_valedit_cancel(hObject);
elseif ~isequal(hObject, handles.valshow)
    % callback called from elsewhere (module, menu, button) - init editing
    udvalshow.en = local_disable(hObject,'Enable');
    udvalshow.key = [];
    figure(handles.cfg_ui);
    set(handles.valshow,'enable','on','Userdata',udvalshow, 'Min',0, ...
                      'Max',1);
    set(findobj(handles.cfg_ui, '-regexp', 'Tag','^valshowBtn.*'), ...
        'Enable','on', 'Visible','on', 'Callback',@local_valedit_list);
    uicontrol(handles.valshow);
    uiwait(handles.cfg_ui);
else
    udvalshow.key = [];
    set(handles.valshow, 'Userdata',udvalshow);
end;

%% Repeat edit dialogue
% --------------------------------------------------------------------
function local_valedit_repeat(hObject,varargin)
handles = guidata(hObject);
udvalshow = get(handles.valshow, 'Userdata');
if ((isempty(udvalshow.key) || ...
        strcmpi(udvalshow.key.Key,'return')) && ...
        isequal(hObject, handles.valshow)) || ...
        isequal(hObject, handles.valshowBtnAccept)
    % Mouse selection - no key
    % Keyboard selection - return key
    ccmd = get(hObject,'Value');
    local_valedit_accept(hObject,udvalshow.cmd{ccmd});
elseif (~isempty(udvalshow.key) && strcmpi(udvalshow.key.Key,'escape') && ...
        isequal(hObject, handles.valshow)) || ...
        isequal(hObject, handles.valshowBtnCancel)
    % callback called from handles.valshow, finish editing
    local_valedit_cancel(hObject);
elseif ~isequal(hObject, handles.valshow)
    % callback called from elsewhere (module, menu, button)
    % if there is just one action, do it, else init editing
    if numel(udvalshow.cmd) == 1
        local_valedit_accept(hObject,udvalshow.cmd{1});
    else
        udvalshow.en = local_disable(hObject,'Enable');
        udvalshow.key = [];
        figure(handles.cfg_ui);
        set(handles.valshow,'enable','on','Userdata',udvalshow, 'Min',0, ...
            'Max',1);
        set(findobj(handles.cfg_ui, '-regexp', 'Tag','^valshowBtn.*'), ...
            'Enable','on', 'Visible','on', 'Callback',@local_valedit_repeat);
        uicontrol(handles.valshow);
        uiwait(handles.cfg_ui);
    end
else
    udvalshow.key = [];
    set(handles.valshow, 'Userdata',udvalshow);
end;

%% Set changed values
% --------------------------------------------------------------------
function local_setvaledit(obj, val,varargin)
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
    udmodlist.modified = true;
    set(handles.modlist,'userdata',udmodlist);
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
    cmid = udmodlist.defid(cmod);
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
    local_setvaledit(hObject, val);
end;

% --------------------------------------------------------------------
function local_valedit_EditValue(hObject)
handles = guidata(hObject);
value = get(handles.module, 'Value');
udmodule = get(handles.module, 'Userdata');
switch udmodule.contents{5}{value}
    case {'cfg_choice','cfg_menu'}
        local_valedit_list(hObject);
    case 'cfg_repeat'
        local_valedit_repeat(hObject);
    otherwise
        local_valedit_edit(hObject);
end;

% --------------------------------------------------------------------
function local_valedit_ClearVal(hObject)
local_setvaledit(hObject, {});

% --------------------------------------------------------------------
function local_valedit_AddDep(hObject)
handles = guidata(hObject);
udmodule = get(handles.module, 'Userdata');
citem = get(handles.module, 'Value');
sout = local_showvaledit_deps(hObject);
str = {sout.sname};
[val sts] = listdlg('Name',udmodule.contents{1}{citem}, 'ListString',str);
if sts
    local_setvaledit(hObject, sout(val));
end;

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
local_setmenu(handles.cfg_ui, [], @local_addtojob, true);

% Check udmodlist
udmodlist = get(handles.modlist, 'userdata');
if isempty(udmodlist) || ~(~isempty(udmodlist.cjob) && cfg_util('isjob_id', udmodlist.cjob))
    udmodlist = local_init_udmodlist;
    udmodlist.cjob = cfg_util('initjob');
    set(handles.modlist, 'userdata', udmodlist);
end;

% set initial font
fs = cfg_get_defaults([mfilename '.lfont']);
local_setfont(hObject, fs);

% set ExpertEdit checkbox
set(handles.MenuViewExpertEdit, ...
    'checked',cfg_get_defaults([mfilename '.ExpertEdit']));

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
if isfield(udmodlist, 'modified') && udmodlist.modified
        cmd = questdlg(['The current batch contains unsaved changes. '...
            'Do you want to replace it with another batch?'], ...
            'Unsaved Changes', 'Continue','Cancel', 'Continue');
else
    cmd = 'Continue';
end;
if strcmpi(cmd,'continue')
    udmodlist = get(handles.modlist, 'userdata');
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
if udmodlist.modified
        cmd = questdlg(['The current batch contains unsaved changes. '...
            'Do you want to replace it with another batch?'], ...
            'Unsaved Changes', 'Continue','Cancel', 'Continue');
else
    cmd = 'Continue';
end;
if strcmpi(cmd,'continue')
    [files sts] = cfg_getfile([1 Inf], 'batch', 'Load Job File(s)', {}, udmodlist.wd);
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
    cd(udmodlist.wd);
end;
[file pth idx] = uiputfile({'*.mat','Matlab .mat File';...
                    '*.m','Matlab .m Script File'}, 'Save Job');
cd(opwd);
if isnumeric(file) && file == 0
    return;
end;
local_pointer('watch');
[p n e] = fileparts(file);
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
    cd(udmodlist.wd);
end;
[file pth idx] = uiputfile({'*.m','Matlab .m Script File'},...
    'Script File name');
cd(opwd);
if isnumeric(file) && file == 0
    return;
end;
local_pointer('watch');
[p n e] = fileparts(file);
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

function MenuFileRun_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_pointer('watch');
udmodlist = get(handles.modlist, 'userdata');
try
    cfg_util('run',udmodlist(1).cjob);
catch
    l = lasterror;
    errordlg(l.message,'Error in job execution', 'modal');
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
if udmodlist.modified
        cmd = questdlg(['The current batch contains unsaved changes. '...
            'Adding a new application will discard this batch.'], ...
            'Unsaved Changes', 'Continue','Cancel', 'Continue');
else
    cmd = 'Continue';
end;
if strcmpi(cmd,'continue')
    [file sts] = cfg_getfile([1 1], '.*\.m$', 'Load Application Configuration');
    if sts
        udmodlist = get(handles.modlist, 'userdata');
        if ~isempty(udmodlist.cmod)
            cfg_util('deljob',udmodlist(1).cjob);
        end;
        [p fun e] = fileparts(file{1});
        addpath(p);
        cfg_util('addapp', fun);
        local_setmenu(handles.cfg_ui, [], @local_addtojob, true);
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

%% Edit Menu Callbacks
% --------------------------------------------------------------------
function MenuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MenuViewUpdateView_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewUpdateView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function seems to be called on startup without guidata - do nothing
% there
if ~isempty(handles) 
    local_setmenu(handles.cfg_ui, [], @local_addtojob, true);
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
        case {'cfg_repeat'},
            local_valedit_repeat(hObject);
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

% --------------------------------------------------------------------
function udmodlist = local_init_udmodlist
% Initialise udmodlist to empty struct
% Don't initialise defid field - this will be added by defaults editor
udmodlist = struct('cjob',[],'cmod',[],'id',[],'sout',[],'modified',false,'wd','');

% --------------------------------------------------------------------
function udvalshow = local_init_udvalshow
% Initialise udvalshow to empty struct
udvalshow = struct('cval',[],'en',[],'key',[]);

% --------------------------------------------------------------------
function ltop = local_getListboxTop(obj, val, maxval)
% Get a safe value for ListboxTop property while keeping previous settings
% if possible.
% obj     handle of Listbox object
% val     new Value property
% maxval  new number of lines in obj
oltop = get(obj, 'ListboxTop');
ltop  = min([max(oltop,1), max(val-1,1), maxval]);

% --- Executes when user attempts to close cfg_ui.
function cfg_ui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to cfg_ui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
udmodlist = get(handles.modlist,'userdata');
if udmodlist.modified
    cmd = questdlg(['The current batch contains unsaved changes. Do you want to quit ' ...
                    'anyway or do you want to hide the batch window ' ...
                    'instead?'], 'Unsaved Changes', 'Quit','Cancel','Hide', ...
                   'Quit');
else
    cmd = 'Quit';
end;
switch lower(cmd)
    case 'quit'
        if ~isempty(udmodlist.cjob)
            cfg_util('deljob', udmodlist.cjob);
        end;
        set(hObject,'Visible','off');
        udmodlist = local_init_udmodlist;
        udmodlist.cjob = cfg_util('initjob');
        set(handles.modlist,'userdata',udmodlist);
    case 'hide'
        set(hObject,'Visible','off');
end;


% --- Executes on button press in valshowBtnAccept.
function valshowBtnAccept_Callback(hObject, eventdata, handles)
% hObject    handle to valshowBtnAccept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in valshowBtnCancel.
function valshowBtnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to valshowBtnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuViewFontSize_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewFontSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fs = uisetfont(cfg_get_defaults([mfilename '.lfont']));
if isstruct(fs)
    local_setfont(hObject,fs);
    MenuViewUpdateView_Callback(hObject, eventdata, handles);
end;

% --- Executes when cfg_ui is resized.
function cfg_ui_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to cfg_ui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% this is just "Update View"
MenuViewUpdateView_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function local_setfont(obj,fs)
handles = guidata(obj);
cfg_get_defaults([mfilename '.lfont'], fs);
% construct argument list for set
fn = fieldnames(fs);
fs = struct2cell(fs);
fnfs = [fn'; fs'];
set(handles.modlist, fnfs{:});
set(handles.module, fnfs{:});
set(handles.valshow, fnfs{:});
set(handles.helpbox, fnfs{:});

% --------------------------------------------------------------------
function local_pointer(ptr)
shh = get(0,'showhiddenhandles');
set(0,'showhiddenhandles','on');
set(get(0,'Children'),'Pointer',ptr);
drawnow;
set(0,'showhiddenhandles',shh);


% --------------------------------------------------------------------
function CmValEditVal_Callback(hObject, eventdata, handles)
% hObject    handle to CmValEditVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_EditValue(hObject);

% --------------------------------------------------------------------
function CmValSelectFiles_Callback(hObject, eventdata, handles)
% hObject    handle to CmValSelectFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_SelectFiles(hObject);

% --------------------------------------------------------------------
function CmValAddDep_Callback(hObject, eventdata, handles)
% hObject    handle to CmValAddDep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_AddDep(hObject);

% --------------------------------------------------------------------
function CmValClearVal_Callback(hObject, eventdata, handles)
% hObject    handle to CmValClearVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

local_valedit_ClearVal(hObject);


% --------------------------------------------------------------------
function MenuView_Callback(hObject, eventdata, handles)
% hObject    handle to MenuView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuViewShowCode_Callback(hObject, eventdata, handles)
% hObject    handle to MenuViewShowCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

udmodlist = get(handles.modlist, 'userdata');
[un matlabbatch] = cfg_util('harvest', udmodlist.cjob);
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
um1 = uimenu('Label','Copy', 'Callback',@(ob,ev)ShowCode_Copy(ob,ev,ctxt), 'Parent',um);
um1 = uimenu('Label','Select all', 'Callback',@(ob,ev)ShowCode_SelAll(ob,ev,ctxt), 'Parent',um);
um1 = uimenu('Label','Unselect all', 'Callback',@(ob,ev)ShowCode_UnSelAll(ob,ev,ctxt), 'Parent',um);
set(ctxt, 'Max',numel(str), 'UIContextMenu',um, 'Value',[], 'ListboxTop',1);
set(ctxt, 'String',str);

function ShowCode_Copy(ob, ev, ctxt)
str = get(ctxt,'String');
sel = get(ctxt,'Value');
str = str(sel);
clipboard('copy',sprintf('%s\n',str{:}));

function ShowCode_SelAll(ob, ev, ctxt)
set(ctxt,'Value', 1:numel(get(ctxt, 'String')));

function ShowCode_UnSelAll(ob, ev, ctxt)
set(ctxt,'Value', []);