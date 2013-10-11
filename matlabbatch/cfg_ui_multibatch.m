function varargout = cfg_ui_multibatch(varargin)
% CFG_UI_MULTIBATCH MATLAB code for cfg_ui_multibatch.fig
%      CFG_UI_MULTIBATCH, by itself, creates a new CFG_UI_MULTIBATCH or raises the existing
%      singleton*.
%
%      H = CFG_UI_MULTIBATCH returns the handle to a new CFG_UI_MULTIBATCH or the handle to
%      the existing singleton*.
%
%      CFG_UI_MULTIBATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CFG_UI_MULTIBATCH.M with the given input arguments.
%
%      CFG_UI_MULTIBATCH('Property','Value',...) creates a new CFG_UI_MULTIBATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cfg_ui_multibatch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cfg_ui_multibatch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cfg_ui_multibatch

% Last Modified by GUIDE v2.5 24-Sep-2013 15:11:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cfg_ui_multibatch_OpeningFcn, ...
                   'gui_OutputFcn',  @cfg_ui_multibatch_OutputFcn, ...
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


% --- Executes when user attempts to close cfg_ui_multibatch.
function cfg_ui_multibatch_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to cfg_ui_multibatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
try
    cfg_util('deljob', handles.ud.cjob);
end
delete(hObject);

% --- Executes just before cfg_ui_multibatch is made visible.
function cfg_ui_multibatch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cfg_ui_multibatch (see VARARGIN)

% Choose default command line output for cfg_ui_multibatch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cfg_ui_multibatch wait for user response (see UIRESUME)
% uiwait(handles.cfg_ui_multibatch);

% Register fatal errors with cfg_message
cfg_message('error', 'level', 'cfg_ui_multibatch:wrongclass');

% If called with extra arguments, try to init them as job
if ~isempty(varargin)
    local_initjob(handles, varargin)
end

% --- Outputs from this function are returned to the command line.
function varargout = cfg_ui_multibatch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes when cfg_ui_multibatch is resized.
function cfg_ui_multibatch_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to cfg_ui_multibatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in BtnRowAddRow.
function BtnRowAddRow_Callback(hObject, eventdata, handles)
% hObject    handle to BtnRowAddRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_AddRow(handles);

% --- Executes on button press in BtnRowDelRow.
function BtnRowDelRow_Callback(hObject, eventdata, handles)
% hObject    handle to BtnRowDelRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_DelRow(handles);

% --- Executes on button press in BtnValEditVal.
function BtnValEditVal_Callback(hObject, eventdata, handles)
% hObject    handle to BtnValEditVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MenuFileLoad_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
jobfiles = cfg_getfile([1 inf], 'batch', 'Select batch file(s)');
if ~isempty(jobfiles)
    local_initjob(handles, jobfiles);
end

% --------------------------------------------------------------------
function MenuFileSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[job_ids, job_sts] = local_filljobs(handles);
[file, pth, idx] = uiputfile({'*.mat','Matlab .mat File';...
                    '*.m','Matlab .m Script File'}, 'Basename for jobs');
if isnumeric(file) && file == 0
    return;
end;
[p, n, e] = fileparts(file);
if isempty(e) || ~any(strcmp(e,{'.mat','.m'}))
    e1 = {'.mat','.m'};
    e2 = e1{idx};
    file = sprintf('%s%s', n, e);
else
    file = n;
    e2 = e;
end
file = sprintf('%s_%%0%dd', file, floor(log10(numel(job_ids))+1));
try
    for cj = 1:numel(job_ids)
        cfg_util('savejob', job_ids{cj}, fullfile(pth, [sprintf(file, cj) e2]));
    end
catch
    l = lasterror;
    errordlg(l.message,'Error saving job', 'modal');
end

% --------------------------------------------------------------------
function MenuFileRun_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[job_ids, job_sts] = local_filljobs(handles);
job_ids = job_ids(job_sts);
for cj = 1:numel(job_ids)
    cfg_util('run', job_ids{cj});
end

% --------------------------------------------------------------------
function MenuFileQuit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.cfg_ui_multibatch);

% --------------------------------------------------------------------
function MenuEditRowAddRow_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditRowAddRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_AddRow(handles);

% --------------------------------------------------------------------
function MenuEditRowDelRow_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditRowDelRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
local_DelRow(handles);

% --------------------------------------------------------------------
function MenuEditValEditVal_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditValEditVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when selected cell(s) is changed in default_table.
function default_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to default_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if size(eventdata.Indices, 1) == 1
    local_showvaledit_default(handles, eventdata.Indices)
end

% --- Executes on selection change in valshow.
function valshow_Callback(hObject, eventdata, handles)
% hObject    handle to valshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns valshow contents as cell array
%        contents{get(hObject,'Value')} returns selected item from valshow

% --- Executes during object creation, after setting all properties.
function valshow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valshow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in helpbox.
function helpbox_Callback(hObject, eventdata, handles)
% hObject    handle to helpbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns helpbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from helpbox

% --- Executes during object creation, after setting all properties.
function helpbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to helpbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected cell(s) is changed in values_table.
function values_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to values_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if size(eventdata.Indices, 1) == 1
    if isequal(get(handles.cfg_ui_multibatch, 'CurrentCharacter'), char(13))
        % handle RETURN character - edit current cell instead of new one.
        % After that, advance to new cell (this behaviour is unwanted, but
        % MATLAB does not allow to set the selected cell otherwise).
        if ~isempty(handles.ud.values_editcb)
            feval(handles.ud.values_editcb);
            handles = guidata(handles.cfg_ui_multibatch);
        end
        set(handles.cfg_ui_multibatch, 'CurrentCharacter', char(0));
    end
    local_showvaledit_values(handles, eventdata.Indices);
    handles = guidata(hObject);
end
handles.ud.values_indices = eventdata.Indices;
guidata(handles.cfg_ui_multibatch, handles);

% --------------------------------------------------------------------
function local_AddRow(handles)
handles.ud.values_contents = [handles.ud.values_contents; handles.ud.default_contents];
guidata(handles.cfg_ui_multibatch, handles);
local_show_values(handles);

% --------------------------------------------------------------------
function local_DelRow(handles)
if ~isempty(handles.ud.values_indices)
    todel = handles.ud.values_indices(:,1);
    handles.ud.values_contents(todel, :) = [];
    handles.ud.values_editcb = [];
    handles.ud.values_indices = [];
    guidata(handles.cfg_ui_multibatch, handles);
    local_show_values(handles);
end

% --------------------------------------------------------------------
function local_initjob(handles, jobs)
try
    if numel(jobs) == 1 && cfg_util('isjob_id', jobs{1})
        cjob = jobs{1};
    else
        cjob = cfg_util('initjob', jobs);
    end
catch le
    cfg_message(le);
end
[mod_job_idlist, mod_names, mod_item_idx, item_mod_idlists, item_names] = cfg_util('getopeninputs', cjob);
if isempty(mod_job_idlist)
    cfg_util('deljob', cjob);
    cfg_message('cfg_ui_multibatch:noopenitems', 'nothing to fill in');
else
    item_names   = [item_names{:}];
    handles.ud.cjob = cjob;
    handles.ud.mod_job_idlist = mod_job_idlist(mod_item_idx);
    handles.ud.item_mod_idlist = [item_mod_idlists{:}];
    handles.ud.default_contents  = cell(size(item_names));
    handles.ud.values_editcb = [];
    handles.ud.values_indices = [];
    for citem=1:numel(item_names)
        handles.ud.default_contents{citem} = cfg_ui_util('showitem', ...
            {handles.ud.cjob, handles.ud.mod_job_idlist{citem}, ...
            handles.ud.item_mod_idlist{citem}}, false);
        if ~any(strcmp(handles.ud.default_contents{citem}{5}, {'cfg_entry', 'cfg_files', 'cfg_menu'}))
            cfg_util('deljob', cjob);
            cfg_message('cfg_ui_multibatch:wrongclass', 'can not fill in items of class "%s"', handles.ud.default_contents{citem}{5});
        end
    end
    column_names = cellfun(@(mod,item)sprintf('%s|%s', mod, item), mod_names(mod_item_idx), item_names, 'UniformOutput', false);
    set(handles.default_table, 'ColumnName', column_names, 'Data', cell(size(item_names)));
    handles.ud.values_contents = cell(0, size(item_names, 2));
    set(handles.values_table, 'ColumnName', column_names, 'Data', cell(0, size(item_names, 2)));
    guidata(handles.cfg_ui_multibatch, handles);
end

% --------------------------------------------------------------------
function local_show_values(handles)
handles = guidata(handles.cfg_ui_multibatch);
datastr = cell(size(handles.ud.values_contents));
for k = 1:numel(datastr)
    [~, datastr{k}] = cfg_ui_util('showitemstr', handles.ud.values_contents{k}, false);
end
set(handles.values_table, 'Data', datastr);

% --------------------------------------------------------------------
function local_showvaledit_default(handles, index)
ciid = {handles.ud.cjob, ...
    handles.ud.mod_job_idlist{index(2)}, ...
    handles.ud.item_mod_idlist{index(2)}};
contents = handles.ud.default_contents{index(2)};
cfg_ui_util('showvaledit', handles.cfg_ui_multibatch, ciid, contents, [], false, ...
    @(nval)local_set_defaultstable(handles, index, nval), ...
    @()local_update_defaultstable(handles, index));
column_names = get(handles.default_table, 'ColumnName');
set(findobj(handles.cfg_ui_multibatch, '-regexp', 'Tag','.*EditVal$'), ...
    'Callback', @(ob,ev)cfg_ui_util('valedit_editvalue', ciid, ...
    column_names{index(2)}, contents{2}));

% --------------------------------------------------------------------
function local_showvaledit_values(handles, index)
ciid = {handles.ud.cjob, ...
    handles.ud.mod_job_idlist{index(2)}, ...
    handles.ud.item_mod_idlist{index(2)}};
contents = handles.ud.values_contents{index(1), index(2)};
cfg_ui_util('showvaledit', handles.cfg_ui_multibatch, ciid, contents, [], false, ...
    @(nval)local_set_valuestable(handles, index, nval), ...
    @()local_show_values(handles));
column_names = get(handles.default_table, 'ColumnName');
handles.ud.values_editcb = @(ob,ev)cfg_ui_util('valedit_editvalue', ciid, ...
    column_names{index(2)}, contents{2});
set(findobj(handles.cfg_ui_multibatch, '-regexp', 'Tag','.*EditVal$'), ...
    'Callback', handles.ud.values_editcb);
guidata(handles.cfg_ui_multibatch, handles);

% --------------------------------------------------------------------
function local_set_defaultstable(handles, index, val)
switch handles.ud.default_contents{index(2)}{5}
    case {'cfg_entry','cfg_files'}
        handles.ud.default_contents{index(2)}{2} = {val};
    case {'cfg_menu'}
        handles.ud.default_contents{index(2)}{2} = handles.ud.default_contents{index(2)}{4}(val);
end        
handles.ud.default_contents{index(2)}{7} = true;
handles.ud.default_contents{index(2)}{8} = true;
for k = 1:size(handles.ud.values_contents,1)
    handles.ud.values_contents{k,index(2)}{2} = handles.ud.default_contents{index(2)}{2};
    handles.ud.values_contents{k,index(2)}{7} = true;
    handles.ud.values_contents{k,index(2)}{8} = true;
end
guidata(handles.cfg_ui_multibatch, handles);

% --------------------------------------------------------------------
function local_set_valuestable(handles, index, val)
switch handles.ud.default_contents{index(2)}{5}
    case {'cfg_entry','cfg_files'}
        handles.ud.values_contents{index(1), index(2)}{2} = {val};
    case {'cfg_menu'}
        handles.ud.values_contents{index(1), index(2)}{2} = handles.ud.default_contents{index(2)}{4}(val);
end        
handles.ud.values_contents{index(1), index(2)}{7} = true;
handles.ud.values_contents{index(1), index(2)}{8} = true;
guidata(handles.cfg_ui_multibatch, handles);

% --------------------------------------------------------------------
function local_update_defaultstable(handles, index)
handles = guidata(handles.cfg_ui_multibatch);
[~, datastr] = cfg_ui_util('showitemstr', handles.ud.default_contents{index(2)}, false);
dat = get(handles.default_table, 'Data');
dat{index(1), index(2)} = datastr;
set(handles.default_table, 'Data', dat);
local_showvaledit_default(handles, index);
local_show_values(handles);

% --------------------------------------------------------------------
function [job_ids, job_sts] = local_filljobs(handles)
% retrieve data from values_table and fill jobs
job_ids = cell(size(handles.ud.values_contents, 1), 1);
job_sts = false(size(handles.ud.values_contents, 1), 1);
for cj = 1:size(handles.ud.values_contents, 1)
    values = cell(1, size(handles.ud.values_contents, 2));
    for cv = 1:size(handles.ud.values_contents, 2)
        if ~isempty(handles.ud.values_contents{cj, cv}{2})
            values{cv} = handles.ud.values_contents{cj, cv}{2}{1};
        end
    end
    job_ids{cj} = cfg_util('clonejob', handles.ud.cjob);
    job_sts(cj) = cfg_util('filljob', job_ids{cj}, values{:});
end
