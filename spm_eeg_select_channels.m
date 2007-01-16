function varargout = spm_eeg_select_channels(varargin)
% SPM_EEG_SELECT_CHANNELS M-file for spm_eeg_select_channels.fig
%      SPM_EEG_SELECT_CHANNELS, by itself, creates a new SPM_EEG_SELECT_CHANNELS or raises the existing
%      singleton*.
%
%      H = SPM_EEG_SELECT_CHANNELS returns the handle to a new SPM_EEG_SELECT_CHANNELS or the handle to
%      the existing singleton*.
%
%      SPM_EEG_SELECT_CHANNELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPM_EEG_SELECT_CHANNELS.M with the given input arguments.
%
%      SPM_EEG_SELECT_CHANNELS('Property','Value',...) creates a new SPM_EEG_SELECT_CHANNELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_eeg_select_channels_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spm_eeg_select_channels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spm_eeg_select_channels

% Last Modified by GUIDE v2.5 13-Jul-2005 14:18:07

%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_select_channels.m 716 2007-01-16 21:13:50Z karl $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_eeg_select_channels_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_eeg_select_channels_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before spm_eeg_select_channels is made visible.
function spm_eeg_select_channels_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_eeg_select_channels (see VARARGIN)

% Choose default command line output for spm_eeg_select_channels
handles.output = hObject;


% display graph of channels
D = varargin{1};
P = fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf);
load(P)

% correct Cnames (channel template file entries can have more than one
% entry per channel)
for i = 1:Nchannels
    if iscell(Cnames{i})
        Cnames{i} = Cnames{i}{1};
    end
end

figure(handles.figure1);

Xrec = [-0.01 0.01 0.01 -0.01];
Yrec = [-0.01 -0.01 0.01 0.01];

Hpatch = cell(1, Nchannels);
Htext  = cell(1, Nchannels);

Cselect = [1 1 1];
Cdeselect = [0.5 0.5 0.5];

ind = zeros(1,length(D.channels.order));
ind(D.gfx.channels) = 1;
tmp = D.channels.order; % if D.channels.order isn't replaced by this variable -> strange error in 7.04

for i = 1:length(D.channels.order)
    if ~ind(i)
        Hpatch{i} = patch(Xrec+Cpos(1,tmp(i)), Yrec+Cpos(2,tmp(i)), Cdeselect, 'EdgeColor', 'none');
    else
        Hpatch{i} = patch(Xrec+Cpos(1,tmp(i)), Yrec+Cpos(2,tmp(i)), Cselect, 'EdgeColor', 'none');        
    end
    Hpatch{i};
    xy = get(Hpatch{i}, 'Vertices');
    Htext{i} = text(min(xy(:,1)), max(xy(:,2)), Cnames{D.channels.order(i)});
    set(Htext{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
end
axis square
axis off

% listbox
set(handles.listbox1, 'String', Cnames(D.channels.order([1:length(ind)])));
set(handles.listbox1, 'Value', find(ind));

handles.D = D;
handles.ind = ind;
handles.Cnames = Cnames;
handles.Cpos = Cpos;
handles.Nchannels = Nchannels;
handles.Hpatch = Hpatch;
handles.Htext = Htext;
handles.Cselect = Cselect;
handles.Cdeselect = Cdeselect;

for i = 1:length(D.gfx.channels)
    % callbacks for patches and text
    if D.gfx.channels(i) ~=0
        set(Hpatch{i}, 'ButtonDownFcn', {@patch_select, handles, i});
        set(Htext{i}, 'ButtonDownFcn', {@patch_select, handles, i});
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spm_eeg_select_channels wait for user response (see UIRESUME)
uiwait(handles.figure1);

function patch_select(obj, eventdata, handles, i, ind);

s = get(handles.listbox1, 'Value');
h = handles.Hpatch{i};

if ~ismember(i, s)
    set(h, 'FaceColor', handles.Cselect);
    set(handles.listbox1, 'Value', sort([s, i]));
else
    set(h, 'FaceColor', handles.Cdeselect);
    set(handles.listbox1, 'Value', setdiff(s, i));
end

% --- Outputs from this function are returned to the command line.
function varargout = spm_eeg_select_channels_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.output;
close(handles.figure1);


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

s = get(handles.listbox1, 'Value');

for i = 1:length(handles.ind)
    if ismember(i, s)
        set(handles.Hpatch{i}, 'FaceColor', handles.Cselect);
    else
        set(handles.Hpatch{i}, 'FaceColor', handles.Cdeselect);
    end
end

% --- Executes on button press in select.
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Cycle through patch elements and recolour
for i = 1:length(handles.ind)
    set(handles.Hpatch{i}, 'FaceColor', handles.Cselect);
end

% set all listbox entries to selected
set(handles.listbox1, 'Value', [1:length(handles.ind)]);

% --- Executes on button press in deselect.
function deselect_Callback(hObject, eventdata, handles)
% hObject    handle to deselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Cycle through patch elements and recolour
for i = 1:length(handles.ind)
    set(handles.Hpatch{i}, 'FaceColor', handles.Cdeselect);
end

% remove all listbox entries
set(handles.listbox1, 'Value', []);

% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[P1, P2] = uigetfile('*.mat', 'Select file to load from');
load(fullfile(P2, P1), 'Iselectedchannels');

if ~exist('Iselectedchannels', 'var')
    errordlg('This file doesn''t contain channel indices', 'Wrong file?');
else
    if max(Iselectedchannels) > length(handles.ind)
        errordlg('This file doesn''t channel indices for these data', 'Wrong file?');
    else
        set(handles.listbox1, 'Value', Iselectedchannels);
        
        for i = 1:handles.Nchannels
            if isempty(find(Iselectedchannels == i))
                set(handles.Hpatch{i}, 'FaceColor', handles.Cdeselect);
            else
                set(handles.Hpatch{i}, 'FaceColor', handles.Cselect);
            end
        end
    end
end

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Iselectedchannels = get(handles.listbox1, 'Value');

[P1, P2] = uiputfile('*.mat', 'Choose file name to save to');

if spm_matlab_version_chk('7') >= 0
    save(fullfile(P2, P1), '-V6', 'Iselectedchannels');
else
    save(fullfile(P2, P1), 'Iselectedchannels');
end

% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% returns vector of channel indices that are to be displayed
if isempty(get(handles.listbox1, 'Value'))
    h = errordlg('Must select at least one channel!', 'Selection error', 'modal');
else
    % return indices of selected channels
    handles.output = get(handles.listbox1, 'Value');
    guidata(hObject, handles);
    uiresume;
end

% --- Executes on button press in removebad.
function removebad_Callback(hObject, eventdata, handles)
% hObject    handle to removebad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

D = handles.D;

s = get(handles.listbox1, 'Value');

s = setdiff(s, D.channels.Bad);

set(handles.listbox1, 'Value', s);

for i = 1:handles.Nchannels
    if isempty(find(s == i))
        set(handles.Hpatch{i}, 'FaceColor', handles.Cdeselect);
    else
        set(handles.Hpatch{i}, 'FaceColor', handles.Cselect);
    end
end

