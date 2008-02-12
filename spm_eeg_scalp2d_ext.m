function varargout = spm_eeg_scalp2d_ext(varargin)
% SPM_EEG_SCALP2D_EXT M-file for spm_eeg_scalp2d_ext.fig
%      SPM_EEG_SCALP2D_EXT, by itself, creates a new SPM_EEG_SCALP2D_EXT or raises the existing
%      singleton*.
%
%      H = SPM_EEG_SCALP2D_EXT returns the handle to a new SPM_EEG_SCALP2D_EXT or the handle to
%      the existing singleton*.
%
%      SPM_EEG_SCALP2D_EXT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPM_EEG_SCALP2D_EXT.M with the given input arguments.
%
%      SPM_EEG_SCALP2D_EXT('Property','Value',...) creates a new SPM_EEG_SCALP2D_EXT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_eeg_scalp2d_ext_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spm_eeg_scalp2d_ext_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help spm_eeg_scalp2d_ext

% Last Modified by GUIDE v2.5 22-Nov-2005 17:06:20

% Colon removed so that times output to Matlab window	Doris Eckstein
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_eeg_scalp2d_ext_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_eeg_scalp2d_ext_OutputFcn, ...
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


% --- Executes just before spm_eeg_scalp2d_ext is made visible.
function spm_eeg_scalp2d_ext_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_eeg_scalp2d_ext (see VARARGIN)

% Choose default command line output for spm_eeg_scalp2d_ext
handles.output = hObject;

handles.D = varargin{1};
handles.T = varargin{2}; % input in peri-stimulus time (ms)

handles.ms = [-handles.D.events.start:handles.D.events.stop]*1000/handles.D.Radc;

for i = 1:length(handles.T)
	tmp = (handles.T(i) - handles.ms).^2;
	[m, ind(i)] = min(tmp);
end

handles.T = ind;

% locations
CTF = load(fullfile(spm('dir'), 'EEGtemplates', handles.D.channels.ctf));
%CTF.Cpos = CTF.Cpos(:, handles.D.channels.order(handles.D.channels.eeg));
handles.D.gfx.channels = intersect(handles.D.gfx.channels,handles.D.channels.eeg);	% to ensure EOG not included
CTF.Cpos = CTF.Cpos(:, handles.D.channels.order(handles.D.gfx.channels));

handles.x = min(CTF.Cpos(1,:)):0.005:max(CTF.Cpos(1,:));
handles.y = min(CTF.Cpos(2,:)):0.005:max(CTF.Cpos(2,:));

[handles.x1, handles.y1] = meshgrid(handles.x, handles.y);
handles.xp = CTF.Cpos(1,:)';
handles.yp = CTF.Cpos(2,:)';

handles.event = varargin{3};

if length(handles.T) > 1
    % average, remove slider
    set(handles.slider1, 'Visible', 'off');
    set(handles.text2, 'Visible', 'off');
else
    % set slider's range and initial value
    set(handles.slider1, 'min', handles.ms(1));
    set(handles.slider1, 'max', handles.ms(end));
    set(handles.slider1, 'Value', handles.ms(handles.T));
end

% Update handles structure
guidata(hObject, handles);

plot_spatial(hObject, handles);

% UIWAIT makes spm_eeg_scalp2d_ext wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spm_eeg_scalp2d_ext_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

T = get(handles.slider1, 'Value')

tmp = (T - handles.ms).^2;
[m, i] = min(tmp);

handles.T = i;

guidata(hObject, handles);

plot_spatial(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function plot_spatial(hObject, handles)

T = handles.T;
D = handles.D;
event = handles.event;

% data
if length(T) == 1
    d = squeeze(D.data(D.gfx.channels, T, event));
    
    if ~isfield(handles, 'Colourbar')
        handles.CLim1 = min(min(D.data(setdiff(D.gfx.channels, D.channels.Bad), :, event)));
        handles.CLim2 = max(max(D.data(setdiff(D.gfx.channels, D.channels.Bad), :, event)));
    end
else
	d = squeeze(mean(D.data(D.gfx.channels, T, event), 2));
end

%Exclude bad channels
badchan = intersect(D.gfx.channels,D.channels.Bad);
if ~isempty(badchan)
        d(badchan) = NaN;
end

z = griddata(handles.xp, handles.yp, d, handles.x1, handles.y1);

if length(T) == 1
    set(handles.text1, 'String', sprintf('%d ms', round(handles.ms(T))));
else
    set(handles.text1, 'String', 'average');
end

axes(handles.axes1);
cla
surface(handles.x, handles.y, z);
axis off
shading('interp')
hold on
plot3(handles.xp, handles.yp, d, 'k.');

if ~isfield(handles, 'Colourbar')
    set(handles.axes1, 'CLim', [handles.CLim1 handles.CLim2])
    handles.Colourbar = colorbar;
    ylabel(handles.Colourbar, D.units, 'FontSize', 16);
end

guidata(hObject, handles);

drawnow




