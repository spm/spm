function varargout = DEM_demo(varargin)
% DEM_DEMO M-file for DEM_demo.fig
%      DEM_DEMO, by itself, creates a new DEM_DEMO or raises the existing
%      singleton*.
%
%      H = DEM_DEMO returns the handle to a new DEM_DEMO or the handle to
%      the existing singleton*.
%
%      DEM_DEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEM_DEMO.M with the given input arguments.
%
%      DEM_DEMO('Property','Value',...) creates a new DEM_DEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DEM_demo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DEM_demo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DEM_demo

% Last Modified by GUIDE v2.5 23-Aug-2007 15:32:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DEM_demo_OpeningFcn, ...
                   'gui_OutputFcn',  @DEM_demo_OutputFcn, ...
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


% --- Executes just before DEM_demo is made visible.
function DEM_demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DEM_demo (see VARARGIN)

% Choose default command line output for DEM_demo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DEM_demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DEM_demo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function run_demo_Callback(hObject, handles, file)
h    = help(file);
str{1} = [file ':'];
str{2} = '______________________________________________________________ ';
str{2} = ' ';
str{3} = h;
set(handles.help,'String',str);
eval(file)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_GLM')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_PEB')

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_factor_analysis')

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_OU')

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_convolution')

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_EM')

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_DEM')

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_filtering')

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_Lorenz')

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_DFP')

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_hdm')

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEM_demo_double_well')

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DFP_demo_double_well')


