function varargout = Neural_demo(varargin)
% NEURAL_DEMO M-file for Neural_demo.fig
%      NEURAL_DEMO, by itself, creates a new NEURAL_DEMO or raises the existing
%      singleton*.
%
%      H = NEURAL_DEMO returns the handle to a new NEURAL_DEMO or the handle to
%      the existing singleton*.
%
%      NEURAL_DEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURAL_DEMO.M with the given input arguments.
%
%      NEURAL_DEMO('Property','Value',...) creates a new NEURAL_DEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Neural_demo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Neural_demo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Neural_demo

% Last Modified by GUIDE v2.5 09-Apr-2008 14:15:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Neural_demo_OpeningFcn, ...
                   'gui_OutputFcn',  @Neural_demo_OutputFcn, ...
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


% --- Executes just before Neural_demo is made visible.
function Neural_demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Neural_demo (see VARARGIN)

% Choose default command line output for Neural_demo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Neural_demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Neural_demo_OutputFcn(hObject, eventdata, handles) 
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
drawnow

% get graphics handle
spm_figure('GetWin','MFM');
clf
eval(file)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_lfp_demo')

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_ind_demo')

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_mtf_demo')

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_csd_demo')

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_sigmoid_demo')

% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_mfm_demo')

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_nested_oscillations_demo')



