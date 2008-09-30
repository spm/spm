function varargout = spm_api_nmm(varargin)
% SPM_API_NMM M-file for spm_api_nmm.fig
%      SPM_API_NMM, by itself, creates a new SPM_API_NMM or raises the existing
%      singleton*.
%
%      H = SPM_API_NMM returns the handle to a new SPM_API_NMM or the handle to
%      the existing singleton*.
%
%      SPM_API_NMM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPM_API_NMM.M with the given input arguments.
%
%      SPM_API_NMM('Property','Value',...) creates a new SPM_API_NMM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_api_nmm_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spm_api_nmm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spm_api_nmm

% Last Modified by GUIDE v2.5 26-Sep-2008 20:21:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_api_nmm_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_api_nmm_OutputFcn, ...
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


% --- Executes just before spm_api_nmm is made visible.
function spm_api_nmm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_api_nmm (see VARARGIN)

% Choose default command line output for spm_api_nmm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spm_api_nmm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spm_api_nmm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%==========================================================================


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)

[f,p]       = uigetfile('*.mat','please select DCM file'); cd(p)
name        = fullfile(p,f);
DCM         = load(name,'-mat');
DCM         = DCM.DCM;
DCM.name    = name;
handles.DCM = DCM;
guidata(hObject,handles);




% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
try
   [p,file]  = fileparts(handles.DCM.name);
catch
    try
        [p,file] = fileparts(handles.DCM.xY.Dfile);
        file     = ['DCM_' file];
    catch
        file     = ['DCM' date];
    end
end
[file,fpath]     = uiputfile(['DCM*.mat'],'DCM file to save',file);
 
if fpath
    handles.DCM.name = fullfile(fpath,file);
    set(handles.name,'String',file);
    DCM              = handles.DCM;
    save(DCM.name,'DCM')
    cd(fpath)
end
 
% assign in base
%--------------------------------------------------------------------------
assignin('base','DCM',handles.DCM)
guidata(hObject,handles);


% --- Executes on button press in kernels.
function kernels_Callback(hObject, eventdata, handles)
clear spm_erp_L spm_gen_erp
DCM = handles.DCM;
M   = DCM.M;

% inital states
%--------------------------------------------------------------------------
[x N] = spm_x_nmm(M.pE);
M.x   = x;
M.GE  = N.GE;
M.GI  = N.GI;
M.Cx  = N.Cx;
M.f   = 'spm_fx_mfm';
M.G   = 'spm_lx_erp';

% prediction (source space)
%--------------------------------------------------------------------------
x   = feval(M.IS,M.pE,M,DCM.xU);              







