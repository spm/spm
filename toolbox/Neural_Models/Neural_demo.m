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

% Last Modified by GUIDE v2.5 05-Aug-2016 16:41:47

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

% default paper
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Mechanisms%20of%20evoked%20and%20induced%20responses.pdf';

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
h      = help(file);
str{1} = [file ':'];
str{2} = '__________________________________________________________________________ ';
str{3} = ' ';
str{4} = h;
set(handles.help,'String',str);
handles.file = file;
guidata(hObject, handles);

% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
set(handles.pushbutton24,'String','please wait')
drawnow
eval(handles.file)
set(handles.pushbutton24,'String','run demo')


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
try, edit(handles.file); end

% -------------------------------------------------------------------------
function uipanel7_ButtonDownFcn(hObject, eventdata, handles)
try, web(handles.web);   end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/A%20neural%20mass%20model%20for%20MEG.pdf';
run_demo_Callback(hObject, handles, 'spm_lfp_demo')

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Mechanisms%20of%20evoked%20and%20induced%20responses.pdf';
run_demo_Callback(hObject, handles, 'spm_ind_demo')

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20causal%20modeling%20of%20evoked%20responses%20in%20EEG%20and%20MEG.pdf';
run_demo_Callback(hObject, handles, 'spm_mtf_demo')

% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20causal%20models%20of%20steady-state%20responses.pdf';
run_demo_Callback(hObject, handles, 'spm_csd_demo')

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Population%20dynamics%20Variance%20and%20the%20sigmoid%20activation%20function.pdf';
run_demo_Callback(hObject, handles, 'spm_sigmoid_demo')

% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Stochastic%20models%20of%20neuronal%20dynamics.pdf';
run_demo_Callback(hObject, handles, 'spm_mfm_demo')

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Population%20dynamics%20under%20the%20Laplace%20assumption.pdf';
run_demo_Callback(hObject, handles, 'spm_nested_oscillations_demo')

% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_seizure_demo')

% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Mechanisms%20of%20evoked%20and%20induced%20responses.pdf';
run_demo_Callback(hObject, handles, 'spm_induced_demo')

% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, ~, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Mechanisms%20of%20evoked%20and%20induced%20responses.pdf';
run_demo_Callback(hObject, handles, 'spm_induced_optimise')


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'spm_induced_optimise_parameters')


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Stochastic%20models%20of%20neuronal%20dynamics.pdf';
run_demo_Callback(hObject, handles, 'spm_mfa_demo')


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Population%20dynamics%20under%20the%20Laplace%20assumption.pdf';
run_demo_Callback(hObject, handles, 'spm_dcm_prior_responses')


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Post%20hoc%20Bayesian%20model%20selection.pdf';
run_demo_Callback(hObject, handles, 'DEMO_model_reduction_ERP')

% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Analysing%20connectivity%20with%20Granger%20causality%20and%20dynamic%20causal%20modelling.pdf';
run_demo_Callback(hObject, handles, 'spm_dcm_Granger_demo')


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Analysing%20connectivity%20with%20Granger%20causality%20and%20dynamic%20causal%20modelling.pdf';
run_demo_Callback(hObject, handles, 'spm_delays_demo')


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Analysing%20connectivity%20with%20Granger%20causality%20and%20dynamic%20causal%20modelling.pdf';
run_demo_Callback(hObject, handles, 'spm_dcm_Granger_asymmetry_demo')


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_epileptor_demo')


% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'spm_erp2csd_demo')


% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
run_demo_Callback(hObject, handles, 'DEMO_dcm_fmri_nnm')

