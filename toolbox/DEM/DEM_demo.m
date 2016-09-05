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

% Last Modified by GUIDE v2.5 03-Sep-2016 16:56:31

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

% default paper
handles.web    = 'http://www.fil.ion.ucl.ac.uk/~karl/The%20free-energy%20principle%20A%20unified%20brain%20theory.pdf';

% Display PDF image
axes5_CreateFcn(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
imagesc(imread('PDF.jpg')), axis off


% --- Outputs from this function are returned to the command line.
function varargout = DEM_demo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function run_demo_Callback(hObject, handles, file)
if isdeployed
    h = sprintf('%s: MATLAB code and help are not available in SPM Standalone.',file);
else
    h  = help(file);
end
str{1} = [file ':'];
str{2} = '__________________________________________________________________________ ';
str{3} = ' ';
str{4} = h;
set(handles.help,'String',str);
handles.file = file;
guidata(hObject, handles);


% --- Executes on button press in pushbutton131.
function pushbutton131_Callback(hObject, eventdata, handles)
try, web(handles.web,'-browser'); end


% --- Executes on button press in pushbutton51.
function pushbutton51_Callback(hObject, eventdata, handles)

set(handles.pushbutton51,'String','please wait')
drawnow
try
    guidata(1,handles);
catch
    spm_figure('GetWin','DEM');
    guidata(1,handles);
end
eval(handles.file)
handles = set(0,'UserData');
handles = guidata(1);
set(handles.pushbutton51,'String','run demo')

% --- Executes on button press in pushbutton93.
function pushbutton93_Callback(hObject, eventdata, handles)
try
    if isdeployed
        set(hObject,'String','MATLAB code (not available)');
        set(hObject,'Enable','Off');
    else
        edit(handles.file);
    end
end



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Variational%20free%20energy%20and%20the%20Laplace%20approximation.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_GLM')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Variational%20free%20energy%20and%20the%20Laplace%20approximation.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_PEB')

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Hierarchical%20Models%20in%20the%20Brain.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_factor_analysis')

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/DEM%20A%20variational%20treatment%20of%20dynamic%20systems.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_OU')

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Hierarchical%20Models%20in%20the%20Brain.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_convolution')

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/DEM%20A%20variational%20treatment%20of%20dynamic%20systems.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_EM')

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/DEM%20A%20variational%20treatment%20of%20dynamic%20systems.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_DEM')

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/DEM%20A%20variational%20treatment%20of%20dynamic%20systems.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_filtering')

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/DEM%20A%20variational%20treatment%20of%20dynamic%20systems.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_Lorenz')

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20causal%20modelling.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm')

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/DEM%20A%20variational%20treatment%20of%20dynamic%20systems.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_double_well')

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Variational%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_DFP')

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Variational%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DFP_demo_double_well')

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Cortical%20circuits%20for%20perceptual%20inference.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_song_priors')

% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Predictive%20coding%20under%20the%20free-energy%20principle.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_song_inference')

% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Cortical%20circuits%20for%20perceptual%20inference.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_song_omission')

% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Learning%20and%20inference%20in%20the%20brain.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_face_inference')

% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Attractors%20in%20song.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MMN')

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_Gabor')

% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Action%20and%20behavior%20A%20free-energy%20formulation.pdf';
run_demo_Callback(hObject, handles, 'ADEM_visual')

% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Action%20and%20behavior%20A%20free-energy%20formulation.pdf';
run_demo_Callback(hObject, handles, 'ADEM_motor')

% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Reinforcement%20Learning%20or%20Active%20Inference.pdf';
run_demo_Callback(hObject, handles, 'ADEM_learning')

% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Action%20and%20behavior%20A%20free-energy%20formulation.pdf';
run_demo_Callback(hObject, handles, 'ADEM_lorenz')

% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Action%20and%20behavior%20A%20free-energy%20formulation.pdf';
run_demo_Callback(hObject, handles, 'ADEM_reaching')

% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Free%20Energy%20Value%20and%20Attractors.pdf';
run_demo_Callback(hObject, handles, 'ADEM_lorenz_entropy')

% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Free%20Energy%20Value%20and%20Attractors.pdf';
run_demo_Callback(hObject, handles, 'ADEM_mountaincar_loss')

% --- Executes on button press in pushbutton80.
function pushbutton80_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Policies%20and%20Priors.pdf';
run_demo_Callback(hObject, handles, 'ADEM_SHC_demo')

% --- Executes on button press in pushbutton81.
function pushbutton81_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Free%20Energy%20Value%20and%20Attractors.pdf';
run_demo_Callback(hObject, handles, 'ADEM_lorenz_surprise')

% --- Executes on button press in pushbutton83.
function pushbutton83_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Generalised%20Filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_LAP')

% --- Executes on button press in pushbutton84.
function pushbutton84_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Generalised%20Filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm_LAP')

% --- Executes on button press in pushbutton91.
function pushbutton91_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Attention%20uncertainty%20and%20free-energy.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_Posner')

% --- Executes on button press in pushbutton92.
function pushbutton92_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Policies%20and%20Priors.pdf';
run_demo_Callback(hObject, handles, 'ADEM_cost_SHC')

% --- Executes on button press in pushbutton94.
function pushbutton94_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, '%%%% DEM_demo_CM_Lorenz %%%')

% --- Executes on button press in pushbutton96.
function pushbutton96_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Action%20understanding%20and%20active%20inference.pdf';
run_demo_Callback(hObject, handles, 'ADEM_observe')

% --- Executes on button press in pushbutton97.
function pushbutton97_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Network%20discovery%20with%20DCM.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_DCM_LAP')

% --- Executes on button press in pushbutton98.
function pushbutton98_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20modeling%20of%20neuronal%20responses%20in%20fMRI%20using%20cubature%20Kalman%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_convolution_LAP')

% --- Executes on button press in pushbutton99.
function pushbutton99_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20modeling%20of%20neuronal%20responses%20in%20fMRI%20using%20cubature%20Kalman%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_lorenz_LAP')

% --- Executes on button press in pushbutton100.
function pushbutton100_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20modeling%20of%20neuronal%20responses%20in%20fMRI%20using%20cubature%20Kalman%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_doublewell_LAP')

% --- Executes on button press in pushbutton101.
function pushbutton101_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20modeling%20of%20neuronal%20responses%20in%20fMRI%20using%20cubature%20Kalman%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm_SCK')

% --- Executes on button press in pushbutton102.
function pushbutton102_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Action%20understanding%20and%20active%20inference.pdf';
run_demo_Callback(hObject, handles, 'ADEM_writing')

% --- Executes on button press in pushbutton103.
function pushbutton103_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Free%20energy%20and%20dendritic%20self-organization.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_dendrite')

% --- Executes on button press in pushbutton104.
function pushbutton104_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Free-energy%20and%20illusions%20the%20Cornsweet%20effect.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_Cornsweet')

% --- Executes on button press in pushbutton105.
function pushbutton105_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Smooth%20Pursuit%20and%20Visual%20Occlusion.pdf';
run_demo_Callback(hObject, handles, 'ADEM_pursuit')

% --- Executes on button press in pushbutton117.
function pushbutton117_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dopamine%20Affordance%20and%20Active%20Inference.pdf';
run_demo_Callback(hObject, handles, 'ADEM_cued_response')

% --- Executes on button press in pushbutton118.
function pushbutton118_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Attractors%20in%20song.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MMN_deviance')

% --- Executes on button press in pushbutton119.
function pushbutton119_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Bayesian%20State%20Estimation%20Using%20Generalized%20Coordinates.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_contact_lens')

% --- Executes on button press in pushbutton120.
function pushbutton120_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Bayesian%20State%20Estimation%20Using%20Generalized%20Coordinates.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_GF_and_KF')

% --- Executes on button press in pushbutton121.
function pushbutton121_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20agency%20optimal%20control%20without%20cost%20functions.pdf';
run_demo_Callback(hObject, handles, 'spm_MDP_mountain_car')

% --- Executes on button press in pushbutton122.
function pushbutton122_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Perceptions%20as%20hypotheses%20saccades%20as%20experiments.pdf';
run_demo_Callback(hObject, handles, 'ADEM_salience')

% --- Executes on button press in pushbutton123.
function pushbutton123_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Perception%20and%20self-organized%20instability.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_SOC')

% --- Executes on button press in pushbutton124.
function pushbutton124_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Smooth%20Pursuit%20and%20Visual%20Occlusion.pdf';
run_demo_Callback(hObject, handles, 'ADEM_occulomotor_delays')

% --- Executes on button press in pushbutton125.
function pushbutton125_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Smooth%20Pursuit%20and%20Visual%20Occlusion.pdf';
run_demo_Callback(hObject, handles, 'ADEM_occlusion')

% --- Executes on button press in pushbutton126.
function pushbutton126_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Generalised%20Filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_ALAP')

% --- Executes on button press in pushbutton127.
function pushbutton127_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20sensory%20attenuation%20and%20illusions.pdf';
run_demo_Callback(hObject, handles, 'ALAP_demo_attenuation')

% --- Executes on button press in pushbutton128.
function pushbutton128_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Observing%20the%20Observer%20I.pdf';
run_demo_Callback(hObject, handles, 'spm_meta_model')

% --- Executes on button press in pushbutton129.
function pushbutton129_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_evidence_accumulation')

% --- Executes on button press in pushbutton130.
function pushbutton130_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20agency%20optimal%20control%20without%20cost%20functions.pdf';
run_demo_Callback(hObject, handles, 'spm_MDP_offer')

% --- Executes on button press in pushbutton132.
function pushbutton132_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Life%20as%20we%20know%20it.pdf';
run_demo_Callback(hObject, handles, 'FEP_Manifold')

% --- Executes on button press in pushbutton133.
function pushbutton133_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/The%20anatomy%20of%20choice%20active%20inference%20and%20agency.pdf';
run_demo_Callback(hObject, handles, 'spm_MDP_trust')

% --- Executes on button press in pushbutton134.
function pushbutton134_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_induced_fMRI')

% --- Executes on button press in pushbutton142.
function pushbutton142_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/The%20anatomy%20of%20choice%20active%20inference%20and%20agency.pdf';
run_demo_Callback(hObject, handles, 'spm_MDP_urn')

% --- Executes on button press in pushbutton143.
function pushbutton143_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_large_fMRI')

% --- Executes on button press in pushbutton144.
function pushbutton144_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_connectivity_fMRI')

% --- Executes on button press in pushbutton145.
function pushbutton145_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_modes_fMRI')

% --- Executes on button press in pushbutton146.
function pushbutton146_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20modeling%20of%20neuronal%20responses%20in%20fMRI%20using%20cubature%20Kalman%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm_SCK')

% --- Executes on button press in pushbutton147.
function pushbutton147_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Generalised%20Filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_convolution_LAP')

% --- Executes on button press in pushbutton148.
function pushbutton148_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20modeling%20of%20neuronal%20responses%20in%20fMRI%20using%20cubature%20Kalman%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_lorenz_LAP')

% --- Executes on button press in pushbutton149.
function pushbutton149_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Dynamic%20modeling%20of%20neuronal%20responses%20in%20fMRI%20using%20cubature%20Kalman%20filtering.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_doublewell_LAP')

% --- Executes on button press in pushbutton150.
function pushbutton150_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Attention%20uncertainty%20and%20free-energy.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_texture')

% --- Executes on button press in pushbutton151.
function pushbutton151_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Predictive%20coding%20under%20the%20free-energy%20principle.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_duet')

% --- Executes on button press in pushbutton152.
function pushbutton152_Callback(hObject, eventdata, handles)
% --- Executes on button press in pushbutton133.
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/The%20anatomy%20of%20choice%20active%20inference%20and%20agency.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_maze')

% --- Executes on button press in pushbutton153.
function pushbutton153_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Post%20hoc%20Bayesian%20model%20selection.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_Bayesian_Model_Reduction')

% --- Executes on button press in pushbutton154.
function pushbutton154_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Life%20as%20we%20know%20it.pdf';
run_demo_Callback(hObject, handles, 'DEM_morphogenesis')

% --- Executes on button press in pushbutton155.
function pushbutton155_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Policies%20and%20Priors.pdf';
run_demo_Callback(hObject, handles, 'ADEM_eyeblink')

% --- Executes on button press in pushbutton156.
function pushbutton156_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Post%20hoc%20Bayesian%20model%20selection.pdf';
run_demo_Callback(hObject, handles, 'DEMO_SLR')

% --- Executes on button press in pushbutton157.
function pushbutton157_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Post%20hoc%20Bayesian%20model%20selection.pdf';
run_demo_Callback(hObject, handles, 'DEMO_BMR_PEB')

% --- Executes on button press in pushbutton158.
function pushbutton158_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Post%20hoc%20Bayesian%20model%20selection.pdf';
run_demo_Callback(hObject, handles, 'DEMO_GROUP_PEB')

% --- Executes on button press in pushbutton159.
function pushbutton159_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/spm/doc/papers/sjk_aibf.pdf';
run_demo_Callback(hObject, handles, 'DEM_spatial_deconvolution')

% --- Executes on button press in pushbutton160.
function pushbutton160_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_ontology')

% --- Executes on button press in pushbutton161.
function pushbutton161_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_habits')

% --- Executes on button press in pushbutton162.
function pushbutton162_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_fit')

% --- Executes on button press in pushbutton163.
function pushbutton163_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_X')

% --- Executes on button press in pushbutton168.
function pushbutton168_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_search')

% --- Executes on button press in pushbutton169.
function pushbutton169_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_reading')

% --- Executes on button press in pushbutton170.
function pushbutton170_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_rule')

% --- Executes on button press in pushbutton197.
function pushbutton197_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_BMR_PEB')

% --- Executes on button press in pushbutton198.
function pushbutton198_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_GROUP_PEB')

% --- Executes on button press in pushbutton199.
function pushbutton199_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_DCM_PEB')

% --- Executes on button press in pushbutton200.
function pushbutton200_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_DCM_PEB_FIT')

% --- Executes on button press in pushbutton201.
function pushbutton201_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_BAYES_FACTORS')

% --- Executes on button press in pushbutton202.
function pushbutton202_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_Lindley_paradox')

% --- Executes on button press in pushbutton203.
function pushbutton203_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_BMR_PEB')

% --- Executes on button press in pushbutton210.
function pushbutton210_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_fMRI_PEB')

% --- Executes on button press in pushbutton211.
function pushbutton211_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_MDP_decision')


% --- Executes on button press in pushbutton212.
function pushbutton212_Callback(hObject, eventdata, handles)
handles.web = 'http://www.fil.ion.ucl.ac.uk/~karl/Active%20inference%20and%20epistemic%20value.pdf';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_DEM')
