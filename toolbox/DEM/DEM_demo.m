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

% Last Modified by GUIDE v2.5 02-Jan-2020 19:35:25

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
handles.web    = 'The free-energy principle A unified brain theory';

% Display PDF image
axes5_CreateFcn(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
ax = handles.axes5;
imagesc(imread('PDF.jpg'),'Parent',ax), axis(ax,'off');


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
try, 
    url = strcat('http://www.fil.ion.ucl.ac.uk/~karl/',handles.web);
    web(url,'-browser','-notoolbar'); 
end


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
handles.web = 'Variational free energy and the Laplace approximation';
run_demo_Callback(hObject, handles, 'DEM_demo_GLM')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
handles.web = 'Variational free energy and the Laplace approximation';
run_demo_Callback(hObject, handles, 'DEM_demo_PEB')

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
handles.web = 'Hierarchical Models in the Brain';
run_demo_Callback(hObject, handles, 'DEM_demo_factor_analysis')

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
handles.web = 'DEM A variational treatment of dynamic systems';
run_demo_Callback(hObject, handles, 'DEM_demo_OU')

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
handles.web = 'Hierarchical Models in the Brain';
run_demo_Callback(hObject, handles, 'DEM_demo_convolution')

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
handles.web = 'DEM A variational treatment of dynamic systems';
run_demo_Callback(hObject, handles, 'DEM_demo_EM')

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
handles.web = 'DEM A variational treatment of dynamic systems';
run_demo_Callback(hObject, handles, 'DEM_demo_DEM')

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
handles.web = 'DEM A variational treatment of dynamic systems';
run_demo_Callback(hObject, handles, 'DEM_demo_filtering')

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
handles.web = 'DEM A variational treatment of dynamic systems';
run_demo_Callback(hObject, handles, 'DEM_demo_Lorenz')

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic causal modelling';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm')

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
handles.web = 'DEM A variational treatment of dynamic systems';
run_demo_Callback(hObject, handles, 'DEM_demo_double_well')

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
handles.web = 'Variational filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_DFP')

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
handles.web = 'Variational filtering';
run_demo_Callback(hObject, handles, 'DFP_demo_double_well')

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
handles.web = 'Cortical circuits for perceptual inference';
run_demo_Callback(hObject, handles, 'DEM_demo_song_priors')

% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
handles.web = 'Predictive coding under the free-energy principle';
run_demo_Callback(hObject, handles, 'DEM_demo_song_inference')

% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
handles.web = 'Cortical circuits for perceptual inference';
run_demo_Callback(hObject, handles, 'DEM_demo_song_omission')

% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
handles.web = 'Learning and inference in the brain';
run_demo_Callback(hObject, handles, 'DEM_demo_face_inference')

% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
handles.web = 'Attractors in song';
run_demo_Callback(hObject, handles, 'DEM_demo_MMN')

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEM_demo_Gabor')

% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
handles.web = 'Action and behavior A free-energy formulation';
run_demo_Callback(hObject, handles, 'ADEM_visual')

% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
handles.web = 'Action and behavior A free-energy formulation';
run_demo_Callback(hObject, handles, 'ADEM_motor')

% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
handles.web = 'Reinforcement Learning or Active Inference';
run_demo_Callback(hObject, handles, 'ADEM_learning')

% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
handles.web = 'Action and behavior A free-energy formulation';
run_demo_Callback(hObject, handles, 'ADEM_lorenz')

% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
handles.web = 'Action and behavior A free-energy formulation';
run_demo_Callback(hObject, handles, 'ADEM_reaching')

% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
handles.web = 'Free Energy Value and Attractors';
run_demo_Callback(hObject, handles, 'ADEM_lorenz_entropy')

% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
handles.web = 'Free Energy Value and Attractors';
run_demo_Callback(hObject, handles, 'ADEM_mountaincar_loss')

% --- Executes on button press in pushbutton80.
function pushbutton80_Callback(hObject, eventdata, handles)
handles.web = 'Policies and Priors';
run_demo_Callback(hObject, handles, 'ADEM_SHC_demo')

% --- Executes on button press in pushbutton81.
function pushbutton81_Callback(hObject, eventdata, handles)
handles.web = 'Free Energy Value and Attractors';
run_demo_Callback(hObject, handles, 'ADEM_lorenz_surprise')

% --- Executes on button press in pushbutton83.
function pushbutton83_Callback(hObject, eventdata, handles)
handles.web = 'Generalised Filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_LAP')

% --- Executes on button press in pushbutton84.
function pushbutton84_Callback(hObject, eventdata, handles)
handles.web = 'Generalised Filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm_LAP')

% --- Executes on button press in pushbutton91.
function pushbutton91_Callback(hObject, eventdata, handles)
handles.web = 'Attention uncertainty and free-energy';
run_demo_Callback(hObject, handles, 'DEM_demo_Posner')

% --- Executes on button press in pushbutton92.
function pushbutton92_Callback(hObject, eventdata, handles)
handles.web = 'Policies and Priors';
run_demo_Callback(hObject, handles, 'ADEM_cost_SHC')

% --- Executes on button press in pushbutton94.
function pushbutton94_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, '%%%% DEM_demo_CM_Lorenz %%%')

% --- Executes on button press in pushbutton96.
function pushbutton96_Callback(hObject, eventdata, handles)
handles.web = 'Action understanding and active inference';
run_demo_Callback(hObject, handles, 'ADEM_observe')

% --- Executes on button press in pushbutton97.
function pushbutton97_Callback(hObject, eventdata, handles)
handles.web = 'Network discovery with DCM';
run_demo_Callback(hObject, handles, 'DEM_demo_DCM_LAP')

% --- Executes on button press in pushbutton98.
function pushbutton98_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic modeling of neuronal responses in fMRI using cubature Kalman filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_convolution_LAP')

% --- Executes on button press in pushbutton99.
function pushbutton99_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic modeling of neuronal responses in fMRI using cubature Kalman filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_lorenz_LAP')

% --- Executes on button press in pushbutton100.
function pushbutton100_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic modeling of neuronal responses in fMRI using cubature Kalman filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_doublewell_LAP')

% --- Executes on button press in pushbutton101.
function pushbutton101_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic modeling of neuronal responses in fMRI using cubature Kalman filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm_SCK')

% --- Executes on button press in pushbutton102.
function pushbutton102_Callback(hObject, eventdata, handles)
handles.web = 'Action understanding and active inference';
run_demo_Callback(hObject, handles, 'ADEM_writing')

% --- Executes on button press in pushbutton103.
function pushbutton103_Callback(hObject, eventdata, handles)
handles.web = 'Free energy and dendritic self-organization';
run_demo_Callback(hObject, handles, 'DEM_demo_dendrite')

% --- Executes on button press in pushbutton104.
function pushbutton104_Callback(hObject, eventdata, handles)
handles.web = 'Free-energy and illusions the Cornsweet effect';
run_demo_Callback(hObject, handles, 'DEM_demo_Cornsweet')

% --- Executes on button press in pushbutton105.
function pushbutton105_Callback(hObject, eventdata, handles)
handles.web = 'Smooth Pursuit and Visual Occlusion';
run_demo_Callback(hObject, handles, 'ADEM_pursuit')

% --- Executes on button press in pushbutton117.
function pushbutton117_Callback(hObject, eventdata, handles)
handles.web = 'Dopamine Affordance and Active Inference';
run_demo_Callback(hObject, handles, 'ADEM_cued_response')

% --- Executes on button press in pushbutton118.
function pushbutton118_Callback(hObject, eventdata, handles)
handles.web = 'Attractors in song';
run_demo_Callback(hObject, handles, 'DEM_demo_MMN_deviance')

% --- Executes on button press in pushbutton119.
function pushbutton119_Callback(hObject, eventdata, handles)
handles.web = 'Bayesian State Estimation Using Generalized Coordinates';
run_demo_Callback(hObject, handles, 'DEM_demo_contact_lens')

% --- Executes on button press in pushbutton120.
function pushbutton120_Callback(hObject, eventdata, handles)
handles.web = 'Bayesian State Estimation Using Generalized Coordinates';
run_demo_Callback(hObject, handles, 'DEM_demo_GF_and_KF')

% --- Executes on button press in pushbutton121.
function pushbutton121_Callback(hObject, eventdata, handles)
handles.web = 'Active inference and agency optimal control without cost functions';
run_demo_Callback(hObject, handles, 'spm_MDP_mountain_car')

% --- Executes on button press in pushbutton122.
function pushbutton122_Callback(hObject, eventdata, handles)
handles.web = 'Perceptions as hypotheses saccades as experiments';
run_demo_Callback(hObject, handles, 'ADEM_salience')

% --- Executes on button press in pushbutton123.
function pushbutton123_Callback(hObject, eventdata, handles)
handles.web = 'Perception and self-organized instability';
run_demo_Callback(hObject, handles, 'DEM_demo_SOC')

% --- Executes on button press in pushbutton124.
function pushbutton124_Callback(hObject, eventdata, handles)
handles.web = 'Active inference, eye movements and oculomotor delays';
run_demo_Callback(hObject, handles, 'ADEM_occulomotor_delays')

% --- Executes on button press in pushbutton125.
function pushbutton125_Callback(hObject, eventdata, handles)
handles.web = 'Smooth Pursuit and Visual Occlusion';
run_demo_Callback(hObject, handles, 'ADEM_occlusion')

% --- Executes on button press in pushbutton126.
function pushbutton126_Callback(hObject, eventdata, handles)
handles.web = 'Generalised Filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_ALAP')

% --- Executes on button press in pushbutton127.
function pushbutton127_Callback(hObject, eventdata, handles)
handles.web = 'Active inference sensory attenuation and illusions';
run_demo_Callback(hObject, handles, 'ALAP_demo_attenuation')

% --- Executes on button press in pushbutton128.
function pushbutton128_Callback(hObject, eventdata, handles)
handles.web = 'Observing the Observer I';
run_demo_Callback(hObject, handles, 'spm_meta_model')

% --- Executes on button press in pushbutton129.
function pushbutton129_Callback(hObject, eventdata, handles)
handles.web = 'Predictive Coding or Evidence Accumulation';
run_demo_Callback(hObject, handles, 'DEM_evidence_accumulation')

% --- Executes on button press in pushbutton130.
function pushbutton130_Callback(hObject, eventdata, handles)
handles.web = 'The anatomy of choice active inference and agency';
run_demo_Callback(hObject, handles, 'spm_MDP_offer')

% --- Executes on button press in pushbutton132.
function pushbutton132_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'FEP_Manifold')

% --- Executes on button press in pushbutton133.
function pushbutton133_Callback(hObject, eventdata, handles)
handles.web = 'The anatomy of choice active inference and agency';
run_demo_Callback(hObject, handles, 'spm_MDP_trust')

% --- Executes on button press in pushbutton134.
function pushbutton134_Callback(hObject, eventdata, handles)
handles.web = 'A DCM for resting state fMRI';
run_demo_Callback(hObject, handles, 'DEM_demo_induced_fMRI')

% --- Executes on button press in pushbutton142.
function pushbutton142_Callback(hObject, eventdata, handles)
handles.web = 'Active Inference, Evidence Accumulation, and the Urn Task';
run_demo_Callback(hObject, handles, 'spm_MDP_urn')

% --- Executes on button press in pushbutton143.
function pushbutton143_Callback(hObject, eventdata, handles)
handles.web = 'Network discovery with large DCMs';
run_demo_Callback(hObject, handles, 'DEM_demo_large_fMRI')

% --- Executes on button press in pushbutton144.
function pushbutton144_Callback(hObject, eventdata, handles)
handles.web = 'On nodes and modes in resting state fMRI';
run_demo_Callback(hObject, handles, 'DEM_demo_connectivity_fMRI')

% --- Executes on button press in pushbutton145.
function pushbutton145_Callback(hObject, eventdata, handles)
handles.web = 'On nodes and modes in resting state fMRI';
run_demo_Callback(hObject, handles, 'DEM_demo_modes_fMRI')

% --- Executes on button press in pushbutton146.
function pushbutton146_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic modeling of neuronal responses in fMRI using cubature Kalman filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_hdm_SCK')

% --- Executes on button press in pushbutton147.
function pushbutton147_Callback(hObject, eventdata, handles)
handles.web = 'Generalised Filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_convolution_LAP')

% --- Executes on button press in pushbutton148.
function pushbutton148_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic modeling of neuronal responses in fMRI using cubature Kalman filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_lorenz_LAP')

% --- Executes on button press in pushbutton149.
function pushbutton149_Callback(hObject, eventdata, handles)
handles.web = 'Dynamic modeling of neuronal responses in fMRI using cubature Kalman filtering';
run_demo_Callback(hObject, handles, 'DEM_demo_doublewell_LAP')

% --- Executes on button press in pushbutton150.
function pushbutton150_Callback(hObject, eventdata, handles)
handles.web = 'Cerebral Hierarchies, predictive processing, precision and the pulvinar';
run_demo_Callback(hObject, handles, 'DEM_demo_texture')

% --- Executes on button press in pushbutton151.
function pushbutton151_Callback(hObject, eventdata, handles)
handles.web = 'A Duet for one';
run_demo_Callback(hObject, handles, 'DEM_demo_duet')

% --- Executes on button press in pushbutton152.
function pushbutton152_Callback(hObject, eventdata, handles)
% --- Executes on button press in pushbutton133.
handles.web = 'Active inference and epistemic value';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_maze')

% --- Executes on button press in pushbutton153.
function pushbutton153_Callback(hObject, eventdata, handles)
handles.web = 'Post hoc Bayesian model selection';
run_demo_Callback(hObject, handles, 'DEM_demo_Bayesian_Model_Reduction')

% --- Executes on button press in pushbutton154.
function pushbutton154_Callback(hObject, eventdata, handles)
handles.web = 'Knowing ones place';
run_demo_Callback(hObject, handles, 'DEM_morphogenesis')

% --- Executes on button press in pushbutton155.
function pushbutton155_Callback(hObject, eventdata, handles)
handles.web = 'Active Inference and Learning in the Cerebellum';
run_demo_Callback(hObject, handles, 'ADEM_eyeblink')

% --- Executes on button press in pushbutton156.
function pushbutton156_Callback(hObject, eventdata, handles)
handles.web = 'Post hoc Bayesian model selection';
run_demo_Callback(hObject, handles, 'DEMO_SLR')

% --- Executes on button press in pushbutton157.
function pushbutton157_Callback(hObject, eventdata, handles)
handles.web = 'Post hoc Bayesian model selection';
run_demo_Callback(hObject, handles, 'DEMO_BMR_PEB')

% --- Executes on button press in pushbutton158.
function pushbutton158_Callback(hObject, eventdata, handles)
handles.web = 'Bayesian model reduction and empirical Bayes for group (DCM) studies';
run_demo_Callback(hObject, handles, 'DEMO_GROUP_PEB')

% --- Executes on button press in pushbutton159.
function pushbutton159_Callback(hObject, eventdata, handles)
handles.web = 'Anatomically informed basis functions';
run_demo_Callback(hObject, handles, 'DEM_spatial_deconvolution')

% --- Executes on button press in pushbutton160.
function pushbutton160_Callback(hObject, eventdata, handles)
handles.web = 'Computational Nosology';
run_demo_Callback(hObject, handles, 'DEM_demo_ontology')

% --- Executes on button press in pushbutton161.
function pushbutton161_Callback(hObject, eventdata, handles)
handles.web = 'Active inference and learning';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_habits')

% --- Executes on button press in pushbutton162.
function pushbutton162_Callback(hObject, eventdata, handles)
handles.web = 'Computational Phenotyping in Psychiatry';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_fit')

% --- Executes on button press in pushbutton163.
function pushbutton163_Callback(hObject, eventdata, handles)
handles.web = 'Active Inference A Process Theory';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_X')

% --- Executes on button press in pushbutton168.
function pushbutton168_Callback(hObject, eventdata, handles)
handles.web = 'Scene Construction,Visual Foraging and Active Inference';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_search')

% --- Executes on button press in pushbutton169.
function pushbutton169_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_reading')

% --- Executes on button press in pushbutton170.
function pushbutton170_Callback(hObject, eventdata, handles)
handles.web = 'Active Inference Curiosity and Insight';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_rule')

% --- Executes on button press in pushbutton197.
function pushbutton197_Callback(hObject, eventdata, handles)
handles.web = 'Empirical Bayes for DCM A Group Inversion Scheme';
run_demo_Callback(hObject, handles, 'DEMO_BMR_PEB')

% --- Executes on button press in pushbutton198.
function pushbutton198_Callback(hObject, eventdata, handles)
handles.web = 'Empirical Bayes for DCM A Group Inversion Scheme';
run_demo_Callback(hObject, handles, 'DEMO_GROUP_PEB')

% --- Executes on button press in pushbutton199.
function pushbutton199_Callback(hObject, eventdata, handles)
handles.web = 'Bayesian model reduction and empirical Bayes for group (DCM) studies';
run_demo_Callback(hObject, handles, 'DEMO_DCM_PEB')

% --- Executes on button press in pushbutton200.
function pushbutton200_Callback(hObject, eventdata, handles)
handles.web = 'Bayesian model reduction and empirical Bayes for group (DCM) studies';
run_demo_Callback(hObject, handles, 'DEMO_DCM_PEB_FIT')

% --- Executes on button press in pushbutton201.
function pushbutton201_Callback(hObject, eventdata, handles)
handles.web = 'Comparing dynamic causal models';
run_demo_Callback(hObject, handles, 'DEMO_BAYES_FACTORS')

% --- Executes on button press in pushbutton202.
function pushbutton202_Callback(hObject, eventdata, handles)
handles.web = '';
run_demo_Callback(hObject, handles, 'DEMO_Lindley_paradox')

% --- Executes on button press in pushbutton203.
function pushbutton203_Callback(hObject, eventdata, handles)
handles.web = 'Bayesian model selection for group studies';
run_demo_Callback(hObject, handles, 'DEMO_BMR_PEB')

% --- Executes on button press in pushbutton210.
function pushbutton210_Callback(hObject, eventdata, handles)
handles.web = 'Bayesian model reduction and empirical Bayes for group (DCM) studies';
run_demo_Callback(hObject, handles, 'DEM_demo_fMRI_PEB')

% --- Executes on button press in pushbutton211.
function pushbutton211_Callback(hObject, eventdata, handles)
handles.web = 'Active inference and epistemic value';
run_demo_Callback(hObject, handles, 'DEM_MDP_decision')

% --- Executes on button press in pushbutton212.
function pushbutton212_Callback(hObject, eventdata, handles)
handles.web = 'The graphical brain';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_DEM')

% --- Executes on button press in pushbutton213.
function pushbutton213_Callback(hObject, eventdata, handles)
handles.web = 'Computational Phenotyping in Psychiatry';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_fit_fields')

% --- Executes on button press in pushbutton214.
function pushbutton214_Callback(hObject, eventdata, handles)
handles.web = 'Planning and navigation as active inference';
run_demo_Callback(hObject, handles, 'DEMO_MDP_maze')

% --- Executes on button press in pushbutton215.
function pushbutton215_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'FEP_fluctuations')

% --- Executes on button press in pushbutton216.
function pushbutton216_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'FEP_MB_demo')

% --- Executes on button press in pushbutton217.
function pushbutton217_Callback(hObject, eventdata, handles)
handles.web = 'Knowing ones place';
run_demo_Callback(hObject, handles, 'DEM_cells_cells')

% --- Executes on button press in pushbutton218.
function pushbutton218_Callback(hObject, eventdata, handles)
handles.web = 'Knowing ones place';
run_demo_Callback(hObject, handles, 'DEM_cells')

% --- Executes on button press in pushbutton219.
function pushbutton219_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'DEM_self_entropy')

% --- Executes on button press in pushbutton223.
function pushbutton223_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'FEP_self_entropy')

% --- Executes on button press in pushbutton221.
function pushbutton221_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'DEM_HB_and_LE')

% --- Executes on button press in pushbutton222.
function pushbutton222_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'FEP_physics')

% --- Executes on button press in pushbutton224.
function pushbutton224_Callback(hObject, eventdata, handles)
handles.web = 'Post hoc Bayesian model selection';
run_demo_Callback(hObject, handles, 'DEMO_AI_NLSI')

% --- Executes on button press in pushbutton225.
function pushbutton225_Callback(hObject, eventdata, handles)
handles.web = 'Active inference and the anatomy of oculomotion';
run_demo_Callback(hObject, handles, 'MDP_DEM_Oculomotion_demo')

% --- Executes on button press in pushbutton226.
function pushbutton226_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'DEMO_MDP_questions')

% --- Executes on button press in pushbutton227.
function pushbutton227_Callback(hObject, eventdata, handles)
handles.web = 'A Multivariate Analysis of PET Activation Studies';
run_demo_Callback(hObject, handles, 'DEMO_CVA_RSA')

% --- Executes on button press in pushbutton228.
function pushbutton228_Callback(hObject, eventdata, handles)
handles.web = 'Planning and navigation as active inference';
run_demo_Callback(hObject, handles, 'DEMO_niche_construction')

% --- Executes on button press in pushbutton229.
function pushbutton229_Callback(hObject, eventdata, handles)
handles.web = 'Generalised Filtering';
run_demo_Callback(hObject, handles, 'DEM_coupled_oscillators')

% --- Executes on button press in pushbutton230.
function pushbutton230_Callback(hObject, eventdata, handles)
handles.web = 'Generalised Filtering';
run_demo_Callback(hObject, handles, 'KLDemo')

% --- Executes on button press in pushbutton231.
function pushbutton231_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'DEM_FEP_Lorenz')

% --- Executes on button press in pushbutton232.
function pushbutton232_Callback(hObject, eventdata, handles)
handles.web = 'Life as we know it';
run_demo_Callback(hObject, handles, 'DEM_FEP_Least_Action')

% --- Executes on button press in pushbutton233.
function pushbutton233_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'spm_voice')

% --- Executes on button press in pushbutton234.
function pushbutton234_Callback(hObject, eventdata, handles)
handles.web = 'Interoceptive inference';
run_demo_Callback(hObject, handles, 'MDP_Heart_Beat')

% --- Executes on button press in pushbutton239.
function pushbutton239_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'spm_voice')

% --- Executes on button press in pushbutton240.
function pushbutton240_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'spm_voice_P300')

% --- Executes on button press in pushbutton241.
function pushbutton241_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'spm_voice_repeat')

% --- Executes on button press in pushbutton242.
function pushbutton242_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'spm_voice_read')

% --- Executes on button press in pushbutton243.
function pushbutton243_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'DEMO_MDP_voice')

% --- Executes on button press in pushbutton244.
function pushbutton244_Callback(hObject, eventdata, handles)
handles.web = 'Deep temporal models and active inference';
run_demo_Callback(hObject, handles, 'MDP_DEM_Mixed_Models_Movement')

% --- Executes on button press in pushbutton245.
function pushbutton245_Callback(hObject, eventdata, handles)
handles.web = 'A free energy principle for a particular physics';
run_demo_Callback(hObject, handles, 'DEMO_DCM_MB')

% --- Executes on button press in pushbutton247.
function pushbutton247_Callback(hObject, eventdata, handles)
handles.web = 'Active Inference A Process Theory';
run_demo_Callback(hObject, handles, 'DEM_demo_MDP_XX')

% --- Executes on button press in pushbutton248.
function pushbutton248_Callback(hObject, eventdata, handles)
handles.web = 'Planning and navigation as active inference';
run_demo_Callback(hObject, handles, 'DEMO_MDP_maze_X')
