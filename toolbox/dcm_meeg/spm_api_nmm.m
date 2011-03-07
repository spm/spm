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
 
% Last Modified by GUIDE v2.5 17-Oct-2008 16:14:27
 
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
%--------------------------------------------------------------------------
function spm_api_nmm_OpeningFcn(hObject, eventdata, handles, varargin)
set(gcf,'name','Neural mass modelling: review of priors');
try
    handles.DCM = varargin{1};
    unpack_Callback(hObject, eventdata, handles);
end
 
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
%--------------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)
 
[f,p]       = uigetfile('*.mat','please select DCM file'); cd(p)
name        = fullfile(p,f);
DCM         = load(name,'-mat');
DCM         = DCM.DCM;
DCM.name    = name;
handles.DCM = DCM;
set(handles.name,'String',f);
 
% display priors
%--------------------------------------------------------------------------
unpack_Callback(hObject, eventdata, handles);
guidata(hObject,handles);
 
 
% --- Executes on button press in save.
%--------------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
try
    [p,file] = fileparts(handles.DCM.name);
catch
    try
        [p,file] = fileparts(handles.DCM.xY.Dfile);
        file     = ['DCM_' file];
    catch
        file     = ['DCM'  date];
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
%--------------------------------------------------------------------------
function kernels_Callback(hObject, eventdata, handles)
clear spm_erp_L spm_gen_erp
DCM   = handles.DCM;
pst   = str2num(get(handles.pst,'String'));
Hz    = str2num(get(handles.Hz, 'String'));
 
% initialise states
%--------------------------------------------------------------------------
M     = DCM.M;
model = DCM.options.model;
[x f] = spm_dcm_x_neural(M.pE,model);
M.x   = x;
M.f   = f;
M.m   = size(M.pE.C,2);
M.g   = {};
M.ons = 32*ones(M.m,1);

    
 
% Volterra Kernels
%==========================================================================
 
% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L1,L2] = spm_bireduce(M,M.pE);
 
% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
M.ns    = 128;
U.dt    = pst/M.ns/1000;
t       = [1:M.ns]*U.dt*1000;
[K0,K1] = spm_kernels(M0,M1,L1,L2,M.ns,U.dt);
 
axes(handles.kernel)
plot(t,K1(:,:,1))
set(gca,'XLim',[min(t) max(t)])
ylabel('1st-order Volterra kernel')
xlabel('time (ms)')
 
% transfer function (source space)
%--------------------------------------------------------------------------
for i = 1:size(K1,2)
    g(:,i) = fft(K1(:,i));
    g(:,i) = abs(g(:,i).*conj(g(:,i)));
end
i   = [1:M.ns/2];
g   = g(i + 1,:);
w   = i/(M.ns*U.dt);
 
axes(handles.transfer)
plot(w,g)
set(gca,'XLim',[1 Hz])
xlabel('frequency time (Hz)')
ylabel('spectral density')
 
% prediction (source space)
%--------------------------------------------------------------------------
x     = spm_gen_erp(M.pE,M,U);
x     = x{1};
 
axes(handles.erp)
plot(t,x)
set(gca,'XLim',[min(t) max(t)])
xlabel('peristimulus time (ms)')
ylabel('voltage and conductance')
 
% fft (source space)
%--------------------------------------------------------------------------
clear g
for i = 1:size(x,2)
    g(:,i) = fft(x(:,i));
    g(:,i) = abs(g(:,i).*conj(g(:,i)));
end
i   = [1:M.ns/2];
g   = g(i + 1,:);
w   = i/(M.ns*U.dt);
 
axes(handles.fft)
plot(w,g)
set(gca,'XLim',[1 Hz])
xlabel('frequency time (Hz)')
ylabel('spectral density')
 
 
% unpack priors and display
%==========================================================================
function unpack_Callback(hObject, eventdata, handles);
 
% clear previous objects
%--------------------------------------------------------------------------
h = get(gcf,'Children');
for i = 1:length(h)
    if strcmp(get(h(i),'Tag'),'tmp')
        delete(h(i));
    end
end

% Type of model (neuronal)
%--------------------------------------------------------------------------
% 'ERP'    - (linear second order NMM slow)
% 'SEP'    - (linear second order NMM fast)
% 'CMC'    - (linear second order NMM Canonical microcircuit)
% 'LFP'    - (linear second order NMM self-inhibition)
% 'NMM'    - (nonlinear second order NMM first-order moments)
% 'MFM'    - (nonlinear second order NMM second-order moments)
try
    model = handles.DCM.options.model;
catch
    model = get(handles.model,'String');
    model = model{get(handles.model,'Value')};
    DCM.options.model = model;
end
switch model
    case{'ERP'}, set(handles.model,'Value',1);
    case{'SEP'}, set(handles.model,'Value',2);
    case{'CMC'}, set(handles.model,'Value',3);
    case{'LFP'}, set(handles.model,'Value',4);
    case{'NMM'}, set(handles.model,'Value',5);
    case{'MFM'}, set(handles.model,'Value',6);
    otherwise
end

% get priors
%--------------------------------------------------------------------------
try
    handles.DCM.M.pE;
catch
    handles = reset_Callback(hObject, eventdata, handles);
end
pE = handles.DCM.M.pE;
pC = handles.DCM.M.pC;
try
    pC = spm_unvec(diag(pC),pE);
end

% display connection switches later
%--------------------------------------------------------------------------
try, pE = rmfield(pE,'B');,  end
try, pE = rmfield(pE,'SA');, end
try, pE = rmfield(pE,'GE');, end
try, pE = rmfield(pE,'GI');, end

% do not display spatial and spectral priors
%--------------------------------------------------------------------------
try, pE = rmfield(pE,{'Lpos','L','J'});,      end
try, pE = rmfield(pE,{'a','b','c','d'});, end

 
% display fields
%--------------------------------------------------------------------------
color = {[1 1 1],get(handles.name,'BackgroundColor')};
 
x0    = 1/8;
y0    = 1 - 1/8;
sx    = 1/24;
sy    = 1/48;
dx    = 1/256 + sx;
dy    = 1/256 + sy;
F     = fieldnames(pE);
for f = 1:length(F)
    
    P = getfield(pE,F{f});
 
    if iscell(P)
 
        % cell
        %------------------------------------------------------------------
        for k = 1:length(P)
            for i = 1:size(P{k},1);
                for j = 1:size(P{k},2)
                    x   = x0 + j*dx + (size(P{k},2) + 1)*dx*(k - 1);
                    y   = y0 - (i - 1)*dy;
                    str = sprintf('handles.DCM.M.pE.%s{%i}(%i,%i)',F{f},k,i,j);
                    str = ['handles=guidata(gcbo);' str '=str2num(get(gcbo,''String''));guidata(gcbo,handles);'];
                    h = uicontrol(gcf,...
                        'Units','Normalized',...
                        'Position',[x y sx sy],...
                        'Style','edit',...
                        'String',sprintf('%-4.2f',P{k}(i,j)),...
                        'enable','on',...
                        'Tag','tmp',...
                        'BackgroundColor',color{1 + ~eval(sprintf('pC.%s{%i}(%i,%i)',F{f},k,i,j))},...
                        'Callback',str);
                end
            end
        end
 
 
    elseif isvector(P)
        
        % vector
        %------------------------------------------------------------------
        for i   = 1:length(P);
            x   = x0 + i*dx;
            y   = y0;
            str = sprintf('handles.DCM.M.pE.%s(%i)',F{f},i);
            str = ['handles=guidata(gcbo);' str '=str2num(get(gcbo,''String''));guidata(gcbo,handles);'];
            uicontrol(gcf,...
                'Units','Normalized',...
                'Position',[x y sx sy],...
                'Style','edit',...
                'String',sprintf('%-4.2f',full(P(i))),...
                'enable','on',...
                'Tag','tmp',...
                'BackgroundColor',color{1 + ~eval(sprintf('pC.%s(%i)',F{f},i))},...
                'Callback',str);
        end
    else
        
        % matrix
        %------------------------------------------------------------------
        for i = 1:size(P,1);
            for j   = 1:size(P,2)
                x   = x0 + j*dx;
                y   = y0 - (i - 1)*dy;
                str = sprintf('handles.DCM.M.pE.%s(%i,%i)',F{f},i,j);
                str = ['handles=guidata(gcbo);' str '=str2num(get(gcbo,''String''));guidata(gcbo,handles);'];
                uicontrol(gcf,...
                    'Units','Normalized',...
                    'Position',[x y sx sy],...
                    'Style','edit',...
                    'String',sprintf('%-4.2f',full(P(i,j))),...
                    'enable','on',...
                    'Tag','tmp',...
                    'BackgroundColor',color{1 + ~eval(sprintf('pC.%s(%i,%i)',F{f},i,j))},...
                    'Callback',str);
            end
        end
    end
 
    % label
    %----------------------------------------------------------------------
    uicontrol(gcf,...
        'Units','Normalized',...
        'Position',[x0 - dx y0 2*sx sy],...
        'Style','text',...
        'String',sprintf('pE.%s',F{f}),...
        'ForegroundColor',[1 0 0],...
        'HorizontalAlignment','left',...
        'FontSize',12,...
        'Tag','tmp');
    y0 = y - dy - dy/2;
    
end
 
% connectivity matrices
%==========================================================================
pE      = handles.DCM.M.pE;
try, SA = pE.SA; catch SA = []; end
try, GE = pE.GE; catch GE = []; end
try, GI = pE.GI; catch GI = []; end

% SA matrix
%--------------------------------------------------------------------------
x0    = 0.7;
y0    = 0.3;
for i = 1:size(SA,1);
    for j   = 1:size(SA,2)
        x   = x0 + j*dx;
        y   = y0 - (i - 1)*dy;
        str = sprintf('handles.DCM.M.pE.SA(%i,%i)',i,j);
        str = ['handles=guidata(gcbo);' str '=str2num(get(gcbo,''String''));guidata(gcbo,handles);'];
        uicontrol(gcf,...
            'Units','Normalized',...
            'Position',[x y sx sy],...
            'Style','edit',...
            'String',sprintf('%-4.2f',SA(i,j)),...
            'enable','on',...
            'Tag','tmp',...
            'Callback',str);
    end
end
 
% GE matrix
%--------------------------------------------------------------------------
x0    = 0.7;
y0    = 0.2;
for i = 1:size(GE,1);
    for j   = 1:size(GE,2)
        x   = x0 + j*dx;
        y   = y0 - (i - 1)*dy;
        str = sprintf('handles.DCM.M.pE.GE(%i,%i)',i,j);
        str = ['handles=guidata(gcbo);' str '=str2num(get(gcbo,''String''));guidata(gcbo,handles);'];
        uicontrol(gcf,...
            'Units','Normalized',...
            'Position',[x y sx sy],...
            'Style','edit',...
            'String',sprintf('%-4.2f',GE(i,j)),...
            'enable','on',...
            'Tag','tmp',...
            'Callback',str);
    end
end
 
% GI matrix
%--------------------------------------------------------------------------
x0    = 0.7;
y0    = 0.1;
for i = 1:size(GI,1);
    for j   = 1:size(GI,2)
        x   = x0 + j*dx;
        y   = y0 - (i - 1)*dy;
        str = sprintf('handles.DCM.M.pE.GI(%i,%i)',i,j);
        str = ['handles=guidata(gcbo);' str '=str2num(get(gcbo,''String''));guidata(gcbo,handles);'];
        uicontrol(gcf,...
            'Units','Normalized',...
            'Position',[x y sx sy],...
            'Style','edit',...
            'String',sprintf('%-4.2f',GI(i,j)),...
            'enable','on',...
            'Tag','tmp',...
            'Callback',str);
    end
end
 
% onset time
%==========================================================================
try
    ons  = handles.DCM.M.ons;
catch
    ons  = 32;
    handles.DCM.M.ons = ons;
end 
str = 'handles.DCM.M.ons';
str = ['handles=guidata(gcbo);' str '=str2num(get(gcbo,''String''));guidata(gcbo,handles);'];
uicontrol(gcf,...
    'Units','Normalized',...
    'Position',[0.848 0.341 0.071 0.022],...
    'Style','edit',...
    'String',sprintf('%-4.0f',ons),...
    'enable','on',...
    'Tag','tmp',...
    'Callback',str);
 
guidata(hObject,handles);
 
 
% --- Executes on button press in reset.
%--------------------------------------------------------------------------
function handles = reset_Callback(hObject, eventdata, handles)

try
    DCM          = handles.DCM;
catch
    load_Callback(hObject, eventdata, handles)
    return
end
model            = DCM.options.model;
[pE,pC]          = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);
handles.DCM.M.pE = pE;
handles.DCM.M.pC = pC;
guidata(hObject,handles);
unpack_Callback(hObject, eventdata, handles);
 
% --- Executes on button press in conditional.
%--------------------------------------------------------------------------
function conditional_Callback(hObject, eventdata, handles)
 
% conditional moments on parameters
%--------------------------------------------------------------------------
try
    handles.DCM.M.pE = handles.DCM.Ep;
    unpack_Callback(hObject, eventdata, handles);
catch
    warndlg('please invert and load a DCM to obtain conditional estimates')
end

% --- Executes on selection change in model.
function model_Callback(hObject, eventdata, handles)

% model type
%--------------------------------------------------------------------------
model = get(handles.model,       'String');
model = model{get(handles.model, 'Value')};
handles.DCM.options.model = model;
try
    reset_Callback(hObject, eventdata, handles);
catch
    load_Callback(hObject, eventdata, handles)
    return
end



