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

% Last Modified by GUIDE v2.5 01-Oct-2008 19:46:25

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
DCM   = handles.DCM;
M     = DCM.M;
pst   = str2num(get(handles.pst,'String'));
Hz    = str2num(get(handles.Hz, 'String'));

% initalise states
%--------------------------------------------------------------------------
M     = DCM.M;
M.x   = spm_x_nmm(M.pE);

% Volterra Kernels
%==========================================================================

% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L1,L2] = spm_bireduce(M,M.pE);

% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
U.dt    = 4/1000;
M.ns    = pst/U.dt/1000;
t       = [1:M.ns]*U.dt*1000;
[K0,K1] = spm_kernels(M0,M1,L1,L2,M.ns,U.dt);

axes(handles.kernel)
plot(t,K1)
set(gca,'XLim',[min(t) max(t)])
ylabel('1st-order Volterra kernel')
xlabel('time (ms)')

% Transfer function (source space)
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


% --- Executes on button press in response.
%--------------------------------------------------------------------------
function response_Callback(hObject, eventdata, handles)
clear spm_erp_L spm_gen_erp
DCM   = handles.DCM;
pst   = str2num(get(handles.pst,'String'));
Hz    = str2num(get(handles.Hz, 'String'));

% initalise states
%--------------------------------------------------------------------------
M     = DCM.M;
M.x   = spm_x_nmm(M.pE);

U.dt  = 4/1000;
M.ns  = pst/U.dt/1000;
M.ons = 32;
t     = [1:M.ns]*U.dt*1000 - M.ons;

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

% get priors
%--------------------------------------------------------------------------
try
    handles.DCM.M.pE;
catch
    handles = reset_Callback(hObject, eventdata, handles);
end
pE = handles.DCM.M.pE;
try
    pE = rmfield(pE,'B');
end

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
                    uicontrol(gcf,...
                        'Units','Normalized',...
                        'Position',[x y sx sy],...
                        'Style','edit',...
                        'String',sprintf('%-4.0f',P{k}(i,j)),...
                        'enable','on',...
                        'Tag','tmp',...
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
                'String',sprintf('%-4.0f',P(i)),...
                'enable','on',...
                'Tag','tmp',...
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
                    'String',sprintf('%-4.0f',P(i,j)),...
                    'enable','on',...
                    'Tag','tmp',...
                    'Callback',str);
            end
        end
    end

    % label
    %----------------------------------------------------------------------
    uicontrol(gcf,...
        'Units','Normalized',...
        'Position',[x0 - 2*dx y0 2*sx sy],...
        'Style','text',...
        'String',sprintf('Prior: %s',F{f}),...
        'ForegroundColor',[1 0 0],...
        'HorizontalAlignment','left',...
        'FontSize',14,...
        'Tag','tmp');
    y0 = y - dy - dy/2;

end


% --- Executes on button press in reset.
%--------------------------------------------------------------------------
function handles = reset_Callback(hObject, eventdata, handles)

% prior moments on parameters
%--------------------------------------------------------------------------
DCM              = handles.DCM;
[pE,gE,pC,gC]    = spm_nmm_priors(DCM.A,DCM.B,DCM.C,DCM.M.dipfit);
handles.DCM.M.pE = pE;
guidata(hObject,handles);
unpack_Callback(hObject, eventdata, handles);



