function varargout = spm_api_erp(varargin)
% SPM_API_ERP Application M-file for spm_api_erp.fig
%    FIG = SPM_API_ERP launch spm_api_erp GUI.
%    SPM_API_ERP('callback_name', ...) invoke the named callback.
%__________________________________________________________________________

% Last Modified by GUIDE v2.5 06-Sep-2006 15:18:43

if nargin == 0  % LAUNCH GUI

    fig     = openfig(mfilename,'reuse');
    Fgraph  = spm_figure('GetWin','Graphics');
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);

    DCM.name = 'ERP';
    DCM.Y.xy = {};

    handles.DCM = DCM;
    handles.X   = [];
    handles.Fgraph = Fgraph;
    guidata(fig, handles);

    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end

end


% DCM files and directories
%==========================================================================

% -------------------------------------------------------------------------
function varargout = swd_Callback(h, eventdata, handles, varargin)
try
    cd(get(handles.swd,'String'));
    set(handles.swd,'String',pwd)
end
handles.DCM.swd = pwd;
guidata(h,handles);

% -------------------------------------------------------------------------
function varargout = cd_Callback(h, eventdata, handles, varargin)
try
    p = uigetdir(pwd, 'select SWD');

    if p ~= 0
        cd(p);
    else, return, end
    set(handles.swd, 'String', pwd)
end
handles.DCM.swd = pwd;
guidata(h,handles);

% -------------------------------------------------------------------------
function varargout = name_Callback(h, eventdata, handles, varargin)
handles.DCM.name   = get(handles.name,'String');
set(handles.save,'enable','on')
guidata(h,handles);


%-Load and save
%==========================================================================

% -------------------------------------------------------------------------
function varargout = load_Callback(h, eventdata, handles, varargin)

[f,p] = uigetfile('DCM*.mat','please select DCM');
if p ~= 0
    DCM         = load(fullfile(p,f), '-mat');
    DCM         = DCM.DCM;
    handles.DCM = DCM;
    guidata(h, handles)
else, return, end

% enter options from saved options and execute data_ok and spatial_ok
try, set(handles.Y1, 'String', num2str(DCM.options.trials));      end
try, set(handles.T1, 'String', num2str(DCM.options.Tdcm(1)));     end
try, set(handles.T2, 'String', num2str(DCM.options.Tdcm(2)));     end

try, set(handles.Spatial_type,'Value', DCM.options.Spatial_type); end
try, set(handles.Nmodes,      'Value', DCM.options.Nmodes);       end
try, set(handles.h,           'Value', DCM.options.h + 1);        end
try, set(handles.D ,          'Value', DCM.options.D);            end
try, set(handles.onset,       'String',num2str(DCM.M.onset));     end
try, set(handles.design,      'String',num2str(DCM.U.X'));        end
try, set(handles.Uname,       'String',DCM.U.name);               end
try, set(handles.Vlocation,   'String',num2str(DCM.M.dipfit.L.Vpos(1,1))); 
end

handles.DCM = DCM;
guidata(h, handles);
% data_ok
%--------------------------------------------------------------------------
if isfield(DCM.options,'data_ok')
    try
        handles = data_ok_Callback(h, eventdata, handles);
    end
end

try, set(handles.name, 'String', DCM.name);           end
try, set(handles.Sname,'String', strvcat(DCM.Sname)); end

if DCM.options.Spatial_type ~= 3
    try, set(handles.Slocation, 'String', num2str(DCM.M.dipfit.L.pos')); end
end
guidata(h, handles);

% spatial_ok
%--------------------------------------------------------------------------
if isfield(DCM.options, 'spatial_ok')
    try
        handles = spatial_ok_Callback(h, eventdata, handles);
    end
    guidata(h, handles);

    try
        spm_api_erp('connections_Callback', gcbo, [], handles);
        set(handles.estimate,   'Enable', 'on');
        set(handles.initialise, 'Enable', 'on');
    end
end

try
    handles.DCM.F;
    set(handles.results,    'Enable','on');
    set(handles.estimate,   'String','Re-estimate','Foregroundcolor', [0 0 0]);
    set(handles.initialise, 'Enable','on');
catch
    set(handles.results,    'Enable','off');
    set(handles.estimate,   'String','Estimate','Foregroundcolor', [0 0 0]);
    set(handles.initialise, 'Enable','on');
end

% -------------------------------------------------------------------------
function handles = save_Callback(h, eventdata, handles, varargin)

DCM   = handles.DCM;
[f,p] = uiputfile(['DCM' DCM.name '.mat'],'save DCM');

if p ~= 0
    
    % store all the options for loading again
    % ---------------------------------------------------------------------
    DCM.options.trials       = str2num(get(handles.Y1,   'String'));
    DCM.options.Tdcm(1)      = str2num(get(handles.T1,   'String'));
    DCM.options.Tdcm(2)      = str2num(get(handles.T2,   'String'));
    DCM.options.Spatial_type = get(handles.Spatial_type, 'Value');
    DCM.options.Nmodes       = get(handles.Nmodes,       'Value');
    DCM.options.h            = get(handles.h,            'Value') - 1;
    DCM.options.D            = get(handles.D,            'Value');
    
    DCM.M.onset        = str2num(get(handles.onset,      'String'));
    DCM.M.Spatial_type = DCM.options.Spatial_type;
    
    handles.DCM = DCM;

    save(fullfile(p,f),     'DCM')
    set(handles.estimate,   'Enable', 'on')
    set(handles.initialise, 'Enable', 'on');
    set(handles.results,    'Enable', 'off')

    try
        handles.DCM.F;
        set(handles.estimate,'String','re-estimate')
    end
end
guidata(h,handles);


% Data selection and design
%==========================================================================

%-Get trials to model
%--------------------------------------------------------------------------
function Y1_Callback(hObject, eventdata, handles)
% hObject    handle to Y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles.DCM.options.trials = str2num(get(handles.Y1,'String'));
    m  = length(handles.DCM.options.trials);
catch
    return
end
try
    handles.DCM.U.X;
catch
    Xdefault(hObject,handles,m);
    warndlg({'Your design matrix has been initialised','Please check it'})
end
if m ~= size(handles.DCM.U.X,1)
    Xdefault(hObject,handles,m);
    warndlg({'Your design matrix has been re-set','Please check it'})
end
set(handles.design,'enable','on')
set(handles.Uname, 'enable','on')
set(handles.Y,'enable','on')
guidata(hObject,handles);

%--------------------------------------------------------------------------
function Uname_Callback(h, eventdata, handles)
% hObject    handle to Uname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    m  = length(handles.DCM.options.trials);
catch
    warndlg('please select trials')
    Xdefault(h,handles,1);
    return
end
Uname = get(handles.Uname,'String');
try 
    n  = size(handles.DCM.U.X,2);
    handles.DCM.U.name = Uname(1:n);
catch
    handles.DCM.U.name = Uname;
end
set(handles.Uname, 'String',handles.DCM.U.name);
guidata(h,handles);

%--------------------------------------------------------------------------
function design_Callback(h, eventdata, handles)
% hObject    handle to design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    m  = length(handles.DCM.options.trials);
catch
    Xdefault(h,handles,1);
    warndlg('please select trials')
    return
end
try
    X               = str2num(get(handles.design,'String'))';
    handles.DCM.U.X = X(1:m,:);
    set(handles.design, 'String',num2str(handles.DCM.U.X'));
catch
    Xdefault(h,handles,m);
end

guidata(h,handles);

%-Get data
% -------------------------------------------------------------------------
function varargout = Y_Callback(h, eventdata, handles, varargin)
DCM   = handles.DCM;
[f,p] = uigetfile({'*.mat'}, 'please select ERP/ERF.mat file (SPM format)');
if p ~=0
    f           = fullfile(p,f);
    DCM.M.Dfile = f;
    D           = spm_eeg_ldata(f);
else, return, end

% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
DCM.M.Ichannels = setdiff(D.channels.eeg, D.channels.Bad);
DCM.Y.Time      = [-D.events.start:D.events.stop]*1000/D.Radc; % ms
DCM.Y.dt        = 1000/D.Radc;

% if MEG, store grad struct in D.channels
%--------------------------------------------------------------------------
try
    handles.grad = D.channels.grad;
end

try
    trial = DCM.options.trials
    for i = 1:length(trial);
        DCM.Y.xy{i} = D.data(DCM.M.Ichannels,:,trial(i))';
    end
catch
    warndlg('please specify appropriate trials');
    return
end

% show data
%--------------------------------------------------------------------------
try
    spm_dcm_erp_results(DCM, 'Response');
end

% check design
%--------------------------------------------------------------------------
try
    m  = length(DCM.Y.xy);
    X  = DCM.U.X;

    % default design
    %----------------------------------------------------------------------
    if size(X,1) ~= m

        Xdefault(h,handles,m);
        warndlg({'Your design matrix has been re-set','Please check it'})

    end
end
set(handles.design,'enable','on')
handles.DCM = DCM;
guidata(h,handles);

% Spatial model specification
%==========================================================================

%--------------------------------------------------------------------------
function Spatial_type_Callback(hObject, eventdata, handles)
% hObject    handle to Spatial_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.Spatial_type,'Value') == 2
    % if MEG ECD
    set(handles.sensorfile,   'Enable', 'off');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation,    'Enable', 'on');
else
    % if EEG ECD
    set(handles.sensorfile,   'Enable', 'on');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation,    'Enable', 'on');
end
guidata(hObject,handles);


%---------------------------------------------------------------------
function sensorfile_Callback(h, eventdata, handles)
% hObject    handle to sensorfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Spatial_type = get(handles.Spatial_type, 'Value');
DCM          = handles.DCM;

% EEG
%---------------------------------------------------------------------
if Spatial_type == 1
    [f,p] = uigetfile({'*.pol';'*.mat'}, 'please select sensor location file');
    if p == 0, return, end
    f                       = fullfile(p,f);
    DCM.M.dipfit.sensorfile = f;
    handles.DCM             = DCM;
end
guidata(h, handles);


%---------------------------------------------------------------------
function handles = spatial_ok_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% spatial model
%------------------------
DCM = handles.DCM;
tmp = deblank(get(handles.Sname, 'String'));
if size(tmp, 1) == 0
    warndlg('Please specify the source names'); return
end
k = 1;
for i = 1:size(tmp, 1)
    if ~isempty(deblank(tmp(i,:)))
        Sname{k} = deblank(tmp(i,:));
        k = k + 1;
    end
end
Nareas    = length(Sname);
Nmodes    = get(handles.Nmodes, 'Value');
DCM.Sname = Sname;

% switch for spatial forward model
Spatial_type             = get(handles.Spatial_type, 'Value');
DCM.options.Spatial_type = Spatial_type;
DCM.M.Spatial_type       = Spatial_type;

% read location coordinates
Slocations = zeros(Nareas, 3);
tmp        = get(handles.Slocation, 'String');
if ~isempty(tmp) & size(tmp,1) == Nareas
    for i = 1:Nareas
        tmp2 = str2num(tmp(i, :));
        if length(tmp2) ~= 3
            errordlg(sprintf('coordinates of area %d not valid',i)); return;
        else
            Slocations(i, :) = tmp2;
        end
    end
    if size(Slocations, 1) ~= Nareas
        errordlg('Number of source names and locations must correspond'); return;
    end
else
    errordlg(sprintf('Please specify %d source locations.', Nareas)); return
end
if Spatial_type == 1
    try
        f = DCM.M.dipfit.sensorfile;
    catch
        errordlg('Please specify sensor location file'); return;
    end
elseif Spatial_type == 2
    try
    DCM.M.grad = handles.grad;
    catch
        warndlg({'Grad. is not defined',...
            'Please ensure this is MEG data'})
        return
    end
end

% prepare forward model
DCM                 = spm_dcm_erp_prepareSpatial(DCM);
DCM.M.dipfit.L.pos  = Slocations';
DCM.M.dipfit.L.mom  = zeros(3,Nareas);
DCM.M.dipfit.L.Vpos = str2num(get(handles.Vlocation,'String'))*ones(3, Nareas);
DCM.M.dipfit.L.Vmom = [8*ones(3, Nareas)];

DCM.M.Lpos = NaN*DCM.M.dipfit.L.pos; %
DCM.M.Lmom = NaN*DCM.M.dipfit.L.mom; %

DCM.options.spatial_ok = 1;

handles.DCM = DCM;
set(handles.Spatial_type,     'Enable', 'off');
set(handles.sensorfile,       'Enable', 'off');
set(handles.spatial_ok,       'Enable', 'off');
set(handles.Sname,            'Enable', 'off');
set(handles.Slocation,        'Enable', 'off');
set(handles.spatial_back,     'Enable', 'off');
set(handles.Vlocation,        'Enable', 'off');
set(handles.connections,       'Enable', 'on');
set(handles.connectivity_back, 'Enable', 'on');

% check whether user has already specified connections earlier
try
    if Nareas ~= handles.Nareas_old
        reset_Callback(hObject, eventdata, handles);
        handles.DCM = rmfield(handles.DCM,{'A','B','C'});
    end
end
guidata(hObject, handles);


% --- Executes on button press in data_ok.
%--------------------------------------------------------------------------
function handles = data_ok_Callback(hObject, eventdata, handles)
% hObject    handle to data_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DCM        = handles.DCM;
Nresponses = length(DCM.Y.xy);
if Nresponses == 0
    errordlg('Please choose some data'); return;
end
T1     = str2num(get(handles.T1, 'String'));
T2     = str2num(get(handles.T2, 'String'));
if T2 <= T1
    warnldg('Second must be greater than first peri-stimulus time')
    return
end
Nchannels = size(DCM.Y.xy{1}, 2);
Nmodes    = get(handles.Nmodes, 'Value');

% show data
%----------------------------------------------------------------------
spm_dcm_erp_results(DCM, 'Response');
DCM.options.data_ok = 1;
handles.DCM = DCM;

% enable next stage, disable data specification
set(handles.Y,             'Enable', 'off');
set(handles.Y1,            'Enable', 'off');
set(handles.T1,            'Enable', 'off');
set(handles.T2,            'Enable', 'off');
set(handles.Nmodes,        'Enable', 'off');
set(handles.h,             'Enable', 'off');
set(handles.D,             'Enable', 'off');
set(handles.design,        'Enable', 'off');
set(handles.data_ok,       'Enable', 'off');
set(handles.Spatial_type,  'Enable', 'on');
set(handles.plot_dipoles,  'Enable', 'on');
% ECD EEG or MEG
if get(handles.Spatial_type,  'Value') == 1
    set(handles.sensorfile,   'Enable', 'on');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation,    'Enable', 'on');
    set(handles.Vlocation,    'Enable', 'on');
elseif get(handles.Spatial_type, 'Value') == 2
    set(handles.sensorfile,   'Enable', 'off');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation,    'Enable', 'on');
    set(handles.Vlocation,    'Enable', 'on');
end

set(handles.spatial_ok,    'Enable', 'on');
set(handles.Sname,         'Enable', 'on');
set(handles.Slocation,     'Enable', 'on');
set(handles.spatial_back,  'Enable', 'on');
set(handles.Vlocation,     'Enable', 'on');
guidata(hObject, handles);


% --- Executes on button press in spatial_back.
%----------------------------------------------------------------------
function spatial_back_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Y,            'Enable', 'on');
set(handles.Y1,           'Enable', 'on');
set(handles.T1,           'Enable', 'on');
set(handles.T2,           'Enable', 'on');
set(handles.Nmodes,       'Enable', 'on');
set(handles.h,            'Enable', 'on');
set(handles.D,            'Enable', 'on');
set(handles.design,       'Enable', 'on');
set(handles.data_ok,      'Enable', 'on');
set(handles.Spatial_type, 'Enable', 'off');
set(handles.sensorfile,   'Enable', 'off');
set(handles.spatial_ok,   'Enable', 'off');
set(handles.Sname,        'Enable', 'off');
set(handles.Slocation,    'Enable', 'off');
set(handles.spatial_back, 'Enable', 'off');
set(handles.plot_dipoles, 'Enable', 'off');
set(handles.Vlocation,    'Enable', 'off');
guidata(hObject, handles);



% --- Executes on button press in plot_dipoles.
%--------------------------------------------------------------------------
function plot_dipoles_Callback(hObject, eventdata, handles)
% hObject    handle to plot_dipoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% read location coordinates
tmp = get(handles.Slocation, 'String');
Slocations = [];
if ~isempty(tmp) 
    for i = 1:size(tmp, 1)
        tmp2 = str2num(tmp(i, :))';
        if length(tmp2) == 3
            Slocations = [Slocations tmp2];
        end
    end
end

Nlocations   = size(Slocations, 2);
sdip.n_seeds = 1;
sdip.n_dip   = Nlocations;
sdip.Mtb     = 1;
sdip.j{1}    = zeros(3*Nlocations, 1);
sdip.loc{1}  = Slocations;
spm_eeg_inv_ecd_DrawDip('Init', sdip)


%-Connectivity
%==========================================================================
function varargout = connections_Callback(h, eventdata, handles, varargin)
n = length(handles.DCM.Sname);    % number of sources
m = size(handles.DCM.U.X,2);      % number of inputs

% onset
%--------------------------------------------------------------------------
handles.DCM.M.onset = str2num(get(handles.onset, 'String'));

% reset connection buttons
%--------------------------------------------------------------------------
reset_Callback(h, eventdata, handles, varargin)

% connection buttons
%--------------------------------------------------------------------------
set(handles.connections,'Units','Normalized')
p  = get(handles.connections,'Position');
x0 = p(1);
y0 = p(2);
sx = 1/32;
sy = 1/64;
for i = 1:n
    for j = 1:n
        for k = 1:3
            x          = x0 + (j - 1)*sx + (k - 1)*(n + 1)*sx;
            y          = y0 - (i + 4)*sy;
            str        = sprintf('data.DCM.A{%i}(%i,%i)',k,i,j);
            str        = ['data=guidata(gcbo);' str '=get(gcbo,''Value'');guidata(gcbo,data)'];
            A{k}(i,j)  = uicontrol(handles.SPM,...
                'Units','Normalized',...
                'Position',[x y sx sy],...
                'Style','radiobutton',...
                'Tag','',...
                'Callback',str);
            if i == j
                set(A{k}(i,j),'Enable','off')
            end
            try
                set(A{k}(i,j),'Value',handles.DCM.A{k}(i,j));
            catch
                handles.DCM.A{k}(i,j) = get(A{k}(i,j),'Value');
            end
        end
        for k = 1:m
            x          = x0 + (j - 1)*sx + (k - 1)*(n + 1)*sx;
            y          = y0 - (i + 4)*sy - (n + 1)*sy;
            str        = sprintf('data.DCM.B{%i}(%i,%i)',k,i,j); 
            str        = ['data=guidata(gcbo);' str '=get(gcbo,''Value'');guidata(gcbo,data)'];
            B{k}(i,j)  = uicontrol(handles.SPM,...
                'Units','Normalized',...
                'Position',[x y sx sy],...
                'Style','radiobutton',...
                'Tag','',...
                'Callback',str);
            if i == j
                set(B{k}(i,j),'Enable','off')
            end
            try
                set(B{k}(i,j),'Value',handles.DCM.B{k}(i,j));
            catch
                handles.DCM.B{k}(i,j) = get(B{k}(i,j),'Value');
            end
        end
    end
end
for i = 1:n
    x          = x0 + (4 - 1)*(n + 1)*sx;
    y          = y0 - (i + 4)*sy;
    str        = sprintf('data.DCM.C(%i)',i);
    str        = ['data=guidata(gcbo);' str '=get(gcbo,''Value'');guidata(gcbo,data)'];
    C(i)       = uicontrol(handles.SPM,...
        'Units','Normalized',...
        'Position',[x y sx sy],...
        'Style','radiobutton',...
        'Tag','',...
        'Callback',str);
    try
        set(C(i),'Value',handles.DCM.C(i));
    catch
        handles.DCM.C(i,1) = get(C(i),'Value');
    end
end

% string labels
%--------------------------------------------------------------------------
constr         = {'A forward' 'A backward' 'A lateral' 'C input'};
nsx            = (n + 1)*sx;
nsy            = 2*sy;
for k = 1:4
    x          = x0 + (k - 1)*nsx;
    y          = y0 - 3*sy;
    str        = constr{k};
    S(k)       = uicontrol(handles.SPM,...
        'Units','Normalized',...
        'Position',[x y nsx nsy],...
        'HorizontalAlignment','left',...
        'Style','text',...
        'String',str,...
        'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
constr         = handles.DCM.U.name;
for k = 1:m
    x          = x0 + (k - 1)*nsx;
    y          = y0 - 6*sy - 2*(n + 1)*sy;
    str        = ['B ' constr{k}];
    S(4 + k)   = uicontrol(handles.SPM,...
        'Units','Normalized',...
        'Position',[x y nsx nsy],...
        'HorizontalAlignment','left',...
        'Style','text',...
        'String',str,...
        'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
handles.S = S;
handles.A = A;
handles.B = B;
handles.C = C;
set(handles.connections,'Enable','off')
set(handles.reset,      'Enable','on' )
set(handles.estimate,   'Enable', 'on');
set(handles.initialise, 'Enable', 'on');
guidata(h,handles)

% remove existing buttons
%--------------------------------------------------------------------------
function varargout = reset_Callback(h, eventdata, handles, varargin)
try
    for i = 1:length(handles.A)
        for j = 1:length(handles.A{i})
            for k = 1:length(handles.A{i})
                delete(handles.A{i}(j,k));
            end
        end
    end
    for i = 1:length(handles.B)
        for j = 1:length(handles.B{i})
            for k = 1:length(handles.B{i})
                delete(handles.B{i}(j,k));
            end
        end
    end
    for i = 1:length(handles.C)
        delete(handles.C(i));
    end
    for i = 1:length(handles.S)
        delete(handles.S(i));
    end
    handles.DCM = rmfield(handles.DCM,{'A','B','C'});
    handles     = rmfield(handles,{'A','B','C','S'});
    guidata(h,handles)
    
end
set(handles.connections,'Enable','on')

% --- Executes on button press in connectivity_back.
%--------------------------------------------------------------------------
function connectivity_back_Callback(hObject, eventdata, handles)
% hObject    handle to connectivity_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Spatial_type,      'Enable', 'on');

if get(handles.Spatial_type,'Value') == 1
    set(handles.sensorfile,    'Enable', 'on');
else
    set(handles.sensorfile,    'Enable', 'off');
end

set(handles.spatial_ok,        'Enable', 'on');
set(handles.Sname,             'Enable', 'on');
set(handles.Slocation,         'Enable', 'on');
set(handles.spatial_back,      'Enable', 'on');
set(handles.Vlocation,         'Enable', 'on');

set(handles.connections,       'Enable', 'off');
set(handles.connectivity_back, 'Enable', 'off');
set(handles.reset,             'Enable', 'off');

% connection buttons
try
    for i = 1:length(handles.A)
        for j = 1:length(handles.A{i})
            for k = 1:length(handles.A{i})
                set(handles.A{i}(j,k), 'Enable', 'off');
            end
        end
    end
    for i = 1:length(handles.B)
        for j = 1:length(handles.B{i})
            for k = 1:length(handles.B{i})
                set(handles.B{i}(j,k), 'Enable', 'off');
            end
        end
    end
    for i = 1:length(handles.C)
        set(handles.C(i), 'Enable', 'off');
    end
end
% save old nr of areas
handles.Nareas_old = size(handles.A{1}, 1);
guidata(hObject,handles)


%-Estimate, initalise and review
%==========================================================================

% -------------------------------------------------------------------------
function varargout = estimate_Callback(h, eventdata, handles, varargin)
set(handles.estimate,'String','Estimating','Foregroundcolor',[1 0 0])

handles.DCM = spm_dcm_erp(handles.DCM);
save_Callback(h, eventdata, handles, varargin)
set(handles.results,    'Enable','on' )
set(handles.save,       'Enable','off')
set(handles.Y,          'Enable','off')
set(handles.name,       'Enable','off')
set(handles.swd,        'Enable','off')
set(handles.reset,      'Enable','off')

% -------------------------------------------------------------------------
function varargout = results_Callback(h, eventdata, handles, varargin)
Action  = get(handles.results, 'String');
Action  = Action{get(handles.results, 'Value')};
spm_dcm_erp_results(handles.DCM, Action);

% -------------------------------------------------------------------------
function initialise_Callback(h, eventdata, handles)
% hObject    handle to initialise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[f,p]           = uigetfile('DCM*.mat','please select estimated DCM');
DCM             = load(fullfile(p,f), '-mat');
handles.DCM.M.P = DCM.DCM.Ep;
guidata(h, handles);


% Auxillary functions
%==========================================================================
function [X,name] = Xdefault(hObject,handles,m);
% default design matix
% m - number of trials
X       = eye(m);
X(:,1)  = [];
name    = {};
for i = 1:size(X,2)
   name{i,1} = sprintf('effect %i',i);
end
handles.DCM.U.X    = X;
handles.DCM.U.name = name;
set(handles.design,'String',num2str(handles.DCM.U.X'));
set(handles.Uname, 'String',handles.DCM.U.name);
guidata(hObject,handles);
return

