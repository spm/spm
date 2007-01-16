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
    handles  = guihandles(fig);
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


% Callbacks
%==========================================================================

%-DCM files and directories: Load and save
%==========================================================================

% load
% -------------------------------------------------------------------------
function varargout = load_Callback(hObject, eventdata, handles, varargin)
try
    [f,p]       = uigetfile('*.mat','please select DCM file');
    name        = fullfile(p,f);
    DCM         = load(name,'-mat');
    DCM         = DCM.DCM;
    DCM.name    = name;
    handles.DCM = DCM;
    guidata(hObject,handles);
catch
    return;
end

% Filename
%--------------------------------------------------------------------------
try, set(handles.name,'String',f); end

% enter options from saved options and execute data_ok and spatial_ok
%--------------------------------------------------------------------------
try, set(handles.Y1, 'String', num2str(DCM.options.trials));      end
try, set(handles.T1, 'String', num2str(DCM.options.Tdcm(1)));     end
try, set(handles.T2, 'String', num2str(DCM.options.Tdcm(2)));     end

try, set(handles.Spatial_type,'Value', DCM.options.Spatial_type); end
try, set(handles.Nmodes,      'Value', DCM.options.Nmodes);       end
try, set(handles.h,           'Value', DCM.options.h + 1);        end
try, set(handles.D ,          'Value', DCM.options.D);            end
try, set(handles.onset,       'String',num2str(DCM.M.onset));     end
try, set(handles.design,      'String',num2str(DCM.xU.X'));       end
try, set(handles.Uname,       'String',DCM.xU.name);              end

set(handles.data_ok,'enable','on')

handles.DCM = DCM;
guidata(hObject, handles);

% data_ok
%--------------------------------------------------------------------------
if isfield(DCM.options,'data_ok')
    try
        handles = data_ok_Callback(hObject, eventdata, handles);
    end
end

try, set(handles.Sname,'String', strvcat(DCM.Sname)); end

if DCM.options.Spatial_type ~= 3
    try, set(handles.Slocation, 'String', num2str(DCM.M.dipfit.L.pos')); end
end
guidata(hObject, handles);

% spatial_ok
%--------------------------------------------------------------------------
if isfield(DCM.options, 'spatial_ok')
    try
        handles = spatial_ok_Callback(hObject, eventdata, handles);
    end
    guidata(hObject, handles);

    try
        connections_Callback(hObject, eventdata, handles);
        set(handles.estimate,   'Enable', 'on');
        set(handles.initialise, 'Enable', 'on');
    end
end

try
    handles.DCM.F;
    set(handles.results,    'Enable','on');
catch
    set(handles.results,    'Enable','off');
end

% -------------------------------------------------------------------------
function handles = save_Callback(hObject, eventdata, handles, varargin)

DCM   = handles.DCM;
[f,p] = uiputfile(['DCM*.mat'],'DCM file to save');

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
    
    DCM.M.onset = str2num(get(handles.onset,      'String'));
    DCM.name    = fullfile(p,f);
    set(handles.name,'String',f);
    
    handles.DCM = DCM;
    save(DCM.name,          'DCM')
    set(handles.estimate,   'Enable', 'on')
    set(handles.initialise, 'Enable', 'on');

end
guidata(hObject,handles);


% Data selection and design
%==========================================================================

%-Get trials and data
%--------------------------------------------------------------------------
function Y1_Callback(hObject, eventdata, handles)

try
    handles.DCM.options.trials = str2num(get(handles.Y1,'String'));
    m  = length(handles.DCM.options.trials);
catch
    m  = 1;
    set(handles.Y1,'String','1')
end
handles = Xdefault(hObject,handles,m);

DCM = handles.DCM;

%-Get new trial date from file
%--------------------------------------------------------------------------
[f,p] = uigetfile({'*.mat'}, 'please select ERP/ERF.mat file (SPM format)');
DCM.xY.Dfile = fullfile(p,f);
DCM          = spm_dcm_erp_data(DCM);

% Display data
%--------------------------------------------------------------------------
try
    spm_dcm_erp_results(DCM, 'Response');
end

set(handles.design,'enable','on')
set(handles.Uname, 'enable','on')
set(handles.Y,     'enable','on')
handles.DCM = DCM;
guidata(hObject,handles);
warndlg({'Your design matrix has been re-set'})

%-Display data
%--------------------------------------------------------------------------
function Y_Callback(hObject, eventdata, handles)
try
    handles.DCM = spm_dcm_erp_data(handles.DCM);
    spm_dcm_erp_results(handles.DCM,'Response');
    guidata(hObject,handles);
catch
    warndlg('please specify trials to load');
end


%--------------------------------------------------------------------------
function Uname_Callback(hObject, eventdata, handles)
try
    m = length(handles.DCM.options.trials);
catch
    warndlg('please select trials')
    handles = Xdefault(hObject,handles,1);
    return
end
Str   = get(handles.Uname,'String');
n     = size(handles.DCM.xU.X,2);
for i = 1:n
    try
        Uname{i} = Str{i};
    catch
        Uname{i} = sprintf('effect %i',i);
    end
end
set(handles.Uname,'string',Uname(:))
handles.DCM.xU.name = Uname;
guidata(hObject,handles);

%--------------------------------------------------------------------------
function design_Callback(hObject, eventdata, handles)
try
    m  = length(handles.DCM.options.trials);
catch
    handles = Xdefault(hObject,handles,1);
    warndlg('please select trials')
    return
end
try
    X  = str2num(get(handles.design,'String'))';
    handles.DCM.xU.X = X(1:m,:);
    set(handles.design, 'String',num2str(handles.DCM.xU.X'));
catch
    handles = Xdefault(hObject,handles,m);
end
n     = size(handles.DCM.xU.X,2);
Uname = {};
for i = 1:n
    try
        Uname{i} = handles.DCM.xU.name{i};
    catch
        Uname{i} = sprintf('effect %i',i);
    end
end
set(handles.Uname,'string',Uname(:))
handles.DCM.xU.name = Uname;
guidata(hObject,handles);


% Spatial model specification
%==========================================================================

%--------------------------------------------------------------------------
function Spatial_type_Callback(hObject, eventdata, handles)

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


%--------------------------------------------------------------------------
function sensorfile_Callback(hObject, eventdata, handles)
Spatial_type = get(handles.Spatial_type, 'Value');

if Spatial_type == 1
    [f,p] = uigetfile({'*.pol';'*.mat'}, 'please select sensor location file');
    handles.DCM.M.dipfit.sensorfile = fullfile(p,f);
end
guidata(hObject, handles);


%--------------------------------------------------------------------------
function handles = spatial_ok_Callback(hObject, eventdata, handles)

% spatial model
%--------------------------------------------------------------------------
DCM = handles.DCM;
tmp = deblank(get(handles.Sname, 'String'));
if size(tmp,1) == 0
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
%--------------------------------------------------------------------------
Spatial_type             = get(handles.Spatial_type, 'Value');
DCM.options.Spatial_type = Spatial_type;

% read location coordinates
%--------------------------------------------------------------------------
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
        DCM.M.dipfit.sensorfile;
    catch
        warndlg('Please specify sensor location file'); return;
    end
elseif Spatial_type == 2
    try
        handles.DCM.xY.grad;
    catch
        warndlg({'Grad. is not defined','Please ensure this is MEG data'})
        return
    end
end

% forward model (spatial)
%--------------------------------------------------------------------------
DCM = spm_dcm_erp_dipfit(DCM);

% set prior expectations about locations
%--------------------------------------------------------------------------
DCM.M.dipfit.L.pos  = Slocations';

DCM.options.spatial_ok = 1;

handles.DCM = DCM;
set(handles.Spatial_type,     'Enable', 'off');
set(handles.sensorfile,       'Enable', 'off');
set(handles.spatial_ok,       'Enable', 'off');
set(handles.Sname,            'Enable', 'off');
set(handles.Slocation,        'Enable', 'off');
set(handles.spatial_back,     'Enable', 'off');
set(handles.connections,      'Enable', 'on');
set(handles.connectivity_back,'Enable', 'on');

% [re]-set connections
%--------------------------------------------------------------------------
handles = connections_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

%-Executes on button press in data_ok.
%--------------------------------------------------------------------------
function handles = data_ok_Callback(hObject, eventdata, handles)
try
    handles.DCM = spm_dcm_erp_data(handles.DCM);
catch
    errordlg('Please choose some trials'); return;
end

% show data
%--------------------------------------------------------------------------
spm_dcm_erp_results(handles.DCM, 'Response');
handles.DCM.options.data_ok = 1;

% enable next stage, disable data specification
%--------------------------------------------------------------------------
set(handles.Y,             'Enable', 'off');
set(handles.Y1,            'Enable', 'off');
set(handles.T1,            'Enable', 'off');
set(handles.T2,            'Enable', 'off');
set(handles.Nmodes,        'Enable', 'off');
set(handles.h,             'Enable', 'off');
set(handles.D,             'Enable', 'off');
set(handles.design,        'Enable', 'off');
set(handles.Uname,         'Enable', 'off');
set(handles.data_ok,       'Enable', 'off');
set(handles.Spatial_type,  'Enable', 'on');
set(handles.plot_dipoles,  'Enable', 'on');
% ECD EEG or MEG
if get(handles.Spatial_type,  'Value') == 1
    set(handles.sensorfile,   'Enable', 'on');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation,    'Enable', 'on');
elseif get(handles.Spatial_type, 'Value') == 2
    set(handles.sensorfile,   'Enable', 'off');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation,    'Enable', 'on');
end

set(handles.spatial_ok,    'Enable', 'on');
set(handles.Sname,         'Enable', 'on');
set(handles.Slocation,     'Enable', 'on');
set(handles.spatial_back,  'Enable', 'on');
guidata(hObject, handles);


%-Executes on button press in spatial_back.
%----------------------------------------------------------------------
function spatial_back_Callback(hObject, eventdata, handles)

set(handles.Y,            'Enable', 'on');
set(handles.Y1,           'Enable', 'on');
set(handles.T1,           'Enable', 'on');
set(handles.T2,           'Enable', 'on');
set(handles.Nmodes,       'Enable', 'on');
set(handles.h,            'Enable', 'on');
set(handles.D,            'Enable', 'on');
set(handles.design,       'Enable', 'on');
set(handles.Uname,        'Enable', 'on');
set(handles.data_ok,      'Enable', 'on');
set(handles.Spatial_type, 'Enable', 'off');
set(handles.sensorfile,   'Enable', 'off');
set(handles.spatial_ok,   'Enable', 'off');
set(handles.Sname,        'Enable', 'off');
set(handles.Slocation,    'Enable', 'off');
set(handles.spatial_back, 'Enable', 'off');
set(handles.plot_dipoles, 'Enable', 'off');
guidata(hObject, handles);

%-Executes on button press in plot_dipoles.
%--------------------------------------------------------------------------
function plot_dipoles_Callback(hObject, eventdata, handles)

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
function handles = connections_Callback(hObject, eventdata, handles)

% delete previous connection buttons
%--------------------------------------------------------------------------
reset_Callback(hObject, eventdata, handles);

DCM = handles.DCM;
n   = length(DCM.Sname);    % number of sources
m   = size(DCM.xU.X,2);     % number of inputs


% no changes in coupling
%--------------------------------------------------------------------------
if ~m, B = {}; DCM.B = {}; end

% check DCM.A, DCM.B, ...
%--------------------------------------------------------------------------
try
    if length(DCM.A{1}) ~= n, DCM = rmfield(DCM,'A'); end
    if length(DCM.B)    ~= m, DCM = rmfield(DCM,'B'); end
    if length(DCM.C)    ~= n, DCM = rmfield(DCM,'C'); end
end

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
                set(A{k}(i,j),'Value',DCM.A{k}(i,j));
            catch
                DCM.A{k}(i,j) = get(A{k}(i,j),'Value');
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
%             % intrinsic modulation
%             if i == j
%                 set(B{k}(i,j),'Enable','on')
%             end
            try
                set(B{k}(i,j),'Value',DCM.B{k}(i,j));
            catch
                DCM.B{k}(i,j) = get(B{k}(i,j),'Value');
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
        set(C(i),'Value',DCM.C(i));
    catch
        DCM.C(i,1) = get(C(i),'Value');
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
constr         = DCM.xU.name;
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
handles.S   = S;
handles.A   = A;
handles.B   = B;
handles.C   = C;
handles.DCM = DCM;
set(handles.estimate,   'Enable','on');
set(handles.initialise, 'Enable','on');
guidata(hObject,handles)

% remove existing buttons
%--------------------------------------------------------------------------
function reset_Callback(hObject, eventdata, handles)
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
end
set(handles.connections,'Enable','on')

%-Executes on button press in connectivity_back.
%--------------------------------------------------------------------------
function connectivity_back_Callback(hObject, eventdata, handles)

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

set(handles.connections,       'Enable', 'off');
set(handles.connectivity_back, 'Enable', 'off');

% connection buttons
%--------------------------------------------------------------------------
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

%-Estimate, initalise and review
%==========================================================================

% -------------------------------------------------------------------------
function varargout = estimate_Callback(hObject, eventdata, handles, varargin)
set(handles.estimate,'String','Estimating','Foregroundcolor',[1 0 0])

% initialise if required
% -------------------------------------------------------------------------
try
    Ep  = handles.DCM.Ep;
    Str = questdlg('intialise with revious esimates');
    if strcmp(Str,'Yes')
        handles.DCM.M.P = Ep;
    elseif strcmp(Str,'No')
        handles.DCM.M.P = [];
    elseif strcmp(Str,'Cancel')
        return
    end
end

% get filename and save
% -------------------------------------------------------------------------
handles     = save_Callback(hObject, eventdata, handles);

% invert and save
% -------------------------------------------------------------------------
handles.DCM = spm_dcm_erp(handles.DCM);

assignin('base','DCM',handles.DCM)
guidata(hObject, handles);

set(handles.results,  'Enable','on' )
set(handles.save,     'Enable','on')
set(handles.estimate, 'String','Estimated','Foregroundcolor',[0 0 0])

% -------------------------------------------------------------------------
function varargout = results_Callback(hObject, eventdata, handles, varargin)
Action  = get(handles.results, 'String');
Action  = Action{get(handles.results, 'Value')};
spm_dcm_erp_results(handles.DCM, Action);

% -------------------------------------------------------------------------
function initialise_Callback(hObject, eventdata, handles)
[f,p]           = uigetfile('DCM*.mat','please select estimated DCM');
DCM             = load(fullfile(p,f), '-mat');
handles.DCM.M.P = DCM.DCM.Ep;
guidata(hObject, handles);


% Auxillary functions
%==========================================================================
function handles = Xdefault(hObject,handles,m);
% default design matix
% m - number of trials
X       = eye(m);
X(:,1)  = [];
name    = {};
for i = 1:size(X,2)
   name{i,1} = sprintf('effect %i',i);
end
handles.DCM.xU.X    = X;
handles.DCM.xU.name = name;
set(handles.design,'String',num2str(handles.DCM.xU.X'));
set(handles.Uname, 'String',handles.DCM.xU.name);
return

