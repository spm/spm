function varargout = spm_api_erp(varargin)
% SPM_API_ERP Application M-file for spm_api_erp.fig
%    FIG = SPM_API_ERP launch spm_api_erp GUI.
%    SPM_API_ERP('callback_name', ...) invoke the named callback.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_api_erp.m 1702 2008-05-21 13:55:11Z vladimir $

if nargin == 0 || nargin == 1  % LAUNCH GUI

    fig     = openfig(mfilename,'reuse');
    Fgraph  = spm_figure('GetWin','Graphics');
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

    % Generate a structure of handles to pass to callbacks, and store it.
    handles  = guihandles(fig);
    handles.Fgraph = Fgraph;
    guidata(fig, handles);
    
    if nargin == 1
        load_Callback(fig, [], handles, varargin{1})
    end

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
% --- Executes during object creation, after setting all properties.
function data_ok_CreateFcn(hObject, eventdata, handles)

set(hObject, 'enable', 'off');


%-DCM files and directories: Load and save
%==========================================================================

% --- Executes on button press in load.
% -------------------------------------------------------------------------
function varargout = load_Callback(hObject, eventdata, handles, varargin)
try
    DCM         = varargin{1};
    [p,f]       = fileparts(DCM.name);
catch
    [f,p]       = uigetfile('*.mat','please select DCM file'); 
    cd(p)
    name        = fullfile(p,f);
    DCM         = load(name,'-mat');
    DCM         = DCM.DCM;
    DCM.name    = name;
    handles.DCM = DCM;
    guidata(hObject,handles);
end


% Type of model
%--------------------------------------------------------------------------
% 'ERP'    - (linear NMM slow)
% 'SEP'    - (linear NMM fast)
% 'NMM'    - (bilinear NMM)
% 'MFM'    - (bilinear MFM)
% 'SSR'    - (linear NMM for steady-state responses)
% 'IND'    - (linear NMM for induced responses)

try
    model = DCM.options.model;
catch
    model = get(handles.ERP,'String');
    model = model{get(handles.ERP,'Value')};
end
switch model
    case{'ERP'},     set(handles.ERP,'Value',1);
    case{'SEP'},     set(handles.ERP,'Value',2);
    case{'NMM'},     set(handles.ERP,'Value',3);
    case{'MFM'},     set(handles.ERP,'Value',4);
    case{'SSR'},     set(handles.ERP,'Value',5);
    case{'IND'},     set(handles.ERP,'Value',6);
    otherwise
end
handles = ERP_Callback(hObject, eventdata, handles);

% Filename
%--------------------------------------------------------------------------
try, set(handles.name,'String',f);  end

% Source locations
%--------------------------------------------------------------------------
try, DCM.Lpos = DCM.M.dipfit.Lpos; end

% enter options from saved options and execute data_ok and spatial_ok
%--------------------------------------------------------------------------
try, set(handles.Y1, 'String', num2str(DCM.options.trials,'%7.0f')); end
try, set(handles.T1, 'String', num2str(DCM.options.Tdcm(1)));        end
try, set(handles.T2, 'String', num2str(DCM.options.Tdcm(2)));        end
try, set(handles.Hz1,'String', num2str(DCM.options.Fdcm(1)));        end
try, set(handles.Hz2,'String', num2str(DCM.options.Fdcm(2)));        end
try, set(handles.Rft,'String', num2str(DCM.options.Rft));            end
try, set(handles.Spatial_type,'Value', DCM.options.type);            end
try, set(handles.Nmodes,      'Value', DCM.options.Nmodes);          end
try, set(handles.h,           'Value', DCM.options.h);               end
try, set(handles.han,         'Value', DCM.options.han);             end
try, set(handles.D,           'Value', DCM.options.D);               end
try, set(handles.lock,        'Value', DCM.options.lock);            end
try, set(handles.design,      'String',num2str(DCM.xU.X','%7.2f'));  end
try, set(handles.Uname,       'String',DCM.xU.name);                 end
try, set(handles.Sname,       'String',DCM.Sname);                   end
try, set(handles.onset,       'String',num2str(DCM.options.onset));  end
try, set(handles.Slocation,   'String',num2str(DCM.Lpos','%4.0f'));  end


% catch for backwards compatibility
%--------------------------------------------------------------------------
if get(handles.Spatial_type,'Value') > length(get(handles.Spatial_type,'String'))
    set(handles.Spatial_type,'Value',1);
end

if get(handles.Spatial_type,'Value') == 2
    set(handles.render, 'Enable','on' )
    set(handles.Imaging,'Enable','on' )
else
    set(handles.render, 'Enable','off' )
    set(handles.Imaging,'Enable','off' )
end
set(handles.data_ok,'enable','on')

handles.DCM = DCM;
guidata(hObject, handles);

% data specification
%--------------------------------------------------------------------------
try
    handles = data_ok_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
catch
    return
end

% spatial model specification
%--------------------------------------------------------------------------
try
    handles = spatial_ok_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
catch
    return
end

% connections specification
%--------------------------------------------------------------------------
try
    connections_Callback(hObject, eventdata, handles);
    set(handles.estimate,   'Enable', 'on');
    set(handles.initialise, 'Enable', 'on');
    guidata(hObject, handles);
catch
    return
end

% estimation and results
%--------------------------------------------------------------------------
try
    handles.DCM.F;
    set(handles.results,    'Enable','on');
    if handles.DCM.options.type == 2
        set(handles.render, 'Enable','on');
        set(handles.Imaging,'Enable','on');

    end
catch
    set(handles.results,    'Enable','off');
end

guidata(hObject, handles);

% --- Executes on button press in save.
% -------------------------------------------------------------------------
function handles = save_Callback(hObject, eventdata, handles)

handles = reset_Callback(hObject, eventdata, handles);
try
   [p,f]  = fileparts(handles.DCM.name);
catch
    try
        [p,f] = fileparts(handles.DCM.xY.Dfile);
        f = ['DCM_' f];
    catch
        f = ['DCM' date];
    end
end
[f,p]     = uiputfile(['DCM*.mat'],'DCM file to save',f);

if p
    handles.DCM.name = fullfile(p,f);
    set(handles.name,'String',f);
    DCM              = handles.DCM;
    save(DCM.name,'DCM')
    set(handles.estimate,   'Enable', 'on')
    set(handles.initialise, 'Enable', 'on');
    cd(p)
end

% assign in base
%--------------------------------------------------------------------------
assignin('base','DCM',handles.DCM)
guidata(hObject,handles);

% store selections in DCM
% -------------------------------------------------------------------------
function handles = reset_Callback(hObject, eventdata, handles)

handles.DCM.options.trials   = str2num(get(handles.Y1,    'String'));
handles.DCM.options.Tdcm(1)  = str2num(get(handles.T1,    'String'));
handles.DCM.options.Tdcm(2)  = str2num(get(handles.T2,    'String'));
handles.DCM.options.Fdcm(1)  = str2num(get(handles.Hz1,   'String'));
handles.DCM.options.Fdcm(2)  = str2num(get(handles.Hz2,   'String'));
handles.DCM.options.Rft      = str2num(get(handles.Rft,   'String'));
handles.DCM.options.onset    = str2num(get(handles.onset, 'String'));
handles.DCM.options.Nmodes   = get(handles.Nmodes,        'Value');
handles.DCM.options.h        = get(handles.h,             'Value');
handles.DCM.options.han      = get(handles.han,           'Value');
handles.DCM.options.D        = get(handles.D,             'Value');
handles.DCM.options.type     = get(handles.Spatial_type,  'Value');
handles.DCM.options.lock     = get(handles.lock,          'Value');


% Check model type
%--------------------------------------------------------------------------
model = get(handles.ERP,     'String');
model = model{get(handles.ERP,'Value')};
handles.DCM.options.model  = model;

% Check data modality
%--------------------------------------------------------------------------
try
    switch handles.DCM.xY.modality

        case {'EEG','MEG'}
            set(handles.Slocation,     'Visable', 'on');
            set(handles.plot_dipoles,  'Enable',  'on');

        otherwise
            set(handles.Slocation,     'Visable', 'off');
            set(handles.plot_dipoles,  'Enable',  'off');
    end
end

guidata(hObject,handles);


% Data selection and design
%==========================================================================

% --- Executes on button press in Datafile.
%--------------------------------------------------------------------------
function Datafile_Callback(hObject, eventdata, handles)
Y1_Callback(hObject, eventdata, handles);


%-Get trials and data
%--------------------------------------------------------------------------
function Y1_Callback(hObject, eventdata, handles)
try
    trials = str2num(get(handles.Y1,'String'));
    handles.DCM.options.trials = trials;
    set(handles.Y1,'String',num2str(trials,'%7.0f'))
    m  = length(handles.DCM.options.trials);
catch
    m  = 1;
    set(handles.Y1,'String','1')
end

if isempty(m) || m == 0
    m  = 1;
    set(handles.Y1,'String','1');
end

handles = Xdefault(hObject,handles,m);

%-Get new trial data from file
%--------------------------------------------------------------------------
try
    D = spm_eeg_load(handles.DCM.xY.Dfile);
catch
    [f,p] = uigetfile({'*.mat'}, 'please select data file');

    if f == 0
        return;
    end

    handles.DCM.xY.Dfile = fullfile(p,f);
    
    D = spm_eeg_load(handles.DCM.xY.Dfile);    
end

[ok, D] = check(D, 'dcm');

if ~ok
    if check(D, 'basic')
        warndlg(['The requested file is not ready for DCM.'...
            'Use prep to specify sensors and fiducials or LFP channels.']);
    else
        warndlg('The meeg file is corrupt or incomplete');
    end
    
    handles.DCM.xY.Dfile = [];
    set(handles.data_ok, 'enable', 'off');
    guidata(hObject,handles);    
    return
end

set(handles.data_ok, 'enable', 'on');

% Assemble and display data
%--------------------------------------------------------------------------
handles     = reset_Callback(hObject, eventdata, handles);
handles.DCM = spm_dcm_erp_data(handles.DCM,handles.DCM.options.h);

set(handles.design,'enable', 'on')
set(handles.Uname, 'enable', 'on')
set(handles.Y,     'enable', 'on')

handles     = reset_Callback(hObject, eventdata, handles);
guidata(hObject,handles);


% --- Executes on button press in Y to display data
%--------------------------------------------------------------------------
function Y_Callback(hObject, eventdata, handles)
handles  = reset_Callback(hObject, eventdata, handles);
try
    handles.DCM = spm_dcm_erp_data(handles.DCM,handles.DCM.options.h);
    spm_dcm_erp_results(handles.DCM,'Data');
    set(handles.dt, 'String',sprintf('bins: %.1fms',handles.DCM.xY.dt*1000))
    set(handles.dt, 'Visible','on')
catch
    warndlg('please specify data and trials');
end
guidata(hObject,handles);

% --- Executes on button press in Uname.
%--------------------------------------------------------------------------
function Uname_Callback(hObject, eventdata, handles)
try
    m = length(handles.DCM.options.trials);
catch
    warndlg('please select trials')
    handles = Xdefault(hObject,handles,1);
    return
end
Str   = cellstr(get(handles.Uname,'String'));
n     = size(handles.DCM.xU.X,2);
for i = 1:n
    try
        Uname{i} = Str{i};
    catch
        Uname{i} = sprintf('effect %i',i);
    end
end
set(handles.Uname,'string',Uname)
handles.DCM.xU.name = Uname;
guidata(hObject,handles);

% --- Executes on button press in design.
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
    set(handles.design, 'String',num2str(handles.DCM.xU.X','%7.2f'));
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
set(handles.Uname,'string',Uname)
handles.DCM.xU.name = Uname;
guidata(hObject,handles);


% Spatial model specification
%==========================================================================
function handles = spatial_ok_Callback(hObject, eventdata, handles)
handles = reset_Callback(hObject, eventdata, handles);

% spatial model - source names
%--------------------------------------------------------------------------
Sname     = cellstr(get(handles.Sname,'String'));
Nareas    = length(Sname);
Nmodes    = get(handles.Nmodes,'Value');


% switch for spatial forward model (for EEG or MEG)
%--------------------------------------------------------------------------
DCM       = handles.DCM;
DCM.Sname = Sname;
DCM.options.type = get(handles.Spatial_type,'Value');
switch DCM.xY.modality
    
    case {'EEG','MEG'}

        % read location coordinates
        %------------------------------------------------------------------
        Slocation = zeros(Nareas, 3);
        tmp       = get(handles.Slocation, 'String');
        if ~isempty(tmp) & size(tmp,1) == Nareas
            for i = 1:Nareas
                tmp2 = str2num(tmp(i, :));
                if length(tmp2) ~= 3
                    errordlg(sprintf('coordinates of area %d not valid',i));
                    return;
                else
                    Slocation(i, :) = tmp2;
                end
            end
            if size(Slocation, 1) ~= Nareas
                errordlg('Number of source names and locations must correspond');
                return;
            end
        else
            errordlg(sprintf('Please specify %d source locations.', Nareas));
            return
        end

        % set prior expectations about locations
        %------------------------------------------------------------------
        DCM.Lpos = Slocation';
        
        % forward model (spatial)
        %--------------------------------------------------------------------------
        DCM = spm_dcm_erp_dipfit(DCM);

    case{'LFP'}
        
        % for LFP
        %------------------------------------------------------------------
        if length(DCM.xY.Ic) ~= Nareas
            errordlg('Number of source names and channels must correspond');
            return;
        end
        DCM.Lpos = zeros(3,0);
        set(handles.Slocation, 'String', '');    
     
    otherwise
        warndlg('Unknown data modlaity')
        return
        

end

handles.DCM = DCM;
set(handles.Spatial_type,     'Enable', 'off');
set(handles.spatial_ok,       'Enable', 'off');
set(handles.onset,            'Enable', 'off');
set(handles.Sname,            'Enable', 'off');
set(handles.Slocation,        'Enable', 'off');
set(handles.spatial_back,     'Enable', 'off');

set(handles.con_reset,        'Enable', 'on');
set(handles.connectivity_back,'Enable', 'on');
set(handles.Hz1,              'Enable', 'on');
set(handles.Hz2,              'Enable', 'on');
set(handles.Rft,              'Enable', 'on');


% [re]-set connections
%--------------------------------------------------------------------------
handles = connections_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

% --- Executes on button press in pos.
%--------------------------------------------------------------------------
function pos_Callback(hObject, eventdata, handles)
[f,p]     = uigetfile('*.mat','source (n x 3) location file');
Slocation = load(fullfile(p,f));
name      = fieldnames(Slocation);
Slocation = getfield(Slocation, name{1});
set(handles.Slocation,'String',num2str(Slocation,'%4.0f'));


%-Executes on button press in data_ok.
%--------------------------------------------------------------------------
function handles = data_ok_Callback(hObject, eventdata, handles)
handles      = reset_Callback(hObject, eventdata, handles);
handles.DCM  = spm_dcm_erp_data(handles.DCM,handles.DCM.options.h);


% enable next stage, disable data specification
%--------------------------------------------------------------------------
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
set(handles.onset,         'Enable', 'on');
set(handles.Sname,         'Enable', 'on');
set(handles.Slocation,     'Enable', 'on');
set(handles.spatial_back,  'Enable', 'on');
set(handles.spatial_ok,    'Enable', 'on');

% assume source and channel names are the same for LPF data
%--------------------------------------------------------------------------
if strcmp(handles.DCM.xY.modality,'LFP')
    Ic = handles.DCM.xY.Ic;
    if length(cellstr(get(handles.Sname,'String'))) ~= length(Ic);
        Sname = handles.DCM.xY.name;
        set(handles.Sname,'String',Sname)
        handles.DCM.Sname = Sname;
    end
end

guidata(hObject, handles);


% --- Executes on button press in spatial_back.
%----------------------------------------------------------------------
function spatial_back_Callback(hObject, eventdata, handles)

set(handles.Spatial_type, 'Enable', 'off');
set(handles.spatial_ok,   'Enable', 'off');
set(handles.onset,        'Enable', 'off');
set(handles.Sname,        'Enable', 'off');
set(handles.Slocation,    'Enable', 'off');
set(handles.spatial_back, 'Enable', 'off');

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

guidata(hObject, handles);

% --- Executes on button press in plot_dipoles.
%--------------------------------------------------------------------------
function plot_dipoles_Callback(hObject, eventdata, handles)

% read location coordinates
tmp = get(handles.Slocation, 'String');
Slocation = [];
if ~isempty(tmp) 
    for i = 1:size(tmp, 1)
        tmp2 = str2num(tmp(i, :))';
        if length(tmp2) == 3
            Slocation = [Slocation tmp2];
        end
    end
end

Nlocations   = size(Slocation, 2);
sdip.n_seeds = 1;
sdip.n_dip   = Nlocations;
sdip.Mtb     = 1;
sdip.j{1}    = zeros(3*Nlocations, 1);
sdip.loc{1}  = Slocation;
spm_eeg_inv_ecd_DrawDip('Init', sdip)

%-Connectivity
%==========================================================================

% Draw buttons
%--------------------------------------------------------------------------
function handles = connections_Callback(hObject, eventdata, handles)
handles = reset_Callback(hObject, eventdata, handles);
DCM     = handles.DCM;
n       = length(DCM.Sname);           % number of sources
m       = size(DCM.xU.X,2);            % number of experimental inputs
l       = length(DCM.options.onset);   % number of peristimulus inputs

% remove previous objects
%--------------------------------------------------------------------------
h = get(handles.SPM,'Children');
for i = 1:length(h)
    if strcmp(get(h(i),'Tag'),'tmp')
        delete(h(i));
    end
end


% no changes in coupling
%--------------------------------------------------------------------------
if ~m, B = {}; DCM.B = {}; end

% check DCM.A, DCM.B, ...
%--------------------------------------------------------------------------
try, if length(DCM.A{1}) ~= n, DCM = rmfield(DCM,'A'); end, end
try, if length(DCM.B{1}) ~= n, DCM = rmfield(DCM,'B'); end, end
try, if length(DCM.B)    ~= m, DCM = rmfield(DCM,'B'); end, end
try, if size(DCM.C,1)    ~= n, DCM = rmfield(DCM,'C'); end, end
try, if size(DCM.C,2)    ~= l, DCM = rmfield(DCM,'C'); end, end

% connection buttons
%--------------------------------------------------------------------------
set(handles.con_reset,'Units','Normalized')
p  = get(handles.con_reset,'Position');
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
                'Tag','tmp',...
                'Callback',str);
            
            % within region, between frequency coupling for 'Induced'
            %--------------------------------------------------------------
            if i == j
                set(A{k}(i,j),'Enable','off')
            end
            
            % allow nonlinear self-connections (induced responses)
            %--------------------------------------------------------------
            if get(handles.ERP,'value') == 6 && k == 2
                set(A{k}(i,j),'Enable','on')
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
                         'Tag','tmp',...
                         'Callback',str);
            
            % intrinsic modulation of H_e
            %--------------------------------------------------------------
            if i == j
                set(B{k}(i,j),'Enable','on')
            end
            try
                set(B{k}(i,j),'Value',DCM.B{k}(i,j));
            catch
                DCM.B{k}(i,j) = get(B{k}(i,j),'Value');
            end
        end
    end
end
for i = 1:n
    for j = 1:l
        x          = x0 + (4 - 1)*(n + 1)*sx +(j - 1)*sx;
        y          = y0 - (i + 4)*sy;
        str        = sprintf('data.DCM.C(%i,%i)',i,j);
        str        = ['data=guidata(gcbo);' str '=get(gcbo,''Value'');guidata(gcbo,data)'];
        C(i,j)     = uicontrol(handles.SPM,...
            'Units','Normalized',...
            'Position',[x y sx sy],...
            'Style','radiobutton',...
            'Tag','tmp',...
            'Callback',str);
        try
            set(C(i,j),'Value',DCM.C(i,j));
        catch
            DCM.C(i,j) = get(C(i,j),'Value');
        end
    end
end

% string labels
%--------------------------------------------------------------------------
if get(handles.ERP,'Value') ~= 6
    constr = {'A forward' 'A backward' 'A lateral' 'C input'};
else
    constr = {'A linear' 'A nonlinear' '(not used)' 'C input'};
end
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
        'Tag','tmp',...
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
        'Tag','tmp',...
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

% remove existing buttons and set DCM.A,.. to zero
%--------------------------------------------------------------------------
function con_reset_Callback(hObject, eventdata, handles)

h = get(handles.SPM,'Children');
for i = 1:length(h)
    if strcmp(get(h(i),'Tag'),'tmp')
        delete(h(i));
    end
end
try
    for i = 1:length(handles.DCM.A)
        handles.DCM.A{i}(:) = 0;
    end
    for i = 1:length(handles.DCM.B)
        handles.DCM.B{i}(:) = 0;
    end
    handles.DCM.C(:) = 0;
end
handles = connections_Callback(hObject, eventdata, handles);
guidata(hObject,handles)

% --- Executes on button press in connectivity_back.
%--------------------------------------------------------------------------
function connectivity_back_Callback(hObject, eventdata, handles)

set(handles.con_reset,         'Enable', 'off');
set(handles.connectivity_back, 'Enable', 'off');
set(handles.Hz1,               'Enable', 'off');
set(handles.Hz2,               'Enable', 'off');
set(handles.Rft,               'Enable', 'off');

set(handles.Spatial_type,      'Enable', 'on');
set(handles.spatial_ok,        'Enable', 'on');
set(handles.onset,             'Enable', 'on');
set(handles.Sname,             'Enable', 'on');
set(handles.Slocation,         'Enable', 'on');
set(handles.spatial_back,      'Enable', 'on');

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

% --- Executes on button press in estimate.
% -------------------------------------------------------------------------
function varargout = estimate_Callback(hObject, eventdata, handles, varargin)
set(handles.estimate,'String','Estimating','Foregroundcolor',[1 0 0])
handles = reset_Callback(hObject, eventdata, handles);

% initialise if required
% -------------------------------------------------------------------------
try
    Ep  = handles.DCM.Ep;
    Str = questdlg('initialize with previous estimates');
    if strcmp(Str,'Yes')
        handles.DCM.M.P = Ep;
    elseif strcmp(Str,'No')
        handles.DCM.M.P = [];
    elseif strcmp(Str,'Cancel')
        return
    end
end

% invert and save
%--------------------------------------------------------------------------
switch handles.DCM.options.model

    % conventional neural-mass and mean-field models
    %----------------------------------------------------------------------
    case{'ERP','SEP','NMM','MFM'}
        handles.DCM = spm_dcm_erp(handles.DCM);

    % Cross-spectral density model (steady-state responses)
    %----------------------------------------------------------------------
    case{'SSR'}
        handles.DCM = spm_dcm_ssr(handles.DCM);

    % Induced responses
    %----------------------------------------------------------------------
    case{'IND'}
        handles.DCM = spm_dcm_ind(handles.DCM);

    otherwise
        warndlg('unknown model type')
        return
end

handles = ERP_Callback(hObject, eventdata, handles);

set(handles.results,    'Enable','on' )
set(handles.save,       'Enable','on')
set(handles.estimate,   'String','Estimated','Foregroundcolor',[0 0 0])
if get(handles.Spatial_type,'Value') == 2
    set(handles.render, 'Enable','on' )
    set(handles.Imaging,'Enable','on' )
end

guidata(hObject, handles);


% --- Executes on button press in results.
% -------------------------------------------------------------------------
function varargout = results_Callback(hObject, eventdata, handles, varargin)
Action  = get(handles.results, 'String');
Action  = Action{get(handles.results, 'Value')};

switch handles.DCM.options.model

    % conventional neural-mass and mean-field models
    %----------------------------------------------------------------------
    case{'ERP','SEP','NMM','MFM'}
        spm_dcm_erp_results(handles.DCM, Action);

    % Cross-spectral density model (steady-state responses)
    %----------------------------------------------------------------------
    case{'SSR'}
        spm_dcm_ssr_results(handles.DCM, Action);

    % Induced responses
    %----------------------------------------------------------------------
    case{'IND'}
        spm_dcm_ind_results(handles.DCM, Action);

    otherwise
        warndlg('unknown model type')
        return
end


% --- Executes on button press in initialise.
% -------------------------------------------------------------------------
function initialise_Callback(hObject, eventdata, handles)

[f,p]           = uigetfile('DCM*.mat','please select estimated DCM');
DCM             = load(fullfile(p,f), '-mat');
handles.DCM.M.P = DCM.DCM.Ep;
guidata(hObject, handles);

% --- Executes on button press in render.
% -------------------------------------------------------------------------
function render_Callback(hObject, eventdata, handles)

spm_eeg_inv_visu3D_api(handles.DCM.xY.Dfile)

% --- Executes on button press in Imaging.
% -------------------------------------------------------------------------
function Imaging_Callback(hObject, eventdata, handles)

spm_eeg_inv_imag_api(handles.DCM.xY.Dfile)


%==========================================================================
function handles = Xdefault(hObject,handles,m)
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
set(handles.design,'String',num2str(handles.DCM.xU.X','%7.2f'));
set(handles.Uname, 'String',handles.DCM.xU.name);
return

% --- Executes on button press in BMC.
%--------------------------------------------------------------------------
function BMC_Callback(hObject, eventdata, handles)
spm_api_bmc

% --- Executes on selection change in ERP.
%--------------------------------------------------------------------------
function handles = ERP_Callback(hObject, eventdata, handles)

% get model type
%--------------------------------------------------------------------------
handles = reset_Callback(hObject, eventdata, handles);
switch handles.DCM.options.model

    % conventional neural-mass and mean-field models
    %----------------------------------------------------------------------
    case{'ERP','SEP','NMM','MFM'}
        Action = {
            'ERPs (mode)',
            'ERPs (sources)',
            'coupling (A)',
            'coupling (B)',
            'coupling (C)',
            'trial-specific effects',
            'Input',
            'Response',
            'Response (image)',
            'Dipoles',
            'Spatial overview'};
        try
            set(handles.Nmodes, 'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes, 'Value', 8);
        end
        set(handles.Spatial_type, 'String', {'ECD','Imaging'});
        set(handles.Wavelet,      'Enable','off','String','-');

    % Cross-spectral density model (steady-state responses)
    %----------------------------------------------------------------------
    case{'SSR'}
        Action = {
            'Data',
            'coupling (A)',
            'coupling (B)',
            'coupling (C)',
            'trial-specific effects',
            'Input',
            'Cross-spectral density',
            'Dipoles'};
        try
            set(handles.Nmodes,   'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes,   'Value', 4);
        end
        set(handles.Spatial_type, 'Value', 1);
        set(handles.Spatial_type, 'String', {'ECD'});
        set(handles.Wavelet,      'Enable','on','String','Spectral density');

    % Induced responses
    %----------------------------------------------------------------------
    case{'IND'}
        Action = {
            'Frequency modes'
            'Time-modes'
            'Time-frequency'
            'Coupling (A - Hz)'
            'Coupling (B - Hz)'
            'Coupling (A - modes)'
            'Coupling (B - modes)'
            'Input (C - Hz)'
            'Input (u - ms)'
            'Dipoles'};
        try
            set(handles.Nmodes,   'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes,   'Value', 4);
        end
        set(handles.Spatial_type, 'Value', 1);
        set(handles.Spatial_type, 'String', {'ECD'});
        set(handles.Wavelet,      'Enable','on','String','Wavelet transform');

    otherwise
        warndlg('unknown model type')
        return
end

set(handles.results,'Value',1);
set(handles.results,'String',Action);
handles = reset_Callback(hObject, eventdata, handles);
guidata(hObject,handles);

% --- Executes on button press in Wavelet.
function Wavelet_Callback(hObject, eventdata, handles)

% get transform
%--------------------------------------------------------------------------
handles     = reset_Callback(hObject, eventdata, handles);

switch handles.DCM.options.model

    case{'SSR'}
        
        % cross-spectral density (if DCM.M.U (eigenspace) exists
        %------------------------------------------------------------------
        try
            handles.DCM = spm_dcm_ssr_data(handles.DCM);
        end

        % and display
        %------------------------------------------------------------------
        spm_dcm_ssr_results(handles.DCM,'Data');


    case{'IND'}
        
        % wavelet tranform
        %------------------------------------------------------------------
        handles.DCM = spm_dcm_ind_data(handles.DCM);

        % and display
        %------------------------------------------------------------------
        spm_dcm_ind_results(handles.DCM,'Wavelet');

end
guidata(hObject,handles);



