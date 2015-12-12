function varargout = spm_api_erp(varargin)
% SPM_API_ERP Application M-file for spm_api_erp.fig
%    FIG = SPM_API_ERP launch spm_api_erp GUI.
%    SPM_API_ERP('callback_name', ...) invoke the named callback.
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_api_erp.m 6644 2015-12-12 14:53:37Z karl $
 

%-Launch GUI
%--------------------------------------------------------------------------
if nargin == 0 || nargin == 1
 
    fig     = openfig(mfilename,'reuse');
    S0      = spm('WinSize','0',1);
    set(fig,'units','pixels'); Fdim = get(fig,'position');
    set(fig,'position',[S0(1) S0(2) 0 0] + Fdim);
    Fgraph  = spm_figure('GetWin','Graphics');
    
    % Use system color scheme for figure:
    %----------------------------------------------------------------------
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
 
    % Generate a structure of handles to pass to callbacks, and store it
    %----------------------------------------------------------------------
    handles  = guihandles(fig);
    handles.Fgraph = Fgraph;
    guidata(fig, handles);
    
    if nargin == 1
        load_Callback(fig, [], handles, varargin{1})
    end
 
    if nargout > 0
        varargout{1} = fig;
    end

%-Invoke named subfunction or callback
%--------------------------------------------------------------------------
elseif ischar(varargin{1})
    try
        if nargout
            [varargout{1:nargout}] = feval(varargin{:});
        else
            feval(varargin{:});
        end
    catch
        disp(lasterr);
    end
end
 
 
%-DCM files and directories: Load and save
%==========================================================================
 
% --- Executes on button press in load.
% -------------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles, varargin)
try
    DCM     = varargin{1};
    f       = spm_file(DCM.name,'filename');
catch
    name    = spm_select(1,'mat','please select DCM file');
    
    f       = spm_file(name,'filename');
    DCM     = load(f,'-mat');
    try
        DCM = DCM.DCM;
    catch
        error('File "%s" does not contain a DCM structure.',f);
    end
    DCM.name    = name;
    handles.DCM = DCM;
    cd(spm_file(name,'path'))
    guidata(hObject,handles);
end
 
% Type of analysis
%--------------------------------------------------------------------------
% 'ERP'    - Event related responses
% 'CSD'    - Cross-spectral density
% 'TFM'    - Time-frequency responses
% 'IND'    - Induced responses
% 'PHA'    - (for phase coupling)
% 'NFM'    - Neural field model

try
    model = DCM.options.analysis;
catch
    model = get(handles.ERP,'String');
    model = model{get(handles.ERP,'Value')};
    DCM.options.analysis = model;
end
switch model
    case{'ERP'}, set(handles.ERP,  'Value',1);
    case{'CSD'}, set(handles.ERP,  'Value',2);
    case{'TFM'}, set(handles.ERP,  'Value',3); 
                 set(handles.model,'Value',7);
    case{'IND'}, set(handles.ERP,  'Value',4);
    case{'PHA'}, set(handles.ERP,  'Value',5);
    case{'NFM'}, set(handles.ERP,  'Value',6);
                 set(handles.model,'Value',10);
    otherwise
end
handles = ERP_Callback(hObject, eventdata, handles);
 
 
% Type of model (neuronal)
%--------------------------------------------------------------------------
% 'ERP'    - (linear second order NMM slow)
% 'SEP'    - (linear second order NMM fast)
% 'CMC'    - (linear second order NMM Canonical microcircuit)
% 'LFP'    - (linear second order NMM self-inhibition)
% 'NMM'    - (nonlinear second order NMM first-order moments)
% 'MFM'    - (nonlinear second order NMM second-order moments)
% 'CMM'    - (nonlinear second order NMM Canonical microcircuit)
% 'DEM'    - (functional architecture based on a DEM scheme)
% 'NMDA'    - (nonlinear second order NMM first-order moments with NMDA receptors)
% 'CMM_NMM'- (nonlinear first order Canonical microcircuit with NMDA receptors)

try
    model = DCM.options.model;
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
    case{'CMM'}, set(handles.model,'Value',7);
    case{'NMDA'}, set(handles.model,'Value',8);
    case{'CMM_NMDA'}, set(handles.model,'Value',9);
    case{'NFM'}, set(handles.model,'Value',10);

    otherwise
end

% Type of model (spatial)
%--------------------------------------------------------------------------
% 'ECD'    - Equivalent current dipole
% 'IMG'    - Imaging
% 'LFP'    - Local field potentials
% 'ITR'    - Intra-Laminar recording
try
    model = DCM.options.spatial;
catch
    model = get(handles.Spatial,'String');
    model = model{get(handles.Spatial,'Value')};
    DCM.options.spatial = model;
end

if ismember(DCM.options.analysis, {'IND', 'PHA'})
    switch model
        case{'IMG'}, set(handles.Spatial,'Value',1);
        case{'ECD'}, set(handles.Spatial,'Value',1);
        case{'LFP'}, set(handles.Spatial,'Value',2);
        case{'ITR'}, set(handles.Spatial,'Value',4);
        otherwise
    end
elseif ismember(DCM.options.analysis, {'NFM'})
    switch model
        case{'LFP'}, set(handles.Spatial,'Value',1);
        otherwise
    end
else
    switch model
        case{'IMG'}, set(handles.Spatial,'Value',1);
        case{'ECD'}, set(handles.Spatial,'Value',2);
        case{'LFP'}, set(handles.Spatial,'Value',3);
        case{'ITR'}, set(handles.Spatial,'Value',4);
        otherwise
    end
end

% Filename
%--------------------------------------------------------------------------
try, set(handles.name,'String',f);  end
 
% Source location
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
try, set(handles.Nmodes,       'Value', DCM.options.Nmodes);         end
try, set(handles.h,            'Value', ...
        find(str2double(get(handles.h,'String')) == DCM.options.h)); end
try, set(handles.han,         'Value', DCM.options.han);             end
try, set(handles.D,           'Value', DCM.options.D);               end
try, set(handles.lock,        'Value', DCM.options.lock);            end
try, set(handles.multiC,      'Value', DCM.options.multiC);          end
try, set(handles.location,    'Value', DCM.options.location);        end
try, set(handles.symmetry,    'Value', DCM.options.symmetry);        end
try, set(handles.design,      'String',num2str(DCM.xU.X','%7.2f'));  end
try, set(handles.Uname,       'String',DCM.xU.name);                 end
try, set(handles.Sname,       'String',DCM.Sname);                   end
try, set(handles.onset,       'String',num2str(DCM.options.onset));  end
try, set(handles.dur,         'String',num2str(DCM.options.dur));    end
try, set(handles.Slocation,   'String',num2str(DCM.Lpos','%4.0f'));  end
 
% Imaging
%--------------------------------------------------------------------------
switch DCM.options.spatial
    case{'IMG'}
        set(handles.Imaging,'Enable','on' )
    otherwise
        set(handles.Imaging,'Enable','off' )
end
 
handles.DCM = DCM;
guidata(hObject, handles);

% estimation and results
%--------------------------------------------------------------------------
try
    handles.DCM.F;
    set(handles.results,    'Enable','on');
    switch handles.DCM.options.spatial
        case{'IMG'}
            set(handles.Imaging,'Enable','on');
        otherwise
            set(handles.Imaging,'Enable','off');
    end
catch
    set(handles.results,    'Enable','off');
end
guidata(hObject, handles);

% data & design specification
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
 

 
% --- Executes on button press in save.
% -------------------------------------------------------------------------
function handles = save_Callback(hObject, eventdata, handles)

handles = reset_Callback(hObject, eventdata, handles);
try
   file      = spm_file(handles.DCM.name,'filename');
catch
    try
        file = spm_file(handles.DCM.xY.Dfile,'filename');
        file = ['DCM_' file];
    catch
        file = ['DCM' date];
    end
end
[file,fpath] = uiputfile('DCM*.mat','DCM file to save',file);
 
if fpath
    handles.DCM.name = fullfile(fpath,file);
    set(handles.name,'String',file);
    DCM              = handles.DCM;
    save(DCM.name,'DCM', spm_get_defaults('mat.format'))
    set(handles.estimate,   'Enable', 'on')
    set(handles.initialise, 'Enable', 'on');
    cd(fpath)
end
 
% assign in base
%--------------------------------------------------------------------------
assignin('base','DCM',handles.DCM)
guidata(hObject,handles);
 
% store selections in DCM
% -------------------------------------------------------------------------
function handles = reset_Callback(hObject, eventdata, handles)
 
% analysis type
%--------------------------------------------------------------------------
model = get(handles.ERP,       'String');
model = model{get(handles.ERP, 'Value')};
handles.DCM.options.analysis = model;
 
if isequal(model, 'ERP')
    set(handles.han, 'Enable', 'on');
else
    set(handles.han, 'Value', 0);
    set(handles.han, 'Enable', 'off');
end

if isequal(model, 'PHA')
    set(handles.text20, 'String', 'sub-trials');
else
    set(handles.text20, 'String', 'modes');
end

% model type
%--------------------------------------------------------------------------
model = get(handles.model,         'String');
model = model{get(handles.model,   'Value')};
handles.DCM.options.model     = model;
        
% spatial type
%--------------------------------------------------------------------------
model = get(handles.Spatial,       'String');
model = model{get(handles.Spatial, 'Value')};
handles.DCM.options.spatial  = model;

handles.DCM.options.trials   = str2num(get(handles.Y1,    'String'));
handles.DCM.options.Tdcm(1)  = str2num(get(handles.T1,    'String'));
handles.DCM.options.Tdcm(2)  = str2num(get(handles.T2,    'String'));
handles.DCM.options.Fdcm(1)  = str2num(get(handles.Hz1,   'String'));
handles.DCM.options.Fdcm(2)  = str2num(get(handles.Hz2,   'String'));
handles.DCM.options.Rft      = str2num(get(handles.Rft,   'String'));
handles.DCM.options.onset    = str2num(get(handles.onset, 'String'));
handles.DCM.options.dur      = str2num(get(handles.dur,   'String'));
handles.DCM.options.Nmodes   = get(handles.Nmodes,        'Value');
detrend_val                  = str2double(get(handles.h,  'String'));
handles.DCM.options.h        = detrend_val(get(handles.h, 'Value'));
handles.DCM.options.han      = get(handles.han,           'Value');
handles.DCM.options.D        = get(handles.D,             'Value');
handles.DCM.options.lock     = get(handles.lock,          'Value');
handles.DCM.options.multiC   = get(handles.multiC,        'Value');
handles.DCM.options.location = get(handles.location,      'Value');
handles.DCM.options.symmetry = get(handles.symmetry,      'Value');
 
guidata(hObject,handles);
 
 
% Data selection and design
%==========================================================================
 
% --- Executes on button press in Datafile.
%--------------------------------------------------------------------------
function Datafile_Callback(hObject, eventdata, handles)
 
%-Get trials and data
%--------------------------------------------------------------------------
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
[f,p] = uigetfile({'*.mat'}, 'please select data file');
if f  == 0, return; end

handles.DCM.xY = [];
handles.DCM.xY.Dfile = fullfile(p,f);
D = spm_eeg_load(handles.DCM.xY.Dfile);

[D, ok] = check(D, 'dcm');

if ~ok
    warndlg(['The requested file is not ready for DCM.'...
        'Use prep to specify sensors and fiducials or LFP channels.']);
    
    handles.DCM.xY.Dfile = [];
    set(handles.data_ok, 'enable', 'off');
    guidata(hObject,handles);
    return
end

[mod, list] = modality(D, 0, 1);

if ismember('MEGCOMB', list)
    list = setdiff(list, 'MEGCOMB');
    if isempty(list)
        errordlg('MEGCOMB modality cannot be used for DCM', 'Error');
        return;
    elseif numel(list) == 1
        mod = list{1};
    end
end

if isequal(mod, 'Multimodal')
    qstr = 'Only one modality can be modelled at a time. Please select.';
    if numel(list) < 4
        
        % Nice looking dialog. Will usually be OK
        %------------------------------------------------------------------
        options = [];
        options.Default = list{1};
        options.Interpreter = 'none';
        handles.DCM.xY.modality = questdlg(qstr, 'Select modality', list{:}, options);
        
    else
        % Ugly but can accomodate more buttons
        %------------------------------------------------------------------
        ind = menu(qstr, list);
        handles.DCM.xY.modality = list{ind};
        
    end
else
    handles.DCM.xY.modality = mod;
end

if isequal(handles.DCM.xY.modality, 'LFP')
    set(handles.Spatial, 'Value', find(strcmp('LFP', get(handles.Spatial, 'String'))));
end  

if isequal(handles.DCM.xY.modality, 'ILAM')
    set(handles.Spatial, 'Value', find(strcmp('ITR', get(handles.Spatial, 'String'))));
end 

% Assemble and display data
%--------------------------------------------------------------------------
handles = reset_Callback(hObject, eventdata, handles);
try
    handles.DCM  = spm_dcm_erp_data(handles.DCM);
    if isfield(handles.DCM.xY, 'y') && ~isfield(handles.DCM.xY, 'xf')
        spm_dcm_erp_results(handles.DCM, 'Data');
    else
        spm_dcm_ind_results(handles.DCM, 'Wavelet');
    end  
    set(handles.dt, 'String',sprintf('bins: %.1fms', handles.DCM.xY.dt*1000))
    set(handles.dt, 'Visible','on')
    set(handles.data_ok, 'enable', 'on'); 
    guidata(hObject,handles);
catch
    errordlg({'please ensure trial selection and data are consistent';
             'data have not been changed'});
    set(handles.data_ok, 'enable', 'off'); 
end

set(handles.design,'enable', 'on')
set(handles.Uname, 'enable', 'on')
set(handles.Y,     'enable', 'on')
 
guidata(hObject,handles);
 
 
% --- Executes on button press in Y to display data
%--------------------------------------------------------------------------
function Y_Callback(hObject, eventdata, handles)
handles = reset_Callback(hObject, eventdata, handles);
try
    DCM  = spm_dcm_erp_data(handles.DCM);
    if isfield(handles.DCM.xY, 'y') && ~isfield(handles.DCM.xY, 'xf')
        spm_dcm_erp_results(DCM, 'Data');
    else
        spm_dcm_ind_results(DCM, 'Wavelet');
    end  
    set(handles.dt, 'String',sprintf('bins: %.1fms', DCM.xY.dt*1000))
    set(handles.dt, 'Visible','on')
    set(handles.data_ok, 'enable', 'on'); 
    guidata(hObject,handles);
    try
        str = ['trials (' num2str(DCM.xY.nt) ')'];
        set(handles.text22, 'String', str);
    end
    
catch
    errordlg({'please ensure trial selection and data are consistent';
             'data have not been changed'});
    set(handles.data_ok, 'enable', 'off'); 
end

 
% --- Executes on button press in Uname.
%--------------------------------------------------------------------------
function  Uname_Callback(hObject, eventdata, handles)
handles = reset_Callback(hObject, eventdata, handles);
try
    m   = length(handles.DCM.options.trials);
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
handles = reset_Callback(hObject, eventdata, handles);
try
    m   = length(handles.DCM.options.trials);
catch
    handles = Xdefault(hObject,handles,1);
    warndlg('please select trials')
    return
end
try
    X = str2num(get(handles.design,'String'))';
    if isempty(X)
        X = sparse(m,0);
    end
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
 
 
 
%-Executes on button press in data_ok.
%--------------------------------------------------------------------------
function handles = data_ok_Callback(hObject, eventdata, handles)
 
%-assemble and display trials.
%--------------------------------------------------------------------------
Y_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
switch get(handles.data_ok,'enable')
    case{'off'}
        return
end

%-check trial-specific effects
%--------------------------------------------------------------------------
design_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

% enable next stage, disable data specification
%--------------------------------------------------------------------------
set(handles.ERP,           'Enable', 'off');
set(handles.model,         'Enable', 'off');
set(handles.Datafile,      'Enable', 'off');
set(handles.Y1,            'Enable', 'off');
set(handles.T1,            'Enable', 'off');
set(handles.T2,            'Enable', 'off');
set(handles.Nmodes,        'Enable', 'off');
set(handles.h,             'Enable', 'off');
set(handles.D,             'Enable', 'off');
set(handles.design,        'Enable', 'off');
set(handles.Uname,         'Enable', 'off');
set(handles.data_ok,       'Enable', 'off');
 
set(handles.Spatial,       'Enable', 'on');
set(handles.plot_dipoles,  'Enable', 'on');
set(handles.Sname,         'Enable', 'on');
set(handles.Slocation,     'Enable', 'on');
set(handles.spatial_back,  'Enable', 'on');
set(handles.spatial_ok,    'Enable', 'on');

switch handles.DCM.options.analysis
    case{'CSD'}
        set(handles.onset, 'Enable', 'off');
        set(handles.dur,   'Enable', 'off');
    otherwise
        set(handles.onset, 'Enable', 'on');
        set(handles.dur,   'Enable', 'on');
end

 
guidata(hObject, handles);
 
% spatial model specification
%==========================================================================
function handles = spatial_ok_Callback(hObject, eventdata, handles)
handles = reset_Callback(hObject, eventdata, handles);
 
% spatial model - source names
%--------------------------------------------------------------------------
Sname     = cellstr(get(handles.Sname,'String'));
Nareas    = length(Sname);
Nchannels = length(handles.DCM.xY.Ic);
 
 
% switch for spatial forward model (for EEG or MEG)
%--------------------------------------------------------------------------
DCM       = handles.DCM;
DCM.Sname = Sname;
switch DCM.options.spatial
    
    case {'ECD','IMG'}
 
        % read location coordinates
        %------------------------------------------------------------------
        Slocation = zeros(Nareas, 3);
        tmp       = get(handles.Slocation, 'String');
        if ~isempty(tmp) && size(tmp,1) == Nareas
            for i = 1:Nareas
                tmp2 = str2num(tmp(i, :));
                if length(tmp2) ~= 3
                    errordlg(sprintf('coordinates of source %d invalid',i));
                    return;
                else
                    Slocation(i, :) = tmp2;
                end
            end
            if size(Slocation, 1) ~= Nareas
                errordlg('Number of sources and locations must correspond');
                return;
            end
        else
            errordlg(sprintf('Please specify %d source locations.', Nareas));
            return
        end
 
        % set prior expectations on location
        %------------------------------------------------------------------
        DCM.Lpos = Slocation';
        
        % forward model (spatial)
        %------------------------------------------------------------------
        try
            DCM = spm_dcm_erp_dipfit(DCM);
        end
        set(handles.plot_dipoles,'enable','on')
        
 
    case{'LFP','ITR'}
          
        if ~isdeployed, 
            addpath(fullfile(spm('Dir'),'toolbox','Neural_Models'));
        end
        
        % for LFP
        %------------------------------------------------------------------
        DCM.Lpos = zeros(3,0);
        
        if numel(Sname) < Nchannels
            warndlg('There are more LFP channels than sources')
            return
        end
        
        set(handles.Slocation, 'String', sprintf('%i electrodes',Nchannels));
        set(handles.plot_dipoles,'enable','off')              
       
     
    otherwise
        warndlg('Unknown data modality')
        return
        
end
 
handles.DCM = DCM;
set(handles.Spatial,          'Enable', 'off');
set(handles.data_ok,          'Enable', 'off');
set(handles.spatial_ok,       'Enable', 'off');
set(handles.onset,            'Enable', 'off');
set(handles.dur,              'Enable', 'off');
set(handles.Sname,            'Enable', 'off');
set(handles.Slocation,        'Enable', 'off');
set(handles.spatial_back,     'Enable', 'off');
 
set(handles.con_reset,        'Enable', 'on');
set(handles.save_spec,        'Enable', 'on');
set(handles.reduce,           'Enable', 'on');
set(handles.priors,           'Enable', 'on');
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
 
 
% --- Executes on button press in spatial_back.
%----------------------------------------------------------------------
function spatial_back_Callback(hObject, eventdata, handles)
 
set(handles.Spatial,      'Enable', 'off');
set(handles.spatial_ok,   'Enable', 'off');
set(handles.onset,        'Enable', 'off');
set(handles.dur,          'Enable', 'off');
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
set(handles.ERP,          'Enable', 'on');
set(handles.Datafile,     'Enable', 'on');

guidata(hObject, handles);

ERP_Callback(hObject, eventdata, handles);   
 
% --- Executes on button press in plot_dipoles.
%--------------------------------------------------------------------------
function plot_dipoles_Callback(hObject, eventdata, handles)
 
% read location coordinates
%--------------------------------------------------------------------------
tmp       = get(handles.Slocation, 'String');
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

% numbers of buttons
%--------------------------------------------------------------------------
n     = length(DCM.Sname);              % number of sources
m     = size(DCM.xU.X,2);               % number of experimental inputs
l     = length(DCM.options.onset);      % number of peristimulus inputs
nk    = 3;                              % number of connection types
nj    = ones(nk,1)*n;                   % number of sources

if strcmpi(DCM.options.analysis,'CSD'), l = 0; end


% remove previous objects
%--------------------------------------------------------------------------
h     = get(handles.SPM,'Children');
for i = 1:length(h)
    if strcmp(get(h(i),'Tag'),'tmp')
        delete(h(i));
    end
end
 
 
% no changes in coupling
%--------------------------------------------------------------------------
if ~m, B = {}; DCM.B = {}; end
if ~l, C = {}; DCM.C = {}; end
 
% check DCM.A, DCM.B, ...
%--------------------------------------------------------------------------
try, if size(DCM.A{1},1) ~= n, DCM = rmfield(DCM,'A'); end, end
try, if size(DCM.B{1},1) ~= n, DCM = rmfield(DCM,'B'); end, end
try, if numel(DCM.B)     ~= m, DCM = rmfield(DCM,'B'); end, end
try, if size(DCM.C,1)    ~= n, DCM = rmfield(DCM,'C'); end, end
try, if size(DCM.C,2)    ~= l, DCM = rmfield(DCM,'C'); end, end


 
% connection buttons (A)
%--------------------------------------------------------------------------
set(handles.con_reset,'Units','Normalized')
p  = get(handles.con_reset,'Position');
x0 = 0.1;
y0 = 0.44;
sx = 1/36;
sy = 1/72;
for i = 1:n
    for k = 1:nk
        for j = 1:nj(k)
            x          = x0 + (j - 1)*sx + (sum(nj(1:(k - 1))) + k - 1)*sx;
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

            % allow nonlinear self-connections when appriopriate
            %--------------------------------------------------------------
            if strcmpi(DCM.options.analysis,'IND') && k == 2
                set(A{k}(i,j),'Enable','on')
            end
            if strcmpi(DCM.options.model,'DEM')
                set(A{k}(i,j),'Enable','on')
            end
            if strcmpi(DCM.options.model,'CMC') && k ==3
                set(A{k}(i,j),'Enable','on')
            end
            
            % Fill in values if specified
            %--------------------------------------------------------------
            try
                set(A{k}(i,j),'Value',DCM.A{k}(i,j));
            catch
                DCM.A{k}(i,j) = get(A{k}(i,j),'Value');
            end
        end
    end
end

% connection buttons (B)
%--------------------------------------------------------------------------
for i = 1:n
    for k = 1:m
        for j = 1:n
            x          = x0 + (j - 1)*sx + ((k - 1)*n + k - 1)*sx;
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

% connection buttons (C)
%--------------------------------------------------------------------------
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
switch DCM.options.model
    case{'NMM','MFM','NMDA'}
        constr = {'Excit.' 'Inhib.'    'Mixed'       'input'};
    case{'CMC'}
        constr = {'forward' 'back'     'Modulatory'  'input'};
    case{'DEM'}
        constr = {'States' ' '         ' '           ' '    };
    otherwise
        constr = {'forward' 'back'     'lateral'     'input'};
end
switch DCM.options.analysis
    case{'CSD'}
        constr{4} = ' ';
    case{'IND'}
        constr = {'linear' 'nonlinear' '(not used)'  'input'};
    case{'PHA'}
        constr = {'endog' '(not used)' '(not used)'  '(not used)' };
    otherwise
end
nsx            = (n + 1)*sx;
nsy            = 2*sy;
for k = 1:4
    x          = x0 + (k - 1)*nsx;
    y          = y0 - 4*sy;
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
set(handles.save_spec,         'Enable', 'off');
set(handles.reduce,            'Enable', 'off');
set(handles.connectivity_back, 'Enable', 'off');
set(handles.Hz1,               'Enable', 'off');
set(handles.Hz2,               'Enable', 'off');
set(handles.Rft,               'Enable', 'off');
 
set(handles.Spatial,           'Enable', 'on');
set(handles.spatial_ok,        'Enable', 'on');
set(handles.onset,             'Enable', 'on');
set(handles.dur,               'Enable', 'on');
set(handles.Sname,             'Enable', 'on');
set(handles.Slocation,         'Enable', 'on');
set(handles.spatial_back,      'Enable', 'on');

switch handles.DCM.options.analysis
    case{'CSD'}
        set(handles.onset,     'Enable', 'off');
        set(handles.dur,       'Enable', 'off');
    otherwise
        set(handles.onset,     'Enable', 'on');
        set(handles.dur,       'Enable', 'on');
end
 
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
 
%-Estimate, initialise and review
%==========================================================================
 
% --- Executes on button press in estimate.
% -------------------------------------------------------------------------
function varargout = estimate_Callback(hObject, eventdata, handles, varargin)
set(handles.estimate,'String','Estimating','Foregroundcolor',[1 0 0])
handles = reset_Callback(hObject, eventdata, handles);
 
% initialise with posteriors if required
% -------------------------------------------------------------------------
try
    Ep  = handles.DCM.Ep;
    Str = questdlg('initialise with previous posteriors');
    if strcmp(Str,'Yes')
        handles.DCM.M.P = Ep;
    elseif strcmp(Str,'No')
        handles.DCM.M.P = [];
    elseif strcmp(Str,'Cancel')
        return
    end
end

% previous priors
% -------------------------------------------------------------------------
try
    handles.DCM.M.pE;
    Str = questdlg('use previous priors');
    if strcmp(Str,'No')
        handles.DCM.M = rmfield(handles.DCM.M,{'pE','pC'});
        try
            handles.DCM.M = rmfield(handles.DCM.M,{'gE','gC'});
        end
    elseif strcmp(Str,'Cancel')
        return
    end
end

% previous hyperpriors
% -------------------------------------------------------------------------
try
    handles.DCM.M.hE;
    Str = questdlg('use previous hyperpriors');
    if strcmp(Str,'No')
        handles.DCM.M = rmfield(handles.DCM.M,{'hE','hC'});
    elseif strcmp(Str,'Cancel')
        return
    end
end
 
% invert and save
%--------------------------------------------------------------------------
switch handles.DCM.options.analysis

    % conventional neural-mass and mean-field models
    %----------------------------------------------------------------------
    case{'ERP'}

        switch handles.DCM.options.model
            
            % DEM
            %----------------------------------------------------------------------
            case{'DEM'}
                handles.DCM = spm_dcm_dem(handles.DCM);

            otherwise
                handles.DCM = spm_dcm_erp(handles.DCM);
        end
        
    % cross-spectral density model (complex)
    %----------------------------------------------------------------------
    case{'CSD'}
        handles.DCM = spm_dcm_csd(handles.DCM);

    % cross-spectral density model (steady-state responses)
    %----------------------------------------------------------------------
    case{'TFM'}
        handles.DCM = spm_dcm_tfm(handles.DCM);

    % induced responses
    %----------------------------------------------------------------------
    case{'IND'}
        handles.DCM = spm_dcm_ind(handles.DCM);

    % phase coupling
    %----------------------------------------------------------------------
    case{'PHA'}
        handles.DCM = spm_dcm_phase(handles.DCM);

    % cross-spectral density model (steady-state responses)
    %----------------------------------------------------------------------
    case{'NFM'}
        handles.DCM = spm_dcm_nfm(handles.DCM);

    otherwise
        warndlg('unknown analysis type')
        return
end

handles = ERP_Callback(hObject, eventdata, handles);

set(handles.results,    'Enable','on' )
set(handles.save,       'Enable','on')
set(handles.estimate,   'String','Estimated','Foregroundcolor',[0 0 0])
if get(handles.Spatial, 'Value') == 1
    set(handles.Imaging,'Enable','on' )
end
 
guidata(hObject, handles);
 
 
% --- Executes on button press in results.
% -------------------------------------------------------------------------
function varargout = results_Callback(hObject, eventdata, handles, varargin)
Action  = get(handles.results, 'String');
Action  = Action{get(handles.results, 'Value')};
 
switch handles.DCM.options.analysis
 
    % conventional neural-mass and mean-field models
    %----------------------------------------------------------------------
    case{'ERP'}
        spm_dcm_erp_results(handles.DCM, Action);
        
    % Cross-spectral density model (complex)
    %----------------------------------------------------------------------
    case{'CSD'}
        spm_dcm_csd_results(handles.DCM, Action);
 
    % Cross-spectral density model (steady-state responses)
    %----------------------------------------------------------------------
    case{'TFM'}
        spm_dcm_tfm_results(handles.DCM, Action);
 
    % Induced responses
    %----------------------------------------------------------------------
    case{'IND'}
        spm_dcm_ind_results(handles.DCM, Action);

    % phase coupling
    %----------------------------------------------------------------------
    case{'PHA'}
        spm_dcm_phase_results(handles.DCM, Action); 
 
    % Cross-spectral density model (complex)
    %----------------------------------------------------------------------
    case{'NFM'}
        spm_dcm_csd_results(handles.DCM, Action);
        
    otherwise
        warndlg('unknown analysis type')
        return
end
 
 
% --- Executes on button press in initialise.
% -------------------------------------------------------------------------
function initialise_Callback(hObject, eventdata, handles)
 
[f,p]           = uigetfile('DCM*.mat','please select estimated DCM');
DCM             = load(fullfile(p,f), '-mat');
handles.DCM.M.P = DCM.DCM.Ep;
guidata(hObject, handles);
 

% --- Executes on button press in Imaging.
% -------------------------------------------------------------------------
function Imaging_Callback(hObject, eventdata, handles)
 
spm_eeg_inv_imag_api(handles.DCM.xY.Dfile)
 
 
% default design matrix
%==========================================================================
function handles = Xdefault(hObject,handles,m)
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
function BMS_Callback(hObject, eventdata, handles)
%spm_api_bmc
spm_jobman('Interactive','','spm.dcm.bms.inference');
 
 
% --- Executes on selection change in ERP.
%--------------------------------------------------------------------------
function handles = ERP_Callback(hObject, eventdata, handles)
 
% get analysis type
%--------------------------------------------------------------------------
handles = reset_Callback(hObject, eventdata, handles);
switch handles.DCM.options.analysis
 
    % conventional neural-mass and mean-field models
    %----------------------------------------------------------------------
    case {'ERP'}
        Action = {
            'ERPs (mode)',...
            'ERPs (sources)',...
            'coupling (A)',...
            'coupling (B)',...
            'coupling (C)',...
            'trial-specific effects',...
            'Input',...
            'Response',...
            'Response (image)',...
            'Scalp maps',...
            'Dipoles'};
        try
            set(handles.Nmodes, 'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes, 'Value', 8);
        end
        
        set(handles.text20,     'String', 'modes');
        set(handles.model,      'Enable','on');
        set(handles.Spatial,    'String',{'IMG','ECD','LFP','ITR'});
        set(handles.Wavelet,    'Enable','off');
        set(handles.onset,      'Enable','on');
        set(handles.dur,        'Enable','on');
        
    % Cross-spectral density model (complex)
    %----------------------------------------------------------------------
    case {'CSD'}
        Action = {
              'spectral data',...
              'Coupling (A)',...
              'Coupling (B)',...
              'Coupling (C)',...
              'trial-specific effects',...
              'Input',...
              'Transfer functions',...
              'Cross-spectra (sources)',...
              'Cross-spectra (channels)',...
              'Coherence (sources)',...
              'Coherence (channels)',...
              'Covariance (sources)',...
              'Covariance (channels)',...
              'Dipoles'};
        try
            set(handles.Nmodes, 'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes, 'Value', 4);
        end
        
        set(handles.text20, 'String', 'modes');
        set(handles.model,  'Enable','on');              
        set(handles.Spatial,'String',{'IMG','ECD','LFP','ITR'});
        set(handles.Wavelet,'Enable','on');
        set(handles.onset,  'Enable','off');
        set(handles.dur,    'Enable','off');

 
    % Cross-spectral density model (steady-state responses)
    %----------------------------------------------------------------------
    case{'TFM'}
        Action = {
              'induced responses',...
              'induced and evoked responses',...
              'Coupling (A)',...
              'Coupling (B)',...
              'Coupling (C)',...
              'trial-specific effects',...
              'Endogenous input',...
              'Exogenous input',...
              'Transfer functions',...
              'induced predictions',...
              'induced and evoked predictions',...
              'induced predictions - sources',...
              'induced and evoked predictions - sources',...
              'Dipoles'};
 
        try
            set(handles.Nmodes, 'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes, 'Value', 4);
        end
        
        set(handles.text20, 'String','modes');            
        set(handles.Spatial,'String',{'IMG','ECD','LFP','ITR'});
        set(handles.Wavelet,'Enable','on','String','Spectral density');
        set(handles.onset,  'Enable','on');
        set(handles.dur,    'Enable','on');
        set(handles.model,  'Value',3,'Enable','off');
        set(handles.Rft,    'Enable','on');
        set(handles.Hz1,    'Enable','on');
        set(handles.Hz2,    'Enable','on');
        
    % induced responses
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
            'Input (C x u)'
            'Dipoles'
            'Save results as img'};
        try
            set(handles.Nmodes, 'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes, 'Value', 4);
        end
        
        set(handles.text20,     'String', 'modes');
        set(handles.model,      'Enable','off');
        if get(handles.Spatial, 'Value') > 2
            set(handles.Spatial,'Value', 2);
        end
        set(handles.Spatial,    'String',{'ECD','LFP'});
        set(handles.Wavelet,    'Enable','on','String','Wavelet transform');
        set(handles.Imaging,    'Enable','off' )
        set(handles.priors,     'Enable','off' )
        set(handles.onset,      'Enable','on');
        set(handles.dur,        'Enable','on');
        set(handles.Rft,        'Enable','on');
        set(handles.Hz1,        'Enable','on');
        set(handles.Hz2,        'Enable','on');
        
    case{'PHA'}
          Action = {
            'Sin(Data) - Region 1'
            'Sin(Data) - Region 2'
            'Sin(Data) - Region 3'
            'Sin(Data) - Region 4'
            'Coupling (As)'
            'Coupling (Bs)'};
       
        try
            set(handles.Nmodes, 'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes, 'Value', 1);
        end
        
        set(handles.text20,      'String', 'sub-trials');
        set(handles.model,       'Enable','off');
        if get(handles.Spatial,  'Value') > 2
            set(handles.Spatial, 'Value', 2);
        end
        set(handles.Spatial, 'String',{'ECD','LFP'});
        set(handles.Wavelet, 'Enable','on','String','Hilbert transform');
        set(handles.Imaging, 'Enable','off' )
        set(handles.onset,   'Enable','off');
        set(handles.dur,     'Enable','off');
        set(handles.Rft,     'Enable','off');
        set(handles.Hz1,     'Enable','off');
        set(handles.Hz2,     'Enable','off');
           
    % Cross-spectral density model (complex)
    %----------------------------------------------------------------------
    case {'NFM'}
        Action = {
              'spectral data',...
              'Coupling (A)',...
              'Coupling (B)',...
              'Coupling (C)',...
              'trial-specific effects',...
              'Input',...
              'Transfer functions',...
              'Cross-spectra (sources)',...
              'Cross-spectra (channels)',...
              'Coherence (sources)',...
              'Coherence (channels)',...
              'Covariance (sources)',...
              'Covariance (channels)',...
              'Dipoles'};
        try
            set(handles.Nmodes, 'Value', handles.DCM.options.Nmodes);
        catch
            set(handles.Nmodes, 'Value', 4);
        end
        
        set(handles.text20, 'String', 'modes');
        set(handles.model,  'Value',10,'Enable','off');
        set(handles.Spatial,'Value',1,'String',{'LFP'}); 
        set(handles.Wavelet,'Enable','on');
        set(handles.onset,  'Enable','off');
        set(handles.dur,    'Enable','off');
        set(handles.Rft,    'Enable','on');
        set(handles.Hz1,    'Enable','on');
        set(handles.Hz2,    'Enable','on');
        
    otherwise
        warndlg('unknown analysis type')
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
handles = reset_Callback(hObject, eventdata, handles);
 
handles.DCM = spm_dcm_erp_dipfit(handles.DCM, 1);

switch handles.DCM.options.analysis
    
    case{'CSD','NFM'}
        
        % cross-spectral density (if DCM.M.U (eigen-space) exists
        %------------------------------------------------------------------
        try
            handles.DCM = spm_dcm_csd_data(handles.DCM);
        end
 
        % and display
        %------------------------------------------------------------------
        spm_dcm_csd_results(handles.DCM,'spectral data');
        
 
    case{'IND'}
        
        % wavelet tranform
        %------------------------------------------------------------------
        handles.DCM = spm_dcm_ind_data(handles.DCM);
        
        % and display
        %------------------------------------------------------------------
        spm_dcm_ind_results(handles.DCM,'Wavelet');
        
    
    case{'TFM'}
        
        % wavelet tranform
        %------------------------------------------------------------------
        handles.DCM = spm_dcm_tfm_data(handles.DCM);
 
        % and display
        %------------------------------------------------------------------
        spm_dcm_tfm_results(handles.DCM,'induced and evoked responses');
        
        
     case{'PHA'}
        handles.DCM = spm_dcm_phase_data(handles.DCM);
 
end
guidata(hObject,handles);

% --- Executes on button press in priors.
function priors_Callback(hObject, eventdata, handles)
handles = reset_Callback(hObject, eventdata, handles);
spm_api_nmm(handles.DCM)

% --- Executes on button press in save_spec.
function save_spec_Callback(hObject, eventdata, handles)
spm_dcm_bmr;


% --- Executes on button press in reduce.
function reduce_Callback(hObject, eventdata, handles)
handles = reset_Callback(hObject, eventdata, handles);
spm_dcm_post_hoc(handles.DCM);





