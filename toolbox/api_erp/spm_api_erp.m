function varargout = spm_api_erp(varargin)
% SPM_API_ERP Application M-file for spm_api_erp.fig
%    FIG = SPM_API_ERP launch spm_api_erp GUI.
%    SPM_API_ERP('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 24-Oct-2005 11:46:40

if nargin == 0  % LAUNCH GUI

	fig     = openfig(mfilename,'reuse');
    Fgraph  = spm_figure('GetWin','Graphics');
	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    
    DCM.name = 'ERP';
    DCM.Y.xy = {};
    
    handles.DCM    = DCM;
    handles.X = [];
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


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


% --------------------------------------------------------------------
function varargout = swd_Callback(h, eventdata, handles, varargin)
try
    cd(get(handles.swd,'String'));
    set(handles.swd,'String',pwd)
end
handles.DCM.swd = pwd;


% --------------------------------------------------------------------
function varargout = cd_Callback(h, eventdata, handles, varargin)

try
    p = uigetdir(pwd, 'select SWD');
    
    if p ~= 0
        cd(p);
    else, return, end
    set(handles.swd, 'String', pwd)
end

handles.DCM.swd = pwd;


% --------------------------------------------------------------------
function varargout = name_Callback(h, eventdata, handles, varargin)
handles.DCM.name   = get(handles.name,'String');
set(handles.save,'enable','on')
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = Y_Callback(h, eventdata, handles, varargin)

DCM = handles.DCM;
% get y
%---------------------------------------------------------------------
[f,p] = uigetfile({'*.mat'}, 'please select ERP/ERF mat file (SPM format)');


if p ~=0
    f     = fullfile(p,f);
    D = spm_eeg_ldata(f);
else, return, end

% indices of EEG channel (excluding bad channels)
DCM.M.Ichannels = setdiff(D.channels.eeg, D.channels.Bad);
DCM.Y.Time = [-D.events.start:D.events.stop]*1000/D.Radc; % ms
DCM.Y.dt = 1000/D.Radc;

% if MEG, store grad struct in D.channels
try    
    handles.grad = D.channels.grad;
end

try
    DCM.Y.xy{1} = D.data(DCM.M.Ichannels, :, str2num(get(handles.Y1, 'String')))';
catch
    warndlg('First evoked response could not be read');
end

if ~isempty(deblank(get(handles.Y2, 'String')))
    try
        DCM.Y.xy{2} = D.data(DCM.M.Ichannels, :, str2num(get(handles.Y2, 'String')))';
    catch
        warndlg('Second evoked response could not be read');
    end
end

% save original data here
DCM.Y.xy_original = DCM.Y.xy;
DCM.Y.Time_original = DCM.Y.Time;

% update X
%---------------------------------------------------------------------
try
    m = length(DCM.Y.xy);
    spm_dcm_erp_results(DCM, 'Response');
    
    % default design
    %-----------------------------------------------------------------
    X(m) = 1;
    DCM.U.X = X';
    DCM.U.name = {'input 1'};
catch
    warndlg('these data are not compatible')
end

% save DCM
handles.DCM = DCM;
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = h_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = Sname_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = L_Callback(h, eventdata, handles, varargin)
[f,p] = uigetfile({'*.mat';'*.txt';'*.dat'},'please select lead field matrix');

if p ~= 0
    f     = fullfile(p,f);
    try
        L = load(f,'-ascii');
    end
    try
    L = load(f,'-mat');
    while ~isnumeric(L)
        s = fieldnames(L);
        L = getfield(L,s{1});
    end
    end
else, return, end

% check the number of channels is consistent with the lead field
%---------------------------------------------------------------------
try
    if size(L, 1) ~= size(handles.DCM.E, 2)
        warndlg(sprintf('lead field (%d) and channel numbers (%d) are inconsistent', size(L, 1),  size(handles.DCM.E, 2)))
        return
    end
end

% and save in DCM.L
%---------------------------------------------------------------------
try
    handles.DCM.L  = L;
    guidata(h,handles);
    Sname_Callback(h, eventdata, handles);
    reset_Callback(h, eventdata, handles);
catch
    warndlg('these data are not compatible')
end

% --------------------------------------------------------------------
function varargout = connections_Callback(h, eventdata, handles, varargin)
n = size(handles.DCM.L,2);    % number of sources
m = size(handles.DCM.U.X,2);  % number of inputs

% put onset here (later: enable user to change this using GUI)
handles.DCM.M.onset = 15;

% reset connection buttons
%---------------------------------------------------------------------
reset_Callback(h, eventdata, handles, varargin)

% connection buttons
%---------------------------------------------------------------------
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
%-------------------------------------------------------------------------
constr         = {'A forward' 'A backward' 'A lateral' 'C input'};
nsx   = (n + 1)*sx;
nsy   = 2*sy;
for k = 1:4
    x          = x0 + (k - 1)*nsx;
    y          = y0 - 3*sy;
    str        = constr{k};
    S(k)       = uicontrol(handles.SPM,...
        'Units','Normalized',...
        'Position',[x y nsx nsy],...
        'HorizontalAlignment','left',...
        'Style','text',...
        'String',str, 'BackgroundColor', [1 1 1]);
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
                'String',str, 'BackgroundColor', [1 1 1]);
end
handles.S = S;
handles.A = A;
handles.B = B;
handles.C = C;
set(handles.connections,'Enable','off')
set(handles.reset,      'Enable','on' )
set(handles.estimate, 'Enable', 'on');


guidata(h,handles)


% --------------------------------------------------------------------
function handles = save_Callback(h, eventdata, handles, varargin)

DCM = handles.DCM;

[f,p] = uiputfile(['DCM' DCM.name '.mat'],'save DCM');

if p ~= 0
    % store all the options for loading again
    DCM.options.Y1 = str2num(get(handles.Y1, 'String'));
    DCM.options.Y2 = str2num(get(handles.Y2, 'String'));
    DCM.options.T1 = str2num(get(handles.T1, 'String'));
    DCM.options.T2 = str2num(get(handles.T2, 'String'));
    DCM.options.Spatial_type = get(handles.Spatial_type, 'Value');
    DCM.M.Spatial_type = DCM.options.Spatial_type; % store for use as switch in spm_gx_erp
    DCM.options.Nmodes = get(handles.Nmodes, 'Value');
    DCM.options.h = get(handles.h, 'Value')-1;

    handles.DCM = DCM;

    save(fullfile(p, f), 'DCM')
    set(handles.estimate, 'Enable', 'on')
    set(handles.results, 'Enable', 'off')

    try
        handles.DCM.F;
        set(handles.estimate,'String','re-estimate')
    end
end
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = load_Callback(h, eventdata, handles, varargin)

[f,p]       = uigetfile('DCM*.mat','please select DCM');

if p ~= 0
    DCM         = load(fullfile(p,f), '-mat');
    DCM         = DCM.DCM;
    handles.DCM = DCM;
    guidata(h, handles)
else, return, end
% enter options from saved options and try to execute data_ok and spatial_ok

% put back original data and restore settings
try
    DCM.Y.xy = DCM.Y.xy_original;
    DCM.Y.Time = DCM.Y.Time_original;
end

try, set(handles.Y1, 'String', num2str(DCM.options.Y1)); end
try, set(handles.Y2, 'String', num2str(DCM.options.Y2)); end
try, set(handles.T1, 'String', num2str(DCM.options.T1)); end
try, set(handles.T2, 'String', num2str(DCM.options.T2)); end
try, set(handles.Spatial_type, 'Value', DCM.options.Spatial_type); end
try, set(handles.Nmodes, 'Value', DCM.options.Nmodes); end
try, set(handles.h, 'Value', DCM.options.h+1); end

handles.DCM = DCM;
guidata(h, handles);

% data_ok
if isfield(DCM.options, 'data_ok')
    try
        handles = data_ok_Callback(h, eventdata, handles);
    end
end

try, set(handles.name, 'String', DCM.name); end
try, set(handles.Sname, 'String', strvcat(DCM.Sname)); end

if DCM.options.Spatial_type ~= 3
    % transform locations back to MNI space first
    M = [[0 -1 0 0]; [1 0 0 -20]; [0 0 1 -10]; [0 0 0 1]]; % transformation matrix from dip to MNI
    tmp = M*[DCM.M.dipfit.L.pos; ones(1, size(DCM.M.dipfit.L.pos, 2))];

    try, set(handles.Slocation, 'String', num2str(tmp(1:3, :)')); end

end

% saving DCM
guidata(h, handles);

% spatial_ok
if isfield(DCM.options, 'spatial_ok')
    try
        handles = spatial_ok_Callback(h, eventdata, handles);
    end

    % saving DCM
    guidata(h, handles);

    try
        spm_api_erp('connections_Callback', gcbo, [], handles); 
        set(handles.estimate, 'Enable', 'on');
    end
end

% make model specification DCM.M global variable
clear M; global M
M = handles.DCM.M;

try
    handles.DCM.F;
    set(handles.results, 'Enable','on');
    set(handles.estimate,'String','Re-estimate','Foregroundcolor', [0 0 0]);
catch
    % if something has been estimated earlier
    set(handles.results, 'Enable','off');
    set(handles.estimate,'String','Estimate','Foregroundcolor', [0 0 0]);
end

% --------------------------------------------------------------------
function varargout = estimate_Callback(h, eventdata, handles, varargin)
set(handles.estimate,'String','Estimating','Foregroundcolor',[1 0 0])

handles = save_Callback(h, eventdata, handles);

handles.DCM = spm_dcm_erp(handles.DCM);


save_Callback(h, eventdata, handles, varargin)
set(handles.results,    'Enable','on' )
set(handles.save,       'Enable','off')
set(handles.Y,          'Enable','off')
set(handles.L,          'Enable','off')
set(handles.name,       'Enable','off')
set(handles.swd,        'Enable','off')
set(handles.reset,      'Enable','off')

% reset and disable connections for contrast specification
%-------------------------------------------------------------------
% for i = 1:length(handles.A)
%     for j = 1:length(handles.A{i})
%         for k = 1:length(handles.A{i})
%             if ~handles.DCM.A{i}(j,k)
%                 set(handles.A{i}(j,k),'enable','off');
%             end
%             set(handles.A{i}(j,k),'Value',0);
%             handles.DCM.A{i}(j,k)       = 0;
%         end
%     end
% end
% for i = 1:length(handles.B)
%     for j = 1:length(handles.B{i})
%         for k = 1:length(handles.B{i})
%             if ~handles.DCM.B{i}(j,k)
%                 set(handles.B{i}(j,k),'enable','off');
%             end
%             set(handles.B{i}(j,k),'Value',0);
%             handles.DCM.B{i}(j,k)       = 0;
%         end
%     end
% end
% for i = 1:length(handles.C)
%     set(handles.C(i),'Value',0);
%     set(handles.C(i),'enable','off');
%     handles.DCM.C(i)       = 0;
% end
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = contrast_Callback(h, eventdata, handles, varargin)
% This callback has no link to it. I'll put it into results later.
% get contrast
%---------------------------------------------------------------------
[n m] = size(handles.DCM.C);
c     = handles.DCM.Ep*0;
c     = spm_erp_pack(c,m,n);
c.A   = handles.DCM.A;
c.B   = handles.DCM.B;
c.C   = handles.DCM.C;
c     = spm_vec(c);
c     = c/sum(c);

warning off
Ec    = full(c'*handles.DCM.Ep);
Cc    = full(c'*handles.DCM.Cp*c);
x     = Ec + [-32:32]*sqrt(Cc)*4/32;
pdf   =     spm_Npdf(x,Ec,Cc);
PP    = 1 - spm_Ncdf(0,abs(Ec),Cc);
warning on

% plot conditional density
%---------------------------------------------------------------------
figure(handles.Fgraph)
clf
subplot(2,1,1);
plot(x,pdf);
title({'Posterior density of contrast',...
        sprintf('P(contrast > 0) = %.1f%s',PP*100,'%')},...
       'FontSize',12)
xlabel('contrast')
ylabel('probability density')
axis square, grid on

% --------------------------------------------------------------------
function varargout = results_Callback(h, eventdata, handles, varargin)
Action  = get(handles.results,'String');
Action  = Action{get(handles.results,'Value')};
spm_dcm_erp_results(handles.DCM,Action);

% remove existing buttons
%---------------------------------------------------------------------
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


% --- Executes during object creation, after setting all properties.
function Slocation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Slocation_Callback(hObject, eventdata, handles)
% hObject    handle to Slocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Slocation as text
%        str2double(get(hObject,'String')) returns contents of Slocation as a double


% --- Executes during object creation, after setting all properties.
function Spatial_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Spatial_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Spatial_type.
function Spatial_type_Callback(hObject, eventdata, handles)
% hObject    handle to Spatial_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Spatial_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Spatial_type

if get(handles.Spatial_type, 'Value') == 3
    % if fixed lead field
    set(handles.L, 'Enable', 'on');
    set(handles.sensorfile, 'Enable', 'off');
    set(handles.plot_dipoles, 'Enable', 'off');
    set(handles.Slocation, 'Enable', 'off');
elseif get(handles.Spatial_type, 'Value') == 2
    % if MEG ECD
    set(handles.L, 'Enable', 'off');
    set(handles.sensorfile, 'Enable', 'off');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation, 'Enable', 'on');

else
    % if EEG ECD
    set(handles.L, 'Enable', 'off');
    set(handles.sensorfile, 'Enable', 'on');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation, 'Enable', 'on');
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Y1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function Y1_Callback(hObject, eventdata, handles)
% hObject    handle to Y1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Y2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Y2_Callback(hObject, eventdata, handles)
% hObject    handle to Y2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Y2 as text
%        str2double(get(hObject,'String')) returns contents of Y2 as a double


% --- Executes during object creation, after setting all properties.
function Nmodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nmodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Nmodes.
function Nmodes_Callback(hObject, eventdata, handles)
% hObject    handle to Nmodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Nmodes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Nmodes


% --- Executes during object creation, after setting all properties.
function projection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in projection.
function projection_Callback(hObject, eventdata, handles)
% hObject    handle to projection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns projection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from projection

if get(handles.projection, 'Value') == 1
    % if SVD
    set(handles.Nmodes, 'Enable', 'on')
else
    set(handles.Nmodes, 'Enable', 'off')
end
guidata(hObject,handles);


% --- Executes on button press in sensorfile.
function sensorfile_Callback(h, eventdata, handles)
% hObject    handle to sensorfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Spatial_type = get(handles.Spatial_type, 'Value');
DCM = handles.DCM;

if Spatial_type == 1
    % EEG
    [f,p] = uigetfile({'*.pol';'*.mat'}, 'please select sensor location file');
    if p == 0, return, end
    f     = fullfile(p,f);
    
    DCM.M.dipfit.sensorfile = f;
        
    handles.DCM = DCM;
end

guidata(h, handles);


% --- Executes on button press in spatial_ok.
function handles = spatial_ok_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% spatial model
%------------------------

DCM = handles.DCM;

% retrieve non-blank names
tmp = deblank(get(handles.Sname, 'String'));

if size(tmp, 1) == 0
    errordlg('Please specify the source names'); return
end

k = 1;
for i = 1:size(tmp, 1)
    if ~isempty(deblank(tmp(i, :)))
        Sname{k} = deblank(tmp(i, :));
        k = k + 1;
    end
end

Nareas = length(Sname);
Nmodes = size(DCM.E, 1);

DCM.Sname = Sname;

Spatial_type = get(handles.Spatial_type, 'Value');
DCM.options.Spatial_type = Spatial_type;
DCM.M.Spatial_type = Spatial_type; % store for use as switch in spm_gx_erp


if Spatial_type ~= 3
    % user chose an ECD model
    % read location coordinates
    Slocations = zeros(Nareas, 3);
    tmp = get(handles.Slocation, 'String');

    if ~isempty(tmp) & size(tmp, 1) == Nareas
        for i = 1:Nareas
            tmp2 = str2num(tmp(i, :));
            if length(tmp2) ~= 3
                errordlg(sprintf('coordinates of area %d not valid', i)); return;
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
else
    
    % fixed lead field

    set(handles.Slocation, 'String', '');
    
    try
        DCM.L;
    catch
        errordlg('Please specify lead field'); return;
    end
    
    if size(DCM.L, 2) ~= Nareas
        errordlg(sprintf('fixed lead field (%d) and nr of components (%d) don''t correspond.', size(DCM.L, 2), Nareas)); return;
    end

    DCM.L = DCM.E*DCM.L;
    DCM.M.L = DCM.L;
    
end

if Spatial_type == 1
    try
        f = DCM.M.dipfit.sensorfile;
    catch
        errordlg('Please specify sensor location file'); return;
    end
    
elseif Spatial_type == 2
    DCM.M.grad = handles.grad;
    % check whether I need to normalise anything (like in the EEG case)
end

if Spatial_type == 1 || Spatial_type == 2
    
    % prepare forward model
    DCM = spm_dcm_erp_prepareSpatial(DCM);
    
    % specify prior distributions for ECD model parameters
    DCM.L = zeros(Nmodes, Nareas);

    % prior distributions of spatial parameters
    % they are kept here (and not in spm_erp_priors, because I want
    % them to be changeable by (power-)user
    % also: transform to eeglab coordinates
    % transformation matrix from dip to MNI
    M = [[0 -1 0 0]; [1 0 0 -20]; [0 0 1 -10]; [0 0 0 1]];

    Slocations = inv(M)*[Slocations ones(Nareas, 1)]';
    DCM.M.dipfit.L.pos = Slocations(1:3,:);

    % zero prior mean on dipole orientation
    DCM.M.dipfit.L.mom = zeros(3,Nareas);

    % prior mean of one on K
    DCM.M.dipfit.L.K = ones(Nareas, 1);

    % tight priors on location, broad priors on moments
    DCM.M.dipfit.L.Vpos = [8*ones(3, Nareas)];
    DCM.M.dipfit.L.Vmom = [8*ones(3, Nareas)];

    % prior precision on K is a bit redundant because there is already broad
    % prior on moments
    DCM.M.dipfit.L.VK = 1*ones(size(DCM.M.dipfit.L.K));

    DCM.M.E = DCM.E;
    DCM.M.L = zeros(Nmodes, Nareas);
    DCM.M.Lpos = NaN*DCM.M.dipfit.L.pos; % NaN forces spm_gx_erp to compute initial lead field
    DCM.M.Lmom = NaN*DCM.M.dipfit.L.mom; %
end

DCM.options.spatial_ok = 1;

handles.DCM = DCM;
set(handles.Spatial_type, 'Enable', 'off');
set(handles.sensorfile, 'Enable', 'off');
set(handles.L, 'Enable', 'off');
set(handles.spatial_ok, 'Enable', 'off');
set(handles.Sname, 'Enable', 'off');
set(handles.Slocation, 'Enable', 'off');
set(handles.spatial_back, 'Enable', 'off');

set(handles.connections, 'Enable', 'on');
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
function handles = data_ok_Callback(hObject, eventdata, handles)
% hObject    handle to data_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Project data

DCM = handles.DCM;

% check all entries
Nresponses = length(DCM.Y.xy);
if Nresponses == 0
    errordlg('Please choose some data'); return;
end

T1 = str2num(get(handles.T1, 'String'));
T2 = str2num(get(handles.T2, 'String'));

if T1 <= 0
    errordlg('Specify positive first peri-stimulus time'); return;
end

if T2 <= T1
    errorldg('Second must be greater than first peri-stimulus time'); return;
end

[m, T1] = min(abs(DCM.Y.Time-T1));
[m, T2] = min(abs(DCM.Y.Time-T2));

DCM.Y.Time = DCM.Y.Time(T1:T2);

for i = 1:Nresponses
    DCM.Y.xy{i} = DCM.Y.xy{i}(T1:T2, :);
end

Nchannels = size(DCM.Y.xy{1}, 2);

% nr of modes
projection = get(handles.projection, 'Value');
if projection == 1
    % svd
    Nmodes = get(handles.Nmodes, 'Value');
else
    % no projection
    Nmodes = Nchannels;
end

% project data
P.projection = projection;
P.Nmodes = Nmodes;
DCM = spm_dcm_eeg_selectdata(DCM, P);

spm_dcm_erp_results(DCM, 'Response');

s = get(handles.h, 'String');
DCM.Y.h = str2num(s{get(handles.h, 'Value')});

DCM.options.data_ok = 1;

handles.DCM = DCM;

% enable next stage, disable data specification
set(handles.Y, 'Enable', 'off');
set(handles.Y1, 'Enable', 'off');
set(handles.Y2, 'Enable', 'off');
set(handles.T1, 'Enable', 'off');
set(handles.T2, 'Enable', 'off');
set(handles.projection, 'Enable', 'off');
set(handles.Nmodes, 'Enable', 'off');
set(handles.h, 'Enable', 'off');
set(handles.data_ok, 'Enable', 'off');
set(handles.Spatial_type, 'Enable', 'on');
set(handles.plot_dipoles, 'Enable', 'on');

if get(handles.Spatial_type, 'Value') == 1
    % ECD EEG
    set(handles.sensorfile, 'Enable', 'on');
    set(handles.L, 'Enable', 'off');
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation, 'Enable', 'on');
elseif get(handles.Spatial_type, 'Value') == 2
    % ECD MEG
    set(handles.sensorfile, 'Enable', 'off');
    set(handles.L, 'Enable', 'off');    
    set(handles.plot_dipoles, 'Enable', 'on');
    set(handles.Slocation, 'Enable', 'on');    
else
    % fixed
    set(handles.sensorfile, 'Enable', 'off');
    set(handles.L, 'Enable', 'on');    
    set(handles.plot_dipoles, 'Enable', 'off');
    set(handles.Slocation, 'Enable', 'off');
end

set(handles.spatial_ok, 'Enable', 'on');
set(handles.Sname, 'Enable', 'on');
set(handles.Slocation, 'Enable', 'on');
set(handles.spatial_back, 'Enable', 'on');


guidata(hObject, handles);


% --- Executes on button press in spatial_back.
function spatial_back_Callback(hObject, eventdata, handles)
% hObject    handle to spatial_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% restore original data (before projection)
handles.DCM.Y.xy = handles.DCM.Y.xy_original;
handles.DCM.Y.Time = handles.DCM.Y.Time_original;

set(handles.Y, 'Enable', 'on');
set(handles.Y1, 'Enable', 'on');
set(handles.Y2, 'Enable', 'on');
set(handles.T1, 'Enable', 'on');
set(handles.T2, 'Enable', 'on');
set(handles.projection, 'Enable', 'on');
set(handles.Nmodes, 'Enable', 'on');
set(handles.h, 'Enable', 'on');
set(handles.data_ok, 'Enable', 'on');
set(handles.Spatial_type, 'Enable', 'off');
set(handles.sensorfile, 'Enable', 'off');
set(handles.L, 'Enable', 'off');
set(handles.spatial_ok, 'Enable', 'off');
set(handles.Sname, 'Enable', 'off');
set(handles.Slocation, 'Enable', 'off');
set(handles.spatial_back, 'Enable', 'off');
set(handles.plot_dipoles, 'Enable', 'off');

guidata(hObject, handles);



function T1_Callback(hObject, eventdata, handles)
% hObject    handle to T1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T1 as text
%        str2double(get(hObject,'String')) returns contents of T1 as a double


% --- Executes during object creation, after setting all properties.
function T1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T2_Callback(hObject, eventdata, handles)
% hObject    handle to T2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T2 as text
%        str2double(get(hObject,'String')) returns contents of T2 as a double


% --- Executes during object creation, after setting all properties.
function T2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in connectivity_back.
function connectivity_back_Callback(hObject, eventdata, handles)
% hObject    handle to connectivity_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Spatial_type, 'Enable', 'on');

if get(handles.Spatial_type, 'Value') == 1
    set(handles.sensorfile, 'Enable', 'on');
    set(handles.L, 'Enable', 'off');
else
    set(handles.sensorfile, 'Enable', 'off');
    set(handles.L, 'Enable', 'on');    
end

set(handles.spatial_ok, 'Enable', 'on');
set(handles.Sname, 'Enable', 'on');
set(handles.Slocation, 'Enable', 'on');
set(handles.spatial_back, 'Enable', 'on');

set(handles.connections, 'Enable', 'off');
set(handles.connectivity_back, 'Enable', 'off');
set(handles.reset, 'Enable', 'off');

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

% --- Executes on button press in plot_dipoles.
function plot_dipoles_Callback(hObject, eventdata, handles)
% hObject    handle to plot_dipoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% can be called before and after calling spatial_ok_Callback.

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

Nlocations = size(Slocations, 2);

sdip.n_seeds = 1;
sdip.n_dip = Nlocations;
sdip.Mtb = 1;
sdip.j{1} = zeros(3*Nlocations, 1);
sdip.loc{1} = Slocations;

spm_eeg_inv_ecd_DrawDip('Init', sdip)
