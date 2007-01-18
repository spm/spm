function varargout = spm_eeg_inv_visu3D_api(varargin)
% SPM_EEG_INV_VISU3D_API M-file for spm_eeg_inv_visu3D_api.fig
% - FIG = SPM_EEG_INV_VISU3D_API launch spm_eeg_inv_visu3D_api GUI.
% - D   = SPM_EEG_INV_VISU3D_API(D) open with D
% - SPM_EEG_INV_VISU3D_API(filename) where filename is the eeg/meg .mat file
% - SPM_EEG_INV_VISU3D_API('callback_name', ...) invoke the named callback.
%
% Last Modified by GUIDE v2.5 17-Jan-2007 13:15:47
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
% Jeremie Mattout
% $Id: $

% INITIALISATION CODE
%--------------------------------------------------------------------------
if nargin == 0 || nargin == 1  % LAUNCH GUI

    fig         = openfig(mfilename,'new');
    handles     = guihandles(fig);
    handles.fig = fig;
    
    try
       set(handles.DataFile,'String',varargin{1}.fname)
       set(handles.Contrast,'String',num2str(varargin{1}.val));
       handles.D = varargin{1};
    end
    
    try  
        cla(handles.sensors_axes)
        cla(handles.sources_axes)
        cla(handles.pred_axes)
    end
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    set(fig,'Menubar','figure');
    guidata(fig,handles);
    
    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1})

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before spm_eeg_inv_visu3D_api is made visible.
function spm_eeg_inv_visu3D_api_OpeningFcn(hObject, eventdata, handles)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    D = handles.D;
catch
    D = spm_eeg_ldata(spm_select(1, '.mat', 'Select EEG/MEG mat file'));
end

if ~isfield(D,'inv')
    error(sprintf('Please specify and invert a foward model\n'));
end

set(handles.DataFile,'String',D.fname);
set(handles.Contrast,'String',sprintf('contrast: %i',D.val));
set(handles.fig,'name',['Source visualisation -' D.fname])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET RESULTS (default: current or last analysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    val = D.val;
catch
    val = 1;
end
if val < length(D.inv)
    set(handles.next,    'Enable','on','Value',1);
else
    set(handles.next,    'Enable','off','Value',1);
end
if val > 1
    set(handles.previous,'Enable','on','Value',1);
else
    set(handles.previous,'Enable','off','Value',1);
end
if strcmp(D.inv{val}.method,'ECD')
    warndlg('Please create an imaging solution');
    guidata(hObject,handles);
    return
end
set(handles.LogEv,'String',num2str(D.inv{val}.inverse.F));
set(handles.LogEv,'Enable','inactive');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBSERVED ACTIVITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start with response
%--------------------------------------------------------------------------
try
    
    % Load Gain or Lead field matrix
    %----------------------------------------------------------------------
    dimT              = 256;
    dimS              = D.inv{val}.inverse.Nd;
    Cs                = setdiff(D.channels.eeg, D.channels.Bad);
    Is                = D.inv{val}.inverse.Is;
    Ts                = ceil([1:dimT]*length(D.inv{val}.inverse.pst)/dimT);
    L                 = D.inv{val}.inverse.L;
    U                 = D.inv{val}.inverse.U;
    R                 = D.inv{val}.inverse.R;
    T                 = D.inv{val}.inverse.T;
    
    % source data
    %----------------------------------------------------------------------
    set(handles.Activity,'Value',1);
    J                 = sparse(dimS,dimT);
    J(Is,:)           = D.inv{val}.inverse.J*T(Ts,:)';
    handles.dimT      = dimT;
    handles.dimS      = dimS;
    handles.pst       = D.inv{val}.inverse.pst(Ts);
    handles.srcs_data = J;
    
    % sensor data
    %----------------------------------------------------------------------
    Y     = sparse(0);
    c     = find(D.events.code == D.events.types(D.inv{val}.inverse.con));
    for i = 1:length(c)
        y = U'*R*(squeeze(D.data(Cs,:,c(i)))*T);
        Y = Y + y;
    end
    handles.sens_data = U*Y*T(Ts,:)';
    handles.pred_data = U*L*J(Is,:);

catch
    warndlg({'Please invert your model';'inverse.J not found'});
    return
end
    
% case 'windowed response' or contrast'
%--------------------------------------------------------------------------
try
    JW                    = sparse(dimS,1);
    GW                    = sparse(dimS,1);
    JW(Is,:)              = D.inv{val}.contrast.JW;
    GW(Is,:)              = D.inv{val}.contrast.GW;
    handles.woi           = D.inv{val}.contrast.woi;
    handles.fboi          = D.inv{val}.contrast.fboi;
    handles.W             = D.inv{val}.contrast.W(Ts,:);
    handles.srcs_data_w   = JW;
    handles.sens_data_w   = handles.sens_data*handles.W(:,1);
    handles.pred_data_w   = handles.pred_data*handles.W(:,1);
    handles.srcs_data_ev  = GW;
    handles.sens_data_ev  = sum((handles.sens_data*handles.W).^2,2);
    handles.pred_data_ev  = sum((handles.pred_data*handles.W).^2,2);
    set(handles.Activity,'enable','on');
catch
    set(handles.Activity,'enable','off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD CORTICAL MESH (default: Individual)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    vert = D.inv{val}.mesh.tess_mni.vert;
    face = D.inv{val}.mesh.tess_mni.face;
    set(handles.Template,  'Value',1);
    set(handles.Individual,'Value',0);
catch
    try
        vert = D.inv{val}.mesh.tess_ctx.vert;
        face = D.inv{val}.mesh.tess_ctx.face;
        set(handles.Template,  'Value',0);
        set(handles.Individual,'Value',1);
    catch
        warndlg('There is no mesh associated with these data');
        return
    end
end

handles.vert  = vert;
handles.face  = face;
handles.grayc = sqrt(sum((vert.^2)')); handles.grayc = handles.grayc'/max(handles.grayc);
clear vert face

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLIDER INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.slider_transparency,'Min',0,'Max',1,'Value',1,'sliderstep',[0.01 0.05]);
set(handles.slider_srcs_up,     'Min',0,'Max',1,'Value',0,'sliderstep',[0.01 0.05]);
set(handles.slider_srcs_down,   'Min',0,'Max',1,'Value',1,'sliderstep',[0.01 0.05]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL SOURCE LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.sources_axes);
cla; axis off

set(handles.slider_time,  'Enable','on');
set(handles.time_bin,     'Enable','on');
set(handles.slider_time,  'Value',1);
set(handles.time_bin,     'String',num2str(fix(handles.pst(1))));
set(handles.slider_time,  'Min',1,'Max',handles.dimT,'sliderstep',[1/(handles.dimT-1) 2/(handles.dimT-1)]);
set(handles.checkbox_absv,'Enable','on','Value',1);
set(handles.checkbox_norm,'Enable','on','Value',0);

handles.srcs_disp = full(handles.srcs_data(:,1));
handles.fig1 = patch('vertices',handles.vert,'faces',handles.face,'FaceVertexCData',handles.srcs_disp);

% display
%--------------------------------------------------------------------------
set(handles.fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
shading interp
lighting gouraud
camlight
cameramenu
lightangle(0,270);lightangle(270,0),lightangle(0,0),lightangle(90,0);
material([.1 .1 .4 .5 .4]);
view(140,15);
axis image
UpDate_Display_SRCS(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD SENSOR FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    load(fullfile(spm('dir'),'EEGtemplates',D.channels.ctf));
catch
    [f,p] = uigetfile({'*.mat'}, 'EEGtemplate file (Cpos)');
    load(fullfile(p,f))
end
try
    Cpos = Cpos(:,D.channels.order(Cs));
catch
    Cpos = Cpos(:,D.channels.order(D.channels.eeg));
end
xp       = Cpos(1,:)';
yp       = Cpos(2,:)';
x        = min(xp):0.01:max(yp);
y        = min(yp):0.01:max(yp);
[xm,ym]  = meshgrid(x,y);
clear Cnames Rxy Cpos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL SENSOR LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.sensors_axes);
cla; axis off

disp     = full(handles.sens_data(:,1));
imagesc(y,x,griddata(xp,yp,disp,xm,ym));
axis image xy off
handles.sens_coord_x   = x;
handles.sens_coord_y   = y;
handles.sens_coord     = [xp yp];
handles.sens_coord2D_X = xm;
handles.sens_coord2D_Y = ym;
hold on
handles.sensor_loc = plot(handles.sens_coord(:,1),handles.sens_coord(:,2),'o','MarkerFaceColor',[1 1 1]/2,'MarkerSize',6);
set(handles.checkbox_sensloc,'Value',1);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL SENSOR LEVEL DISPLAY - PREDICTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.pred_axes); cla;
disp      = full(handles.pred_data(:,1));
imagesc(griddata(xp,yp,disp,xm,ym));
axis image xy off

guidata(hObject,handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE SOURCE LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpDate_Display_SRCS(hObject,handles)
axes(handles.sources_axes);

% Adjust the threshold
%--------------------------------------------------------------------------
Set_colormap(hObject, [], handles);

if isfield(handles,'fig1')
    ActToDisp = get(handles.Activity,'Value');
    A         = get(handles.checkbox_absv,'Value');
    N         = get(handles.checkbox_norm,'Value');
    switch ActToDisp
        
        % case 1: response (J)
        %------------------------------------------------------------------
        case 1
            TS   = fix(get(handles.slider_time,'Value'));
            handles.srcs_disp = handles.srcs_data(:,TS);
            Nmax = max(abs(handles.srcs_data(:)));
            if A == 0 & N == 0
                handles.Vmin = min(handles.srcs_disp);
                handles.Vmax = max(handles.srcs_disp);
            elseif A == 0 & N == 1
                handles.Vmin = -Nmax;
                handles.Vmax =  Nmax;
                handles.srcs_disp = handles.srcs_disp;
            elseif A == 1 & N == 0
                handles.Vmin = 0;
                handles.Vmax = max(abs(handles.srcs_disp));
                handles.srcs_disp = abs(handles.srcs_disp);
            else % A == 1 & N == 1
                handles.Vmin = 0;
                handles.Vmax = Nmax;
                handles.srcs_disp = abs(handles.srcs_disp);
            end

            
        % case 2: Windowed response (JW)
        %------------------------------------------------------------------
        case 2
                handles.Vmin = min(handles.srcs_data_w);
                handles.Vmax = max(handles.srcs_data_w);
                handles.srcs_disp = handles.srcs_data_w;

        % case 3: Evoked power  (JWWJ)
        %------------------------------------------------------------------
        case 3
                handles.Vmin = min(handles.srcs_data_ev);
                handles.Vmax = max(handles.srcs_data_ev);
                handles.srcs_disp = handles.srcs_data_ev;

        % case 4: Induced power  (JWWJ)
        %------------------------------------------------------------------
        case 4
                handles.Vmin = min(handles.srcs_data_ind);
                handles.Vmax = max(handles.srcs_data_ind);
                handles.srcs_disp = handles.srcs_data_ind;
    end
    set(handles.fig1,'FaceVertexCData',full(handles.srcs_disp));
    set(handles.sources_axes,'CLim',[handles.Vmin handles.Vmax]);
    set(handles.sources_axes,'CLimMode','manual');
    colorbar    
end
drawnow
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE SENSOR LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpDate_Display_SENS(hObject,handles)
TypOfDisp = get(handles.sens_display,'Value');
ActToDisp = get(handles.Activity,'Value');

% topography
%--------------------------------------------------------------------------
if TypOfDisp == 1

    % responses at one pst
    %----------------------------------------------------------------------
    if ActToDisp == 1

        TS = fix(get(handles.slider_time,'Value'));
        sens_disp = handles.sens_data(:,TS);
        pred_disp = handles.pred_data(:,TS);
        
    % contrast
    %----------------------------------------------------------------------
    elseif ActToDisp == 2
        
        sens_disp = handles.sens_data_w;
        pred_disp = handles.pred_data_w;

    % power
    %----------------------------------------------------------------------
    elseif ActToDisp == 3
        sens_disp = handles.sens_data_ev;
        pred_disp = handles.pred_data_ev;
    end
    
    axes(handles.sensors_axes); cla
    disp = griddata(handles.sens_coord(:,1),handles.sens_coord(:,2),full(sens_disp),handles.sens_coord2D_X,handles.sens_coord2D_Y);
    imagesc(handles.sens_coord_y,handles.sens_coord_x,disp);
    axis image xy off

    axes(handles.pred_axes); cla
    disp = griddata(handles.sens_coord(:,1),handles.sens_coord(:,2),full(pred_disp),handles.sens_coord2D_X,handles.sens_coord2D_Y);
    imagesc(handles.sens_coord_y,handles.sens_coord_x,disp);
    axis image xy off

    % add sensor locations
    %----------------------------------------------------------------------
    if get(handles.checkbox_sensloc,'Value') == 1
        axes(handles.sensors_axes)
        hold on
        handles.sensor_loc = plot(handles.sens_coord(:,1),handles.sens_coord(:,2),'o','MarkerFaceColor',[1 1 1]/2,'MarkerSize',6);
        hold off
    end
    guidata(hObject,handles);

% time series
%--------------------------------------------------------------------------
elseif TypOfDisp == 2
    axes(handles.sensors_axes)
    daspect('auto')
    handles.fig2 = ...
        plot(handles.pst,handles.sens_data,'b-.',handles.pst,handles.pred_data,'r:');
    if ActToDisp > 1
        hold on
        Scal = norm(handles.sens_data,1)/norm(handles.W,1);
        plot(handles.pst,handles.W*Scal,'k')
        hold off
    end
    axis on tight;
    axes(handles.pred_axes); cla, axis off
end
drawnow
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataFile_Callback(hObject, eventdata, handles)
S     = get(handles.DataFile,'String');
try
    D = spm_eeg_ldata(S);
catch
    LoadData_Callback(hObject, eventdata, handles);
end

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, eventdata, handles)
S         = spm_select(1, '.mat', 'Select EEG/MEG mat file');
D         = spm_eeg_ldata(S);
handles.D = D;
spm_eeg_inv_visu3D_api_OpeningFcn(hObject, [], handles);

% --- Executes on button press in OpenData.
function OpenData_Callback(hObject, eventdata, handles)
spm_eeg_inv_visu3D_api_OpeningFcn(hObject, [], handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTIVITY TO DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in Activity.
function Activity_Callback(hObject, eventdata, handles)
ActToDisp = get(handles.Activity,'Value');
if ActToDisp == 1
    set(handles.checkbox_absv,   'Enable','on');
    set(handles.checkbox_norm,   'Enable','on');
    set(handles.slider_time,     'Enable','on');
    set(handles.time_bin,        'Enable','on');
else
    set(handles.checkbox_norm,   'Enable','off');
    set(handles.slider_time,     'Enable','off');
    set(handles.time_bin,        'Enable','off');
end
if ActToDisp == 2
    set(handles.checkbox_absv,   'Enable','off');
end

% update displays
%--------------------------------------------------------------------------
UpDate_Display_SRCS(hObject,handles);
UpDate_Display_SENS(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWITCH FROM TEMPLATE MESH TO INDIVIDUAL MESH AND BACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Individual_Callback(hObject, eventdata, handles)
set(handles.Template,'Value',0);
try
    handles.vert = handles.D.inv{handles.D.val}.mesh.tess_ctx.vert;
    set(handles.Template,  'Value',0);
    set(handles.Individual,'Value',1);
end
handles.grayc = sqrt(sum((handles.vert.^2)')); handles.grayc = handles.grayc'/max(handles.grayc);
set(handles.fig1,'vertices',handles.vert,'faces',handles.face);
UpDate_Display_SRCS(hObject,handles);
axes(handles.sources_axes);
axis image;
guidata(hObject,handles);

%--------------------------------------------------------------------------
function Template_Callback(hObject, eventdata, handles)
set(handles.Individual,'Value',0);
try
    handles.vert = handles.D.inv{handles.D.val}.mesh.tess_mni.vert;
    set(handles.Template,  'Value',1);
    set(handles.Individual,'Value',0);
end
handles.grayc = sqrt(sum((handles.vert.^2)')); handles.grayc = handles.grayc'/max(handles.grayc);
set(handles.fig1,'vertices',handles.vert,'faces',handles.face);
UpDate_Display_SRCS(hObject,handles);
axes(handles.sources_axes);
axis image;
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLD SLIDERS - SOURCE LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% upper threshold
% --- Executes on slider movement.
function slider_srcs_up_Callback(hObject, eventdata, handles)
Set_colormap(hObject, eventdata, handles);

%%% lower threshold
% --- Executes on slider movement.
function slider_srcs_down_Callback(hObject, eventdata, handles)
Set_colormap(hObject, eventdata, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSPARENCY SLIDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider_transparency_Callback(hObject, eventdata, handles)
Transparency = get(hObject,'Value');
set(handles.fig1,'facealpha',Transparency);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALISE VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in checkbox_norm.
function checkbox_norm_Callback(hObject, eventdata, handles)
UpDate_Display_SRCS(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE ABSOLUTE VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in checkbox_absv.
function checkbox_absv_Callback(hObject, eventdata, handles)
UpDate_Display_SRCS(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SENSOR LOCATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in checkbox_sensloc.
function checkbox_sensloc_Callback(hObject, eventdata, handles)
try
    if get(hObject,'Value')
        set(handles.sensor_loc,'Visible','on');
    else
        set(handles.sensor_loc,'Visible','off');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME SLIDER - SOURCE & SENSOR LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
ST  = fix(handles.pst(fix(get(hObject,'Value'))));
set(handles.time_bin,'String',num2str(ST));

% Source and sensor space update
%--------------------------------------------------------------------------
UpDate_Display_SRCS(hObject,handles);
UpDate_Display_SENS(hObject,handles);

% --- Callback function
function time_bin_Callback(hObject, eventdata, handles)
[i ST] = min(abs(handles.pst - str2double(get(hObject,'String'))));
set(handles.slider_time,'Value',fix(ST));

% Source and sensor space update
%--------------------------------------------------------------------------
UpDate_Display_SRCS(hObject,handles);
UpDate_Display_SENS(hObject,handles);

% --- Executes on button press in movie.
%--------------------------------------------------------------------------
function movie_Callback(hObject, eventdata, handles)
for t = 1:length(handles.pst)
    set(handles.slider_time,'Value',t);
    ST  = fix(handles.pst(t));
    set(handles.time_bin,'String',num2str(ST));
    UpDate_Display_SRCS(hObject,handles);
end
UpDate_Display_SENS(hObject,handles);


% --- Executes on button press in movie_sens.
%--------------------------------------------------------------------------
function movie_sens_Callback(hObject, eventdata, handles)
for t = 1:length(handles.pst)
    set(handles.slider_time,'Value',t);
    ST  = fix(handles.pst(t));
    set(handles.time_bin,'String',num2str(ST));
    UpDate_Display_SENS(hObject,handles);
end
UpDate_Display_SRCS(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TYPE OF SENSOR LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in sens_display.
function sens_display_Callback(hObject, eventdata, handles)
TypOfDisp = get(handles.sens_display,'Value');

% if time series
%--------------------------------------------------------------------------
if TypOfDisp == 2
    set(handles.checkbox_sensloc,'Value',0);
    set(handles.checkbox_sensloc,'Enable','off');
else
    set(handles.checkbox_sensloc,'Value',1);
    set(handles.checkbox_sensloc,'Enable','on');
end
UpDate_Display_SENS(hObject,handles);


% --- Executes on button press in Exit.
%--------------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
spm_eeg_inv_visu3D_api_OutputFcn(hObject, eventdata, handles);
close(handles.fig);


% --- Executes on button press in Mip.
%--------------------------------------------------------------------------
function Mip_Callback(hObject, eventdata, handles)
ActToDisp = get(handles.Activity,'Value');
if get(handles.Activity,'Value') == 1
    PST = str2num(get(handles.time_bin,'String'));
    spm_eeg_invert_display(handles.D,PST);
else
    spm_eeg_inV_results_display(handles.D);
end

% --- Outputs from this function are returned to the command line.
%--------------------------------------------------------------------------
function varargout = spm_eeg_inv_visu3D_api_OutputFcn(hObject, eventdata, handles)
D = handles.D;
if nargout == 1
    varargout{1} = D;
end


% --- rest threshold
%--------------------------------------------------------------------------
function Set_colormap(hObject, eventdata, handles)
NewMap  = jet(128);

% unsigned values
%--------------------------------------------------------------------------
if get(handles.checkbox_absv,'Value') || get(handles.Activity,'Value') == 3
    
    UpTh    = get(handles.slider_srcs_up,  'Value');
    N       = length(NewMap);
    Low     = fix(N*UpTh);
    Hig     = fix(N - N*UpTh);
    i       = [ones(Low,1); [1:Hig]'*N/Hig];
    NewMap  = NewMap(fix(i),:);
    
% signed values
%--------------------------------------------------------------------------
else

    UpTh    =     get(handles.slider_srcs_up,  'Value');
    DoTh    = 1 - get(handles.slider_srcs_down,'Value');
    N       = length(NewMap)/2;
    Low     = fix(N - N*DoTh);
    Hig     = fix(N - N*UpTh);
    i       = [[1:Low]'*N/Low; ones(N + N - Hig - Low,1)*N; [1:Hig]'*N/Hig + N];
    NewMap  = NewMap(fix(i),:);

end
colormap(NewMap);


% --- Executes on button press in next.
%--------------------------------------------------------------------------
function next_Callback(hObject, eventdata, handles)
if handles.D.val < length(handles.D.inv)
    handles.D.val = handles.D.val + 1;
end
spm_eeg_inv_visu3D_api_OpeningFcn(hObject, eventdata, handles)

% --- Executes on button press in previous.
%--------------------------------------------------------------------------
function previous_Callback(hObject, eventdata, handles)
if handles.D.val > 1
    handles.D.val = handles.D.val - 1;
end
spm_eeg_inv_visu3D_api_OpeningFcn(hObject, eventdata, handles)

        
% --- Executes on button press in VDE.
%--------------------------------------------------------------------------
function Velec_Callback(hObject, eventdata, handles)
axes(handles.sources_axes);
try
    vde = getCursorInfo(handles.location);
catch
    vde = [];
end
if ~length(vde)
    handles.location = datacursormode(handles.figure1);
    set(handles.location,'Enable','on','DisplayStyle','datatip','SnapToDataVertex','on');
    waitforbuttonpress
    datacursormode off
end
vde = getCursorInfo(handles.location);
spm_eeg_invert_display(handles.D,vde.Position)
guidata(hObject,handles);


% --- Executes on button press in Rot.
function Rot_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
axes(handles.sources_axes);
rotate3d


