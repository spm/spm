function varargout = spm_eeg_inv_visu3D_api(varargin)
% SPM_EEG_INV_VISU3D_API M-file for spm_eeg_inv_visu3D_api.fig
% - FIG = SPM_EEG_INV_VISU3D_API launch spm_eeg_inv_visu3D_api GUI.
% - D   = SPM_EEG_INV_VISU3D_API(D) open with D
% - SPM_EEG_INV_VISU3D_API(filename) where filename is the eeg/meg .mat file
% - SPM_EEG_INV_VISU3D_API('callback_name', ...) invoke the named callback.
%
% Last Modified by GUIDE v2.5 18-Feb-2011 14:23:27
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_visu3D_api.m 4818 2012-07-31 14:53:10Z guillaume $

% INITIALISATION CODE
%--------------------------------------------------------------------------
if nargin == 0 || nargin == 1  % LAUNCH GUI
    
    
    % open new api
    %----------------------------------------------------------------------
    fig         = openfig(mfilename,'new');
    handles     = guihandles(fig);
    handles.fig = fig;
    guidata(fig,handles);
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    
    % load D if possible and try to open
    %----------------------------------------------------------------------
    try
        handles.D = spm_eeg_inv_check(varargin{1});
        set(handles.DataFile,'String',handles.D.fname)
        spm_eeg_inv_visu3D_api_OpeningFcn(fig, [], handles)
    end
    
    % return figure handle if necessary
    %----------------------------------------------------------------------
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
        disp(lasterror);
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
    D = spm_eeg_load(spm_select(1, '.mat', 'Select EEG/MEG mat file'));
end

if ~isfield(D,'inv')
    error('Please specify and invert a forward model\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET RESULTS (default: current or last analysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(handles.fig);
axes(handles.sensors_axes);

try, val = D.val; catch,  val = 1; D.val = 1; end
try, con = D.con; catch,  con = 1; D.con = 1; end
if (D.con == 0) || (D.con > length(D.inv{D.val}.inverse.J))
    con = 1; D.con = 1;
end
handles.D = D;

set(handles.DataFile,'String',D.fname);
set(handles.next,'String',sprintf('model %i',val));
set(handles.con, 'String',sprintf('condition %i',con));
set(handles.fig,'name',['Source visualisation -' D.fname])

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
    Is                = D.inv{val}.inverse.Is;
    L                 = D.inv{val}.inverse.L;
    U                 = D.inv{val}.inverse.U;
    T                 = D.inv{val}.inverse.T;
    Y                 = D.inv{val}.inverse.Y{con};
    Ts                = ceil(linspace(1,size(T,1),dimT));
    
    % source data
    %----------------------------------------------------------------------
    set(handles.Activity,'Value',1);
    J                 = sparse(dimS,dimT);
    J(Is,:)           = D.inv{val}.inverse.J{con}*T(Ts,:)';
    handles.dimT      = dimT;
    handles.dimS      = dimS;
    handles.pst       = D.inv{val}.inverse.pst(Ts);
    handles.srcs_data = J;
    handles.Nmax      = max(abs(J(:)));
    handles.Is        = Is;
    
    % sensor data
    %----------------------------------------------------------------------        
    if ~iscell(U)
        U = {U'};
    end
    
    A             = spm_pinv(spm_cat(spm_diag(U))')';    
    
    handles.sens_data = A*Y*T(Ts,:)';
    handles.pred_data = A*L*J(Is,:);
catch
    warndlg({'Please invert your model';'inverse solution not valid'});
    return
end

% case 'windowed response' or contrast'
%--------------------------------------------------------------------------
try
    JW                    = sparse(dimS,1);
    GW                    = sparse(dimS,1);
    JW(Is,:)              = D.inv{val}.contrast.JW{con};
    GW(Is,:)              = D.inv{val}.contrast.GW{con};
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
handles.grayc = sqrt(sum((vert.^2),2)); handles.grayc = handles.grayc'/max(handles.grayc);
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

srcs_disp = full(abs(handles.srcs_data(:,1)));
handles.fig1 = patch('vertices',handles.vert,'faces',handles.face,'FaceVertexCData',srcs_disp);

% display
%--------------------------------------------------------------------------
set(handles.fig1,'FaceColor',[.5 .5 .5],'EdgeColor','none');
shading interp
lighting gouraud
camlight
zoom off
lightangle(0,270);lightangle(270,0),lightangle(0,0),lightangle(90,0);
material([.1 .1 .4 .5 .4]);
view(140,15);
axis image
handles.colorbar = colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD SENSOR FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ic = {};
if iscell(D.inv{val}.inverse.Ic)
    for i = 1:numel(D.inv{val}.inverse.Ic)
        if i == 1
            Ic{i} = 1:length(D.inv{val}.inverse.Ic{i});
        else
            Ic{i} = Ic{i-1}(end)+(1:length(D.inv{val}.inverse.Ic{i}));
        end
    end
else
   Ic{1} = 1:length(D.inv{val}.inverse.Ic); 
end

handles.Ic = Ic;

coor = D.coor2D(full(spm_cat(D.inv{val}.inverse.Ic)));
xp   = coor(1,:)';
yp   = coor(2,:)';

x        = linspace(min(xp),max(xp),64);
y        = linspace(min(yp),max(yp),64);
[xm,ym]  = meshgrid(x,y);
handles.sens_coord = [xp yp];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL SENSOR LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(D.inv{val}.inverse, 'modality')
  set(handles.modality, 'String', D.inv{val}.inverse.modality);
else
   set(handles.modality, 'String', 'MEEG');  % This is for backward compatibility with old DCM-IMG
end

figure(handles.fig)
axes(handles.sensors_axes);
cla; axis off

im       = get(handles.modality, 'Value');

ic         = handles.Ic{im};
disp       = full(handles.sens_data(ic,1));
imagesc(x,y,griddata(xp(ic),yp(ic),disp,xm,ym));
axis image xy off
handles.sens_coord_x   = x;
handles.sens_coord_y   = y;
handles.sens_coord2D_X = xm;
handles.sens_coord2D_Y = ym;
hold on
handles.sensor_loc = plot(handles.sens_coord(ic,1),handles.sens_coord(ic,2),'o','MarkerFaceColor',[1 1 1]/2,'MarkerSize',6);
set(handles.checkbox_sensloc,'Value',1);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL SENSOR LEVEL DISPLAY - PREDICTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.pred_axes); cla;
disp      = full(handles.pred_data(ic,1));
imagesc(x,y,griddata(xp(ic),yp(ic),disp,xm,ym));
axis image xy off
drawnow

guidata(hObject,handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE SOURCE LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpDate_Display_SRCS(hObject,handles)
axes(handles.sources_axes);

if isfield(handles,'fig1')
    ActToDisp = get(handles.Activity,'Value');
    A         = get(handles.checkbox_absv,'Value');
    N         = get(handles.checkbox_norm,'Value');
    switch ActToDisp
        
        % case 1: response (J)
        %------------------------------------------------------------------
        case 1
            TS = fix(get(handles.slider_time,'Value'));
            if A
                srcs_disp = abs(handles.srcs_data(:,TS));
            else
                srcs_disp = handles.srcs_data(:,TS);
            end
            if N
                if A
                    handles.Vmin =  0;
                    handles.Vmax =  handles.Nmax;
                else
                    handles.Vmin = -handles.Nmax;
                    handles.Vmax =  handles.Nmax;
                end
            else
                handles.Vmin = min(srcs_disp);
                handles.Vmax = max(srcs_disp);
            end
            
            
            % case 2: Windowed response (JW)
            %------------------------------------------------------------------
        case 2
            handles.Vmin = min(handles.srcs_data_w);
            handles.Vmax = max(handles.srcs_data_w);
            srcs_disp    = handles.srcs_data_w;
            
            % case 3: Evoked power  (JWWJ)
            %------------------------------------------------------------------
        case 3
            handles.Vmin = min(handles.srcs_data_ev);
            handles.Vmax = max(handles.srcs_data_ev);
            srcs_disp    = handles.srcs_data_ev;
            
            % case 4: Induced power  (JWWJ)
            %------------------------------------------------------------------
        case 4
            handles.Vmin = min(handles.srcs_data_ind);
            handles.Vmax = max(handles.srcs_data_ind);
            srcs_disp    = handles.srcs_data_ind;
    end
    set(handles.fig1,'FaceVertexCData',full(srcs_disp));
    set(handles.sources_axes,'CLim',[handles.Vmin handles.Vmax]);
    set(handles.sources_axes,'CLimMode','manual');
    
end

% Adjust the threshold
%--------------------------------------------------------------------------
Set_colormap(hObject, [], handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE SENSOR LEVEL DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpDate_Display_SENS(hObject,handles)
TypOfDisp = get(handles.sens_display,'Value');
ActToDisp = get(handles.Activity,'Value');

im       = get(handles.modality, 'Value');
ic       = handles.Ic{im};


% topography
%--------------------------------------------------------------------------
if TypOfDisp == 1
    
    % responses at one pst
    %----------------------------------------------------------------------
    if ActToDisp == 1
        
        TS = fix(get(handles.slider_time,'Value'));
        sens_disp = handles.sens_data(ic,TS);
        pred_disp = handles.pred_data(ic,TS);
        
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
    
    axes(handles.sensors_axes);
    disp = griddata(handles.sens_coord(ic,1),handles.sens_coord(ic,2),full(sens_disp),handles.sens_coord2D_X,handles.sens_coord2D_Y);
    imagesc(handles.sens_coord_x,handles.sens_coord_y,disp);
    axis image xy off
    
    % add sensor locations
    %----------------------------------------------------------------------
    try, delete(handles.sensor_loc); end
    hold(handles.sensors_axes, 'on');
    handles.sensor_loc = plot(handles.sensors_axes,...
        handles.sens_coord(ic,1),handles.sens_coord(ic,2),'o','MarkerFaceColor',[1 1 1]/2,'MarkerSize',6);
    hold(handles.sensors_axes, 'off');
    
    axes(handles.pred_axes);
    disp = griddata(handles.sens_coord(ic,1),handles.sens_coord(ic,2),full(pred_disp),handles.sens_coord2D_X,handles.sens_coord2D_Y);
    imagesc(handles.sens_coord_x,handles.sens_coord_y,disp);
    axis image xy off;
    
    checkbox_sensloc_Callback(hObject, [], handles);
    
    % time series
    %--------------------------------------------------------------------------
elseif TypOfDisp == 2
    axes(handles.sensors_axes)
    daspect('auto')
    handles.fig2 = ...
        plot(handles.pst,handles.sens_data(ic, :),'b-.',handles.pst,handles.pred_data(ic, :),'r:');
    if ActToDisp > 1
        hold on
        Scal = norm(handles.sens_data,1)/norm(handles.W,1);
        plot(handles.pst,handles.W*Scal,'k')
        hold off
    end
    axis on tight;
    axes(handles.pred_axes); cla, axis off
end

% Adjust the threshold
%--------------------------------------------------------------------------
Set_colormap(hObject, [], handles);
guidata(hObject,handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataFile_Callback(hObject, eventdata, handles)
S     = get(handles.DataFile,'String');
try
    D = spm_eeg_load(S);
catch
    LoadData_Callback(hObject, eventdata, handles);
end

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, eventdata, handles)
S         = spm_select(1, '.mat', 'Select EEG/MEG mat file');
handles.D = spm_eeg_load(S);
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
    tess_ctx = gifti(handles.D.inv{handles.D.val}.mesh.tess_ctx);
    handles.vert = tess_ctx.vertices;
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
    if get(handles.checkbox_sensloc,'Value')
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
global MOVIE
for t = 1:length(handles.pst)
    set(handles.slider_time,'Value',t);
    ST  = fix(handles.pst(t));
    set(handles.time_bin,'String',num2str(ST));
    UpDate_Display_SRCS(hObject,handles);
    
    % record movie if requested
    %----------------------------------------------------------------------
    if MOVIE, M(t) = getframe(handles.sources_axes); end;
end
UpDate_Display_SENS(hObject,handles);
try
    filename = fullfile(handles.D.path,'SourceMovie');
    movie2avi(M,filename,'compression','Indeo3','FPS',24)
end



% --- Executes on button press in movie_sens.
%--------------------------------------------------------------------------
function movie_sens_Callback(hObject, eventdata, handles)
global MOVIE
for t = 1:length(handles.pst)
    set(handles.slider_time,'Value',t);
    ST  = fix(handles.pst(t));
    set(handles.time_bin,'String',num2str(ST));
    UpDate_Display_SENS(hObject,handles);
    
    % record movie if requested
    %----------------------------------------------------------------------
    if MOVIE, M(t) = getframe(handles.sensors_axes); end;
    
end
UpDate_Display_SRCS(hObject,handles);
try
    filename = fullfile(handles.D.path,'SensorMovie');
    movie2avi(M,filename,'compression','Indeo3','FPS',24)
end

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
    spm_eeg_inv_results_display(handles.D);
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
NewMap  = jet;

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
drawnow


% --- Executes on button press in next.
%--------------------------------------------------------------------------
function next_Callback(hObject, eventdata, handles)
if length(handles.D.inv) == 1
    set(handles.next,'Value',0);
    return
end
handles.D.val = handles.D.val + 1;
handles.D.con = 1;
if handles.D.val > length(handles.D.inv)
    handles.D.val = 1;
end
set(handles.next,'String',sprintf('model %d',handles.D.val),'Value',0);
spm_eeg_inv_visu3D_api_OpeningFcn(hObject, eventdata, handles)

% --- Executes on button press in previous.
%--------------------------------------------------------------------------
function con_Callback(hObject, eventdata, handles)
if length(handles.D.inv{handles.D.val}.inverse.J) == 1
    set(handles.con,'Value',0);
    return
end

handles.D.con = handles.D.con + 1;
if handles.D.con > length(handles.D.inv{handles.D.val}.inverse.J)
    handles.D.con = 1;
end
set(handles.con,'String',sprintf('condition %d',handles.D.con),'Value',0);
spm_eeg_inv_visu3D_api_OpeningFcn(hObject, eventdata, handles)

% --- Executes on button press in VDE.
%--------------------------------------------------------------------------
function Velec_Callback(hObject, eventdata, handles)
axes(handles.sources_axes);
rotate3d off;
datacursormode off;

if isfield(handles, 'velec')
    delete(handles.velec);
    handles = rmfield(handles, 'velec');
end

set(handles.fig1, 'ButtonDownFcn', @Velec_ButtonDown)

guidata(hObject,handles);

function Velec_ButtonDown(hObject, eventdata)
handles     = guidata(hObject);
vert        = handles.vert(handles.Is, :);
coord       = get(handles.sources_axes, 'CurrentPoint');
dist        = sum((vert - repmat(coord(1, :), size(vert, 1), 1)).^2, 2);
[junk, ind] = min(dist);
coord       = vert(ind, :);

axes(handles.sources_axes);
hold on
handles.velec = plot3(coord(1), coord(2), coord(3), 'rv', 'MarkerSize', 10);

spm_eeg_invert_display(handles.D, coord);

set(handles.fig1, 'ButtonDownFcn', '');
guidata(hObject,handles);


% --- Executes on button press in Rot.
function Rot_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
rotate3d(handles.sources_axes)
return


% --- Executes on selection change in modality.
function modality_Callback(hObject, eventdata, handles)
% hObject    handle to modality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modality contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modality
UpDate_Display_SENS(hObject,handles)
