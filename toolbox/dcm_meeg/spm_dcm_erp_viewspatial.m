function varargout = spm_dcm_erp_viewspatial(varargin)
% SPM_DCM_ERP_VIEWSPATIAL M-file for spm_dcm_erp_viewspatial.fig
%      SPM_DCM_ERP_VIEWSPATIAL, by itself, creates a new SPM_DCM_ERP_VIEWSPATIAL or raises the existing
%      singleton*.
%
%      H = SPM_DCM_ERP_VIEWSPATIAL returns the handle to a new SPM_DCM_ERP_VIEWSPATIAL or the handle to
%      the existing singleton*.
%
%      SPM_DCM_ERP_VIEWSPATIAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPM_DCM_ERP_VIEWSPATIAL.M with the given input arguments.
%
%      SPM_DCM_ERP_VIEWSPATIAL('Property','Value',...) creates a new SPM_DCM_ERP_VIEWSPATIAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_dcm_erp_viewspatial_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spm_dcm_erp_viewspatial_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spm_dcm_erp_viewspatial

% Last Modified by GUIDE v2.5 12-Jul-2006 11:06:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_dcm_erp_viewspatial_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_dcm_erp_viewspatial_OutputFcn, ...
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


% --- Executes just before spm_dcm_erp_viewspatial is made visible.
function spm_dcm_erp_viewspatial_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_dcm_erp_viewspatial (see VARARGIN)

% Choose default command line output for spm_dcm_erp_viewspatial
handles.output = hObject;

DCM = varargin{1};
handles.DCM = DCM;
handles.ms  = DCM.xY.pst;
handles.T   = 1;

M = DCM.M;

load(DCM.xY.Dfile); % ----> returns SPM/EEG struct D
handles.D = D;

% locations for plotting
CTF = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));
CTF.Cpos = CTF.Cpos(:, D.channels.order(D.channels.eeg));

handles.x = min(CTF.Cpos(1,:)):0.01:max(CTF.Cpos(1,:));
handles.y = min(CTF.Cpos(2,:)):0.01:max(CTF.Cpos(2,:));

[handles.x1, handles.y1] = meshgrid(handles.x, handles.y);
handles.xp = CTF.Cpos(1,:)';
handles.yp = CTF.Cpos(2,:)';

Nchannels = size(CTF.Cpos, 2);
handles.Nchannels = Nchannels;

handles.y_proj   = spm_cat(DCM.xY.y)*DCM.M.E;
handles.CLim1_yp = min(min(handles.y_proj));
handles.CLim2_yp = max(max(handles.y_proj));

handles.Nt = size(handles.y_proj, 1);

% data and model fit
handles.yd = NaN*ones(handles.Nt, Nchannels);
handles.ym = NaN*ones(handles.Nt, Nchannels);
handles.yd(:, DCM.M.dipfit.Ic) = handles.y_proj*DCM.M.E'; % data (back-projected to channel space)
handles.ym(:, DCM.M.dipfit.Ic) = cat(1,DCM.H{:})*DCM.M.E'; % model fit

handles.CLim1 = min(min([handles.yd handles.ym]));
handles.CLim2 = max(max([handles.yd handles.ym]));

% set slider's range and initial value
set(handles.slider1, 'min', 1);
set(handles.slider1, 'max', handles.Nt);
set(handles.slider1, 'Value', 1);
set(handles.slider1, 'Sliderstep', [1/(handles.Nt-1) 10/(handles.Nt-1)]); % moves slider in dt and 10*dt steps

plot_images(hObject, handles);
plot_modes(hObject, handles);
plot_dipoles(hObject, handles);
plot_components_space(hObject, handles);
try
    plot_components_time(hObject, handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spm_dcm_erp_viewspatial wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spm_dcm_erp_viewspatial_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.T = round(get(handles.slider1, 'Value'));

plot_images(hObject, handles);
plot_modes(hObject, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%--------------------------------------------------------------------
function plot_images(hObject, handles)

T = handles.T;

axes(handles.axes1); cla
z = griddata(handles.xp, handles.yp, handles.yd(T,:), handles.x1, handles.y1);
surface(handles.x, handles.y, z);
axis off
axis square
shading('interp')
hold on
plot3(handles.xp, handles.yp, handles.yd, 'k.');
set(handles.axes1, 'CLim', [handles.CLim1 handles.CLim2])
title('data', 'FontSize', 16);

axes(handles.axes2); cla
z = griddata(handles.xp, handles.yp, handles.ym(T,:), handles.x1, handles.y1);
surface(handles.x, handles.y, z);
axis off
axis square
shading('interp')
hold on
plot3(handles.xp, handles.yp, handles.ym, 'k.');
set(handles.axes2, 'CLim', [handles.CLim1 handles.CLim2])
title('model', 'FontSize', 16);

guidata(hObject, handles);

drawnow
%--------------------------------------------------------------------
function plot_modes(hObject, handles)
% plot temporal expression of modes

DCM     = handles.DCM;
Nt      = size(DCM.xY.xy{1}, 1);
Ntrials = length(DCM.H);


ms_all  = kron(ones(1, Ntrials), handles.ms);

% data and model prediction, cond 1
axes(handles.axes3); cla
plot(ms_all, handles.y_proj);
hold on
plot(ms_all, cat(1, DCM.H{:}), '--');
plot([ms_all(handles.T) ms_all(handles.T)], [handles.CLim1_yp handles.CLim2_yp], 'k', 'LineWidth', 3);
% set(handles.axes3, 'XLim', [0 handles.ms(end)]);
set(handles.axes3, 'YLim', [handles.CLim1_yp handles.CLim2_yp]);
xlabel('ms');
title('Temporal expressions of modes', 'FontSize', 16);
grid on

%--------------------------------------------------------------------
function plot_dipoles(hObject, handles)

DCM = handles.DCM;
Nsources = length(DCM.M.pE.A{1});
% ECD
%--------------------------------------------------------------------
try
    Lpos = handles.DCM.Eg.Lpos;
    Lmom = handles.DCM.Eg.L;
    Lmom = spm_cat(Lmom);
    
% Imaging
%--------------------------------------------------------------------
catch
    Lpos = handles.DCM.M.dipfit.Lpos;
    Lmom = Lpos*0;
end
try
    elc = handles.DCM.M.dipfit.elc;
    % transform sensor locations to MNI-space
    iMt = inv(DCM.M.dipfit.Mmni2polsphere);
    elc = iMt*[elc'; ones(1, size(elc, 1))];
    elc = elc(1:3, :)';
catch
    elc = sparse(0,3);
end



axes(handles.axes5);
for j = 1:Nsources
    % plot dipoles using small ellipsoids
    [x, y, z] = ellipsoid(handles.axes5,Lpos(1,j), Lpos(2,j), Lpos(3,j), 4, 4,4);
    surf(x, y, z, 'EdgeColor', 'none', 'FaceColor', 'g');
    hold on
    
    % plot dipole moments
    plot3([Lpos(1, j) Lpos(1, j) + 5*Lmom(1, j)],...
          [Lpos(2, j) Lpos(2, j) + 5*Lmom(2, j)],...
          [Lpos(3, j) Lpos(3, j) + 5*Lmom(3, j)], 'b', 'LineWidth', 4);
    
    plot3(elc(:,1), elc(:,2), elc(:,3), 'r.','MarkerSize',18);

end
xlabel('x'); ylabel('y');
rotate3d(handles.axes5);

axis equal

%--------------------------------------------------------------------
function plot_components_space(hObject, handles)
% plots spatial expression of each dipole
DCM = handles.DCM;
Nsources = length(DCM.M.pE.A{1});

lf  = NaN*ones(handles.Nchannels, Nsources);
lfo = NaN*ones(handles.Nchannels, Nsources);
x   = [0 kron([zeros(1, 8) 1], ones(1, Nsources))];

% unprojected and projected leadfield
%--------------------------------------------------------------------
lfo(DCM.M.dipfit.Ic,:) = spm_erp_L(DCM.Eg,DCM.M);
lf                     = DCM.M.E*DCM.M.E'*lfo;


% use subplots in new figure
h = figure;
set(h, 'Name', 'Spatial expression of sources');
for j = 1:2
    for i = 1:Nsources
        if j == 1
            % projected
            handles.hcomponents{i} = subplot(2,Nsources, i);
            z = griddata(handles.xp, handles.yp, lf(:, i), handles.x1, handles.y1);
            title(DCM.Sname{i}, 'FontSize', 16);
        
        else
            % not projected
            handles.hcomponents{i} = subplot(2,Nsources, i+Nsources);
            z = griddata(handles.xp, handles.yp, lfo(:, i), handles.x1, handles.y1);
        end            
        surface(handles.x, handles.y, z);
        axis off
        axis square
        shading('interp')
        %    set(handles.axes1, 'CLim', [handles.CLim1 handles.CLim2])
    end
end

drawnow

function plot_components_time(hObject, handles)

DCM = handles.DCM;
Lmom = sqrt(sum(spm_cat(DCM.Eg.L)).^2);
Nsources = length(DCM.M.pE.A{1});

h = figure;
set(h, 'Name', 'Effective source amplitudes');
% multiply with dipole amplitudes
for i = 1:Nsources
    handles.hcomponents_time{2*(i-1)+1} = subplot(Nsources, 2, 2*(i-1)+1);
    plot(handles.ms, Lmom(i)*DCM.K{1}(:, i));
    title([handles.DCM.Sname{i} ', ERP 1'], 'FontSize', 16);
    
    try
    handles.hcomponents_time{2*(i-1)+2} = subplot(Nsources, 2, 2*(i-1)+2);
    plot(handles.ms, Lmom(i)*DCM.K{2}(:, i));
    title([handles.DCM.Sname{i} ', ERP 2'], 'FontSize', 16);
    end
end


