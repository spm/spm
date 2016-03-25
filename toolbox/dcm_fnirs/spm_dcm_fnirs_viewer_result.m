function varargout = spm_dcm_fnirs_viewer_result(varargin)
% GUI for displaying DCM-fNIRS results
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_dcm_fnirs_viewer_result.m 6754 2016-03-25 06:44:58Z will $

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @spm_dcm_fnirs_viewer_result_OpeningFcn, ...
    'gui_OutputFcn',  @spm_dcm_fnirs_viewer_result_OutputFcn, ...
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


% --- Executes just before spm_dcm_fnirs_viewer_result is made visible.
function spm_dcm_fnirs_viewer_result_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for spm_dcm_fnirs_viewer_result
handles.output = hObject;
cla
%--------------------------------------------------------------------------
% Load DCM.mat 
if isempty(varargin)
    [P, sts] = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; end
    load(P);
else
    DCM = varargin{1,1};
end

%--------------------------------------------------------------------------
% calculate DCM parameters using estimated latent variables 
[A, B, C] = spm_dcm_fnirs_params(DCM);

%--------------------------------------------------------------------------
% obtain predicted data (DCM fit y)
M = DCM.M; 
M.l = M.l + DCM.n; 

y = feval(M.IS, DCM.Ep, M, DCM.U); 
time = linspace(0, DCM.v./DCM.Y.P.fs, DCM.v); 

% information about time window 
time_str = {}; 
if strcmpi(DCM.Y.W.type, 'y')
    time_str{1,1} = 'time window [begin end]';
    nw = size(DCM.Y.W.names, 2); 
    t0 = 0; 
    for i = 1:nw 
        time_str = [time_str; sprintf('[%4.2f %4.2f] [sec]: %s', t0, t0+DCM.Y.W.durations(i), DCM.Y.W.names{i})];
        t0 = t0 + DCM.Y.W.durations(i); 
    end
end

%--------------------------------------------------------------------------
% obtain fNIRS measurements 
Y = DCM.Y.y; % optical density changes 
X0 = DCM.Y.X0; 
Y = Y - X0 * (pinv(DCM.Y.X0) * (Y-y(:, 1:DCM.M.l))); 

SSt = sum((Y - ones(DCM.v, 1) * mean(Y)).^2);
SSr = sum((Y - y(:,1:DCM.M.l)).^2);
R2 = 1 - SSr./SSt; % coefficient of determination

nwav = size(DCM.Y.P.wav, 2); % number of wavelengths 
R2 = mean(reshape(R2, DCM.Y.P.nch, nwav), 2); 
indx_ch = find(R2 > 0.3); 

%--------------------------------------------------------------------------
% for gui
dcm_str{1,1} = 'A - Effective Connectivity in the Absence of Input';
dcm_str{2,1} = 'B - Effective Connectivity Modulated by Input'; 
dcm_str{3,1} = 'C - Influence of Input on Regional Activity'; 
dcm_str{4,1} = 'Estimated Neural Response'; 
dcm_str{5,1} = 'DCM Fit to fNIRS Measurements'; 
radios = {'radio_A', 'radio_B', 'radio_C', 'radio_neural', 'radio_fit'}; 

%--------------------------------------------------------------------------
% save variables in handle structure 
handles.DCM = DCM;
handles.A = A; 
handles.B = B; 
handles.C = C; 

handles.y = y; 
handles.Y = Y; 
handles.nwav = nwav; 
handles.indx_ch = indx_ch; 
handles.time = time; 
handles.time_str = time_str; 

handles.dcm_str = dcm_str; 
handles.radios = radios; 

%--------------------------------------------------------------------------
% display connectivity (B) modulated by input 
dim = size(DCM.b); 
indx = find(DCM.b == 1); 
n = size(indx, 1); 
[I, J, K] = ind2sub(dim, indx); 

iname = {}; 
conn_str = {};
for i = 1:n
    iname = [iname; sprintf('%d,%d,%d', I(i), J(i), K(i))]; 
    conn_str = [conn_str; sprintf('%s : connectivity from %s to %s (modulated by %s)', iname{i,1}, DCM.xY(J(i)).name, DCM.xY(I(i)).name, DCM.U.name{K(i)});]; 
end

axes(handles.axes_signal);
bar(B.mean(indx), 'FaceColor', [0.25 0.25 0.25]);
hold on
errorbar(B.mean(indx), B.std(indx), '.', 'color', [1 0 0]);
title(dcm_str{2,1}, 'FontWeight', 'bold'); 
set(handles.axes_signal, 'XTick',1:n, 'XTickLabel',iname)
hold off 

%--------------------------------------------------------------------------
% for gui
dcm_p = 2; 
nr = size(radios, 2); 
h = zeros(nr, 1); h(dcm_p) = 1; 
for i = 1:nr, set(handles.(handles.radios{i}), 'value', h(i)); end 

set(handles.listbox_dcm, 'string', conn_str, 'value', 1);
set(handles.slider_ch, 'enable', 'off');
if isempty(indx_ch), set(handles.radio_fit, 'enable', 'off'); end 

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = spm_dcm_fnirs_viewer_result_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_ch_Callback(hObject, eventdata, handles)
num_s = round(get(handles.slider_ch, 'value')); set(handles.slider_ch, 'value', num_s); 

h = get(handles.radio_neural, 'value'); 
DCM = handles.DCM;
time = handles.time;
y = handles.y; 

axes(handles.axes_signal);
if h % plot neural responses 
    dcm_p = 4; 
    y(:, 1:DCM.M.l) = [];
    plot(time, y(:, num_s), 'LineWidth', 4);
    title(sprintf('%s in %s', handles.dcm_str{dcm_p,1}, DCM.xY(num_s).name), 'FontWeight', 'bold');
else % plot DCM fit
    dcm_p = 5; 
    nwav = handles.nwav;
    Y = reshape(handles.Y, DCM.v, DCM.Y.P.nch, nwav);
    y = reshape(y(:, 1:DCM.M.l), DCM.v, DCM.Y.P.nch, nwav);
    indx_ch = handles.indx_ch;
    
    mcolor{1} = [0 0 0]; mcolor{2} = [0 0 1]; mcolor{3} = [0 1 0];
    for i = 1:nwav
        plot(time, Y(:, indx_ch(num_s), i),'Color', mcolor{i}, 'LineWidth', 1)
        hold on;
        plot(time, y(:, indx_ch(num_s), i), 'Color', [1 0 0], 'LineWidth', 3);
    end
    title(sprintf('%s (Ch %d)', handles.dcm_str{dcm_p,1}, DCM.Y.P.rois(indx_ch(num_s))), 'FontWeight', 'bold');
end
xlabel('Time (s)');
axis tight;
hold off;

% --- Executes during object creation, after setting all properties.
function slider_ch_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radio_neural.
function radio_neural_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
%-Set radio button properties 
dcm_p = 4; 
nr = size(handles.radios, 2); 
h = zeros(nr, 1); h(dcm_p) = 1; 
for i = 1:nr, set(handles.(handles.radios{i}), 'value', h(i)); end 

%--------------------------------------------------------------------------
%- display estimated neural response 
DCM = handles.DCM; 
y = handles.y; y(:, 1:DCM.M.l) = []; 
time = handles.time; 

axes(handles.axes_signal);
plot(time, y(:, 1), 'LineWidth', 4); 
title(sprintf('%s in %s', handles.dcm_str{dcm_p,1}, DCM.xY(1).name), 'FontWeight', 'bold');
xlabel('Time (s)');
axis tight;

% update strings for listbox 
set(handles.listbox_dcm, 'string', handles.time_str, 'value', 1);
set(handles.slider_ch, 'sliderstep', [1/(DCM.n-1), 1/(DCM.n-1)], 'max', DCM.n, 'min', 1, 'value', 1, 'enable', 'on');

% --- Executes on button press in radio_fit.
function radio_fit_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
%-Set radio button properties 
dcm_p = 5; 
nr = size(handles.radios, 2); 
h = zeros(nr, 1); h(dcm_p) = 1; 
for i = 1:nr, set(handles.(handles.radios{i}), 'value', h(i)); end 

%--------------------------------------------------------------------------
%- display DCM fit to fNIRS measurements 
DCM = handles.DCM; 
nwav = handles.nwav; 
Y = reshape(handles.Y, DCM.v, DCM.Y.P.nch, nwav); 
y = reshape(handles.y(:, 1:DCM.M.l), DCM.v, DCM.Y.P.nch, nwav); 
time = handles.time; 
indx_ch = handles.indx_ch; 

mcolor{1} = [0 0 0]; mcolor{2} = [0 0 1]; mcolor{3} = [0 1 0]; 
axes(handles.axes_signal);
for i = 1:nwav
    plot(time, Y(:, indx_ch(1), i),'Color', mcolor{i}, 'LineWidth', 1)
    hold on
    plot(time, y(:, indx_ch(1), i), 'Color', [1 0 0], 'LineWidth', 3);
end
xlabel('Time (s)');
axis tight;
title(sprintf('%s (Ch %d)', handles.dcm_str{dcm_p,1}, DCM.Y.P.rois(indx_ch(1))), 'FontWeight', 'bold');
hold off;

% update strings for listbox 
set(handles.listbox_dcm, 'string', handles.time_str, 'value', 1);
set(handles.slider_ch, 'sliderstep', [1/(length(indx_ch)-1), 1/(length(indx_ch)-1)], 'max', length(indx_ch), 'min', 1, 'value', 1, 'enable', 'on');

% --- Executes on selection change in listbox_dcm.
function listbox_dcm_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function listbox_dcm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_A.
function radio_A_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
%-Set radio button properties 
dcm_p = 1; 
nr = size(handles.radios, 2); 
h = zeros(nr, 1); h(dcm_p) = 1; 
for i = 1:nr, set(handles.(handles.radios{i}), 'value', h(i)); end 

%--------------------------------------------------------------------------
%- display estimated DCM parameter A 
DCM = handles.DCM; 

dim = size(DCM.a); 
indx = find(DCM.a == 1); 
n = size(indx, 1); 
[I, J] = ind2sub(dim, indx); 

axes(handles.axes_signal);
bar(handles.A.mean(indx), 'FaceColor', [0.25 0.25 0.25]);
hold on;
errorbar(handles.A.mean(indx), handles.A.std(indx), '.', 'color', [1 0 0]);
title(handles.dcm_str{dcm_p,1}, 'FontWeight', 'bold'); 

% update strings for listbox 
iname = {}; 
conn_str = {};

for i = 1:n
    iname = [iname; sprintf('%d,%d', I(i), J(i))]; 
    conn_str = [conn_str; sprintf('%s : connectivity from %s to %s', iname{i,1}, DCM.xY(J(i)).name, DCM.xY(I(i)).name)]; 
end

set(handles.axes_signal, 'XTick',1:n, 'XTickLabel',iname)
set(handles.listbox_dcm, 'string', conn_str, 'value', 1);
set(handles.slider_ch, 'enable', 'off');
hold off; 

% --- Executes on button press in radio_B.
function radio_B_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
%-Set radio button properties 
dcm_p = 2; 
nr = size(handles.radios, 2); 
h = zeros(nr, 1); h(dcm_p) = 1; 
for i = 1:nr, set(handles.(handles.radios{i}), 'value', h(i)); end 

%--------------------------------------------------------------------------
%- display estimated DCM parameter B 
DCM = handles.DCM; 

dim = size(DCM.b); 
indx = find(DCM.b == 1); 
n = size(indx, 1); 
[I, J, K] = ind2sub(dim, indx); 

axes(handles.axes_signal);
bar(handles.B.mean(indx), 'FaceColor', [0.25 0.25 0.25]); 
hold on;
errorbar(handles.B.mean(indx), handles.B.std(indx), '.', 'color', [1 0 0]);
title(handles.dcm_str{dcm_p,1}, 'FontWeight', 'bold'); 

% update strings for listbox 
iname = {}; 
conn_str = {};

for i = 1:n
    iname = [iname; sprintf('%d,%d,%d', I(i), J(i), K(i))];
    conn_str = [conn_str; sprintf('%s : connectivity from %s to %s (modulated by %s)', iname{i,1}, DCM.xY(J(i)).name, DCM.xY(I(i)).name, DCM.U.name{K(i)});];
end

set(handles.axes_signal, 'XTick',1:n, 'XTickLabel',iname)
set(handles.listbox_dcm, 'string', conn_str, 'value', 1);
set(handles.slider_ch, 'enable', 'off');
hold off; 


% --- Executes on button press in radio_C.
function radio_C_Callback(hObject, eventdata, handles)
%--------------------------------------------------------------------------
%-Set radio button properties 
dcm_p = 3; 
nr = size(handles.radios, 2); 
h = zeros(nr, 1); h(dcm_p) = 1; 
for i = 1:nr, set(handles.(handles.radios{i}), 'value', h(i)); end 

%--------------------------------------------------------------------------
%- display estimated DCM parameter C 
DCM = handles.DCM; 

dim = size(DCM.c); 
indx = find(DCM.c == 1); 
n = size(indx, 1); 
[I, J] = ind2sub(dim, indx); 

axes(handles.axes_signal);
bar(handles.C.mean(indx), 'FaceColor', [0.25 0.25 0.25]);
hold on;
errorbar(handles.C.mean(indx), handles.C.std(indx), '.', 'color', [1 0 0]);
title(handles.dcm_str{dcm_p,1}, 'FontWeight', 'bold'); 

% update strings for listbox 
iname = {}; 
conn_str = {};

for i = 1:n
    iname = [iname; sprintf('%d,%d', I(i), J(i))];
    conn_str = [conn_str; sprintf('%s : influence of input (%s) on %s', iname{i,1}, DCM.U.name{J(i)}, DCM.xY(I(i)).name)]; 
end

set(handles.axes_signal, 'XTick',1:n, 'XTickLabel',iname)
set(handles.listbox_dcm, 'string', conn_str, 'value', 1);
set(handles.slider_ch, 'enable', 'off');
hold off; 
