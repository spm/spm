function varargout = spm_eeg_firstlevel(varargin)
% SPM_EEG_FIRSTLEVEL M-file for spm_eeg_firstlevel.fig
%      SPM_EEG_FIRSTLEVEL, by itself, creates a new SPM_EEG_FIRSTLEVEL or raises the existing
%      singleton*.
%
%      H = SPM_EEG_FIRSTLEVEL returns the handle to a new SPM_EEG_FIRSTLEVEL or the handle to
%      the existing singleton*.
%
%      SPM_EEG_FIRSTLEVEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPM_EEG_FIRSTLEVEL.M with the given input arguments.
%
%      SPM_EEG_FIRSTLEVEL('Property','Value',...) creates a new SPM_EEG_FIRSTLEVEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_eeg_firstlevel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spm_eeg_firstlevel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help spm_eeg_firstlevel

% Last Modified by GUIDE v2.5 02-Feb-2006 16:29:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_eeg_firstlevel_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_eeg_firstlevel_OutputFcn, ...
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


% --- Executes just before spm_eeg_firstlevel is made visible.
function spm_eeg_firstlevel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_eeg_firstlevel (see VARARGIN)

% Choose default command line output for spm_eeg_firstlevel
handles.output = hObject;

handles.con_types = {};
handles.Nt = 0;
set(handles.listbox1, 'String', {});
set(handles.listbox2, 'String', {});

spm_defaults;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spm_eeg_firstlevel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spm_eeg_firstlevel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addfile.
function addfile_Callback(hObject, eventdata, handles)
% hObject    handle to addfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname, filterindex] = uigetfile('*.img', 'Select an image');

if filename == 0, return, end

[status, Iimg, Limg] = spm_eeg_firstlevel_checkfile(fullfile(pathname, filename), handles.Nt);

if status == 0
    errordlg(sprintf('Can''t add: new file has %d time points, selected have %d', Limg, handles.Nt));
    return
else
    handles.Nt = status;
end

% add image to list
set(handles.listbox1, 'String', [get(handles.listbox1, 'String'); {fullfile(pathname, filename)}]);

set(handles.Dfile, 'Enable', 'on');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in removefile.
function removefile_Callback(hObject, eventdata, handles)
% hObject    handle to removefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(get(handles.listbox1, 'String')) == 0, return, end

ind = get(handles.listbox1, 'Value');
tmp = get(handles.listbox1, 'String');
tmp(ind) = [];
set(handles.listbox1, 'Value', max(1, ind-1));
set(handles.listbox1, 'String', tmp);

% user removed all files?
if length(get(handles.listbox1, 'String')) == 0
    handles.Nt = 0;
    set(handles.Dfile, 'Enable', 'off');
    set(handles.Dfile, 'String', 'Select SPM M/EEG file');
    handles.Dfilename = [];
    set(handles.mstext, 'String', []);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in addcontrast.
function addcontrast_Callback(hObject, eventdata, handles)
% hObject    handle to addcontrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% averages only (add more options soon)
s = inputdlg({'from', 'to'}, 'Specify peri-stimulus time window for average [ms]');

from = str2num(s{1});
to = str2num(s{2});

strname = sprintf('Average from %d to %d ms', from, to);

% checks
if from < handles.ms(1)
    errordlg(sprintf('Start of time window must be later than %d ms.', handles.ms(1)));
    return;
end

if to > handles.ms(2)
    errordlg(sprintf('End of time window must be earlier than %d ms.', handles.ms(2)));
    return;
end

if from > to
    errordlg('Start of time window must be earlier than its end.');
    return;
end

% add contrast to list
set(handles.listbox2, 'String', [get(handles.listbox2, 'String'); {strname}]);

handles.con_types = cat(1, handles.con_types, {'average', from, to});

set(handles.compute, 'Enable', 'on');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in removecontrast.
function removecontrast_Callback(hObject, eventdata, handles)
% hObject    handle to removecontrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(get(handles.listbox2, 'String')) == 0, return, end

ind = get(handles.listbox2, 'Value');
tmp = get(handles.listbox2, 'String');
tmp(ind) = [];
handles.con_types(ind) = [];
set(handles.listbox2, 'Value', max(1, ind-1));
set(handles.listbox2, 'String', tmp);

% user removed all contrasts?
if length(get(handles.listbox2, 'String')) == 0
    handles.con_types = {};
    set(handles.compute, 'Enable', 'off');
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in compute.
function compute_Callback(hObject, eventdata, handles)
% hObject    handle to compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% make sure images are not flipped
global defaults
defaults.analyze.flip = 0;

con_types = handles.con_types;
D = handles.D;

ms = 1000/D.Radc*[-D.events.start:D.events.stop];

% assemble contrast vectors into matrix
for i = 1:size(handles.con_types, 1)
    c = zeros(handles.Nt, 1);
    
    if strcmp(con_types(i, 1), 'average')
        % window in time points
        [tmp, tp(1)] = min((ms-con_types{i, 2}).^2);
        [tmp, tp(2)] = min((ms-con_types{i, 3}).^2);

        c(tp(1):tp(2)) = 1./(tp(2)-tp(1)+1);
        C(:,i) = c;
    end
    
end
fnames = get(handles.listbox1, 'String');

spm('Pointer', 'Watch');drawnow

% compute and write contrast images
for j = 1:size(fnames, 1); % over files
    
    % map file    
    Vbeta = nifti(deblank(fnames{j}));  
    
    % cd to target directory
    [path, f] = fileparts(fnames{j});
    cd(path);
    
    for i = 1:size(C, 2) % over contrasts
              
       % code taken from spm_contrasts
       fprintf('\t%-32s: %-10s%20s', sprintf('file %s, contrast %d', fnames{j}, i),...
           '(spm_add)','...initialising') %-#

       %-Prepare handle for contrast image
       %-----------------------------------------------------------
       if strcmp(con_types{i, 1}, 'average')
          descrip = sprintf('SPM contrast - average from %d to %d ms',...
           con_types{i, 2}, con_types{i, 3});
       else
           descrip = sprintf('SPM contrast - %d: %s', i, con_type{i, 1});
       end
       
       % prepare nifti image (the usual spm_add doesn't seem to work for
       % many input files under windows)
       Vcon = nifti;
       Vcon.descrip = descrip;

       Dcon = file_array;
       Dcon.fname = sprintf('%s_con_%04d.img', f, i);
       Dcon.dtype = spm_type('float32');
       Dcon.offset  = ceil(348/8)*8;
       Dcon.dim = Vbeta.dat.dim(1:3);

       Vcon.dat = Dcon;
       
       %-Write image
       %-----------------------------------------------------------
       fprintf('%s%20s', repmat(sprintf('\b'),1,20),'...computing')%-#

       % replace loop by one line if possible
       d = zeros(Vbeta.dat.dim(1:3));
       for k = 1:Vbeta.dat.dim(4)
           d = d + Vbeta.dat(:,:,:, k)*C(k,i);
       end
       
       if Vbeta.dat.dim(3) == 1
           Dcon(:,:) = d;
       else
           Dcon.dat(:,:,:) = d;
       end
       
       Vcon.dat = Dcon;
       create(Vcon);
       

       fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),sprintf(...
           '...written %s',spm_str_manip(Vcon.dat.fname,'t')))%-#


   end
end

spm('Pointer', 'Arrow');

% --- Executes on button press in Dfile.
function Dfile_Callback(hObject, eventdata, handles)
% hObject    handle to Dfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname, filterindex] = uigetfile('*.mat', 'Select a SPM M/EEG file');

try
    D = spm_eeg_ldata(fullfile(pathname, filename));
catch
    errordlg('Selected file isn''t a valid SPM M/EEG file.')
    return;
end

try
    D.Nsamples;
catch
    errordlg('Selected file isn''t a valid SPM M/EEG file.')
    return;
end

if D.Nsamples ~= handles.Nt
    errordlg(sprintf('Number of time points (%d) different from images (%d).', D.Nsamples, handles.Nt))
    return;
end

% Dfile is valid
handles.Dfilename = fullfile(pathname, filename);
handles.D = D;
set(handles.listbox2, 'Enable', 'on');
set(handles.addcontrast, 'Enable', 'on');
set(handles.removecontrast, 'Enable', 'on');
set(handles.clearcontrasts, 'Enable', 'on');
handles.ms = 1000/D.Radc*[-D.events.start D.events.stop];
set(handles.mstext, 'Enable', 'on');
set(handles.mstext, 'String', ...
    sprintf('%d to %d ms, %d time points, %.2f ms steps', ...
    handles.ms(1), handles.ms(2), handles.Nt, 1000/D.Radc));
set(handles.Dfile, 'String', handles.Dfilename);
    
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in adddir.
function adddir_Callback(hObject, eventdata, handles)
% hObject    handle to adddir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d = uigetdir(pwd, 'Select a directory to search for SPM M/EEG-files');

if d == 0, return, end;

[flist] = spm_eeg_firstlevel_searchdir(d, []);

spm('Pointer', 'Watch'); drawnow;

% check for these files that length (time) is the same
for i = 1:size(flist,1)
    
    [status, Iimg, Limg] = spm_eeg_firstlevel_checkfile(flist(i,:), handles.Nt);

    % suppress error warning that file in list is of unequal length

    if status ~= 0
        % add image to list
        set(handles.listbox1, 'String', [get(handles.listbox1, 'String'); {deblank(flist(i, :))}]);
        set(handles.Dfile, 'Enable', 'on'); drawnow;
        handles.Nt = status;
    end
    
    drawnow;
end

spm('Pointer', 'Arrow');

% Update handles structure
guidata(hObject, handles);





% --- Executes on button press in clearfiles.
function clearfiles_Callback(hObject, eventdata, handles)
% hObject    handle to clearfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.listbox1, 'String', []);

handles.Nt = 0;
set(handles.Dfile, 'Enable', 'off');
set(handles.Dfile, 'String', 'Select SPM M/EEG file');
handles.Dfilename = [];
set(handles.mstext, 'String', []);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in clearcontrasts.
function clearcontrasts_Callback(hObject, eventdata, handles)
% hObject    handle to clearcontrasts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.listbox2, 'String', []);

handles.con_types = {};
set(handles.compute, 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);


