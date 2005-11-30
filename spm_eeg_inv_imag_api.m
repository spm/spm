function varargout = spm_eeg_inv_imag_api(varargin)
% SPM_EEG_INV_IMAG_API M-file for spm_eeg_inv_imag_api.fig
%    FIG = SPM_EEG_INV_IMAG_API launch spm_eeg_inv_imag_api GUI.
%    D = SPM_EEG_INV_IMAG_API(D,newInv) launch spm_eeg_inv_imag_api GUI to
%    analyse data D (default newInv = 1 i.e. new inverse analysis).
%    SPM_EEG_INV_IMAG_API('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 21-Nov-2005 16:33:57

spm_defaults
spm('Clear')

if ~nargin  % LAUNCH GUI

	fig     = openfig(mfilename,'new');
    WS      = spm('WinScale');
    Rect    = spm('WinSize','Menu','raw').*WS;
    set(fig,'units','pixels');
    Fdim    = get(fig,'position');
    set(fig,'position',[Rect(1) Rect(2) Fdim(3) Fdim(4)]);
    handles = guihandles(fig);

    % Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    S = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D = spm_eeg_ldata(S);
    handles.D = D;
    set(handles.DataFile,'String',D.fname);
    [pth,nam,ext] = fileparts(S);
    cd(pth);
    
    guidata(fig,handles);

    spm_eeg_inv_imag_api_OpeningFcn(fig, [], handles);
       
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if nargout
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
	    error(sprintf('Wrong input format\n'));
	end
    
else
    
    error(sprintf('Wrong input format\n'));

end


% Last Modified by GUIDE v2.5 20-Nov-2005 16:17:15

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


% --- Executes just before spm_eeg_inv_imag_api is made visible.
function spm_eeg_inv_imag_api_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_eeg_inv_imag_api (see VARARGIN)

% Choose default command line output for spm_eeg_inv_imag_api
handles.output = hObject;

try
    D         = handles.D;
catch
    S         = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D         = spm_eeg_ldata(S);
    handles.D = D;
end

if ~isfield(D,'inv')
    set(handles.Analysis,'Value',1);
end

cd(D.path);
set(handles.Imaging,'Value',1);
set(handles.DataFile,'String',D.fname);

guidata(hObject,handles);

% UIWAIT makes spm_eeg_inv_imag_api wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spm_eeg_inv_imag_api_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;

try
    D = get(handles.figure1,'UserData');
    if str2num(version('-release'))>=14
        save(fullfile(D.path, D.fname), '-V6', 'D');
    else
      	save(fullfile(D.path, D.fname), 'D');
    end
    varargout{1} = D;
end


% --- Executes on button press in radiobutton Analysis.
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton 'Analysis' (see GCBO)
NewInv = get(handles.Analysis,'Value');
try
    D         = handles.D;
catch
    spm_eeg_inv_imag_api_OpeningFcn(hObject, eventdata, handles);
end
if ~NewInv & ~isfield(D,'inv')
    set(handles.Analysis,'Value',1);
end

guidata(hObject,handles);


% --- Executes on button press in togglebutton Imaging.
function Imaging_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton 'Imaging' (see GCBO)
set(handles.ECD,'Value',0);

% --- Executes on button press in togglebutton ECD.
function ECD_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton 'ECD' (see GCBO)
set(handles.Imaging,'Value',0);

% --- Executes on button press in CreateMeshes.
function CreateMeshes_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton CreateMeshes (see GCBO)

NewInv = get(handles.Analysis,'Value');
Meth   = get(handles.Imaging,'Value');
D      = handles.D;

if ~isfield(D,'inv')
    val = 1;
else
    val = length(D.inv) + NewInv;
end

load('defaults_eeg_inv.mat') % load default structure invdetfs

if Meth
    D.inv{val} = invdefts.imag;
else
    D.inv{val} = invdefts.ECD;
end

if NewInv
    D = set_CommentDate(D);
end

D = spm_eeg_inv_mesh_ui(D);

handles.D = D;
set(handles.Analysis,'Value',0);

guidata(hObject,handles);


% --- Executes on button press in Data Reg.
function DataReg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton DataReg (see GCBO)

NewInv = get(handles.Analysis,'Value');

D      = handles.D;
if ~isfield(D,'inv')
    error(sprintf('You need to compute the Meshes first\n'));
end

val    = length(D.inv);

if NewInv
    val = val + 1;
    load('defaults_eeg_inv.mat') % load default structure invdetfs
    Meth   = get(handles.Imaging,'Value');
    if Meth
        D.inv{val} = invdefts.imag;
    else
        D.inv{val} = invdefts.ECD;
    end
    D = set_CommentDate(D);
end
    
if val > 1
    names_m = fieldnames(D.inv{val-1}.mesh);
    count_m = 0;
    for i = 1:length(names_m)
        count_m = count_m + isempty(getfield(D.inv{val}.mesh,names_m{i}));
    end
    if count_m == length(names_m)
        disp('Copy Mesh files from previous analysis');
        D = spm_eeg_inv_copyfields(D,[2 0 0 0]);
    end
end
                 
D = spm_eeg_inv_datareg_ui(D);

handles.D = D;
set(handles.Analysis,'Value',0);

guidata(hObject,handles);


% --- Executes on button press in Forward Sol.
function Forward_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton Inverse (see GCBO)
NewInv = get(handles.Analysis,'Value');
Meth   = get(handles.Imaging,'Value');

D      = handles.D;
if ~isfield(D,'inv')
    error(sprintf('You need to compute Meshes and Registrate the Data first\n'));
end

val    = length(D.inv);

if NewInv
    val = val + 1;
    load('defaults_eeg_inv.mat') % load default structure invdetfs
    if Meth
        D.inv{val} = invdefts.imag;
    else
        D.inv{val} = invdefts.ECD;
    end
    D = set_CommentDate(D);
end

if val > 1
    names_m = fieldnames(D.inv{val-1}.mesh);
    count_m = 0;
    for i = 1:length(names_m)
        count_m = count_m + isempty(getfield(D.inv{val}.mesh,names_m{i}));
    end
    if count_m == length(names_m)
        disp('Copy Mesh files from previous analysis');
        D = spm_eeg_inv_copyfields(D,[2 0 0 0]);
    end
    names_d = fieldnames(D.inv{val-1}.datareg);
    count_d = 0;
    for i = 1:length(names_d)
        count_d = count_d + isempty(getfield(D.inv{val}.datareg,names_d{i}));
    end
    if count_d == length(names_d)
        disp('Copy Datareg files from previous analysis');
        D = spm_eeg_inv_copyfields(D,[0 2 0 0]);
    end
end
                  
if Meth
    D = spm_eeg_inv_forward_ui(D);
else
    D = spm_eeg_inv_fp_elec_Rsph_ui(D); 
end

handles.D = D;
set(handles.Analysis,'Value',0);

guidata(hObject,handles);


% --- Executes on button press in CreateMeshes.
function Inverse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton CreateMeshes (see GCBO)
NewInv = get(handles.Analysis,'Value');
Meth   = get(handles.Imaging,'Value');

D      = handles.D;
if ~isfield(D,'inv')
    error(sprintf('You need to compute the whole generative model first (Meshes/DataReg/Forward)\n'));
end

val    = length(D.inv);

if NewInv
    val = val + 1;
    load('defaults_eeg_inv.mat') % load default structure invdetfs
    if Meth
        D.inv{val} = invdefts.imag;
    else
        D.inv{val} = invdefts.ECD;
    end
    D = set_CommentDate(D);
end

if val > 1
    names_m = fieldnames(D.inv{val-1}.mesh);
    count_m = 0;
    for i = 1:length(names_m)
        count_m = count_m + isempty(getfield(D.inv{val}.mesh,names_m{i}));
    end
    if count_m == length(names_m)
        disp('Copy Mesh files from previous analysis');
        D = spm_eeg_inv_copyfields(D,[2 0 0 0]);
    end
    names_d = fieldnames(D.inv{val-1}.datareg);
    count_d = 0;
    for i = 1:length(names_d)
        count_d = count_d + isempty(getfield(D.inv{val}.datareg,names_d{i}));
    end
    if count_d == length(names_d)
        disp('Copy Datareg files from previous analysis');
        D = spm_eeg_inv_copyfields(D,[0 2 0 0]);
    end
    names_f = fieldnames(D.inv{1}.forward);
    count_f = 0;
    for i = 1:length(names_f)
        count_f = count_f + isempty(getfield(D.inv{val}.forward,names_f{i}));
    end
    if count_f == length(names_f)
        disp('Copy Forward files from previous analysis');
        D = spm_eeg_inv_copyfields(D,[0 0 2 0]);
    end
end

if Meth
    D = spm_eeg_inv_inverse_ui(D);
else
    D = spm_eeg_inv_ecd_ui(D);
end

handles.D = D;
set(handles.Analysis,'Value',0);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Other_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu 'Other' (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu1.
function Other_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu 'Other' (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

contents = get(hObject,'String');
task     = contents{get(hObject,'Value')};
D        = handles.D;

switch task
    case 'check meshes'
        spm_eeg_inv_checkmeshes(D);
    case 'read polhemus'
        spm_eeg_inv_ReadPolhemus(D);
    case 'check datareg'
        spm_eeg_inv_checkdatareg(D);
    case 'visualisation'
        return
    case 'delete analysis'
        spm_eeg_inv_deletefields(D);
end

handles.D = D;

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function DataFile_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function DataFile_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of DataFile as text
%        str2double(get(hObject,'String')) returns contents of DataFile as a double
try
    S = get(handles.DataFile,'String');
    [pth,nam,ext] = fileparts(S);
    if ~isempty(pth)
        cd(pth)
    end
    D = spm_eeg_ldata(S);
    handles.D = D;        
    guidata(hObject,handles);
catch
    Load_Callback(hObject, eventdata, handles);
end
   

% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

S         = spm_select(1, '.mat', 'Select EEG/MEG mat file');
D         = spm_eeg_ldata(S);
handles.D = D;
set(handles.DataFile,'String',D.fname);
[pth,nam,ext] = fileparts(S);
cd(pth);
guidata(hObject,handles);


% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spm_eeg_inv_imag_api_OutputFcn(hObject, eventdata, handles);
close(handles.figure1);


% --- Set Comment and Date for new inverse analysis
function S = set_CommentDate(D)

val = length(D.inv);

clck = fix(clock);
if clck(5) < 10
    clck = [num2str(clck(4)) ':0' num2str(clck(5))];
else
    clck = [num2str(clck(4)) ':' num2str(clck(5))];
end
D.inv{val}.date = strvcat(date,clck);
D.inv{val}.comment = inputdlg('Comment/Label for this analysis:');

S = D;


% --------------------------------------------------------------------
% --- FRAMES
function bigframe1_Callback(hObject, eventdata, handles)
function bigframe2_Callback(hObject, eventdata, handles)
function bigframe3_Callback(hObject, eventdata, handles)


