function varargout = spm_eeg_inv_imag_api(varargin)
% SPM_EEG_INV_IMAG_API M-file for spm_eeg_inv_imag_api.fig
%    FIG = SPM_EEG_INV_IMAG_API launch spm_eeg_inv_imag_api GUI.
%    D = SPM_EEG_INV_IMAG_API(D,newInv) launch spm_eeg_inv_imag_api GUI to
%    analyse data D (default newInv = 1 i.e. new inverse analysis).
%    SPM_EEG_INV_IMAG_API('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 08-Sep-2006 15:55:36
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_imag_api.m 628 2006-09-18 13:33:26Z karl $

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
    try
        D.modality;
    catch
        D.modality = 'EEG';
        warndlg('Assuming data is EEG')
    end
    handles.D = D;
    set(handles.DataFile,'String',D.fname);
    [pth,nam,ext] = fileparts(S);
    cd(pth);
    
    handles.fig = fig;
    guidata(fig,handles);
    spm_eeg_inv_imag_api_OpeningFcn(fig, [], handles);

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    if nargout
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    else
        feval(varargin{:}); % FEVAL switchyard
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
Reset(hObject, eventdata, handles)




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
    if spm_matlab_version_chk('7.1') >= 0
        save(fullfile(D.path, D.fname), '-V6', 'D');
    else
        save(fullfile(D.path, D.fname), 'D');
    end
    varargout{1} = D;
    assignin('base','D',handles.D)
end


% --- Executes on button press in CreateMeshes.
function CreateMeshes_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton CreateMeshes (see GCBO)

handles.D = spm_eeg_inv_mesh_ui(handles.D);
Reset(hObject, eventdata, handles);

% --- Executes on button press in For2tem.
function Reg2tem_Callback(hObject, eventdata, handles)
% hObject    handle to For2tem (see GCBO)

global defaults
Cdir = [defaults.SWD filesep 'EEGcanonical'];
D    = handles.D;
val  = D.val;

D.inv{val}.mesh.sMRI         = fullfile(Cdir,'smri.img,1');
D.inv{val}.mesh.def          = fullfile(Cdir,'smri_vbm_sn_1.mat');
D.inv{val}.mesh.invdef       = fullfile(Cdir,'smri_vbm_inv_sn_1.mat');
D.inv{val}.mesh.msk_iskull   = fullfile(Cdir,'smri_iskull.img');
D.inv{val}.mesh.msk_scalp    = fullfile(Cdir,'smri_scalp.img');
D.inv{val}.mesh.tess_ctx     = fullfile(Cdir,'smri_CortexMesh_3004.mat');
D.inv{val}.mesh.tess_iskull  = fullfile(Cdir,'smri_iskull_Mesh_2002.mat');
D.inv{val}.mesh.tess_scalp   = fullfile(Cdir,'smri_scalp_Mesh_2002.mat');
D.inv{val}.mesh.CtxGeoDist   = fullfile(Cdir,'smri_CortexGeoDist_3004.mat');
D.inv{val}.mesh.Ctx_Nv       = 3004;
D.inv{val}.mesh.Ctx_Nf       = 6000;
D.inv{val}.mesh.Iskull_Nv    = 3004;
D.inv{val}.mesh.Iskull_Nf    = 4000;
D.inv{val}.mesh.Scalp_Nv     = 3004;
D.inv{val}.mesh.Scalp_Nf     = 4000;

D.inv{val}.datareg.scalpvert = D.inv{val}.mesh.tess_scalp;
D.inv{val}.datareg.fidmri    = fullfile(Cdir,'smri_fids.mat');

handles.D = D;
DataReg_Callback(hObject, eventdata, handles)

% --- Executes on button press in Data Reg.
function DataReg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton DataReg (see GCBO)

handles.D = spm_eeg_inv_datareg_ui(handles.D);
Reset(hObject, eventdata, handles);


% --- Executes on button press in Forward Model.
function Forward_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton Inverse (see GCBO)

method = questdlg('recontruction','Please select','Imaging','ECD','Imaging')
handles.D.inv{handles.D.val}.method = method;

if strcmp(method,'Imaging')
    handles.D = spm_eeg_inv_forward_ui(handles.D);
else
    handles.D = spm_eeg_inv_elec_Rsph_ui(handles.D);
end
Reset(hObject, eventdata, handles);

% --- Executes on button press in Reg2tem.
%--------------------------------------------------------------------------
function For2tem_Callback(hObject, eventdata, handles)
% hObject    handle to Reg2tem (see GCBO)

% get standard head model
%--------------------------------------------------------------------------
global defaults
load(fullfile(defaults.SWD,'EEGcanonical','standard_head_model.mat'))

% Input standard electrode sets
%--------------------------------------------------------------------------
[set_Nel,set_name]   = spm_eeg_inv_electrset;
el_set               = spm_input('Which set of electrodes','+1','m',set_name);
[el_sphc,el_name]    = spm_eeg_inv_electrset(el_set);
flags_el.q_RealLoc   = 0;
flags_el.q_RealNI    = 0;
[electrodes,f_el]    = spm_eeg_inv_model('Elec2Scalp',model.head(end), ...
                       el_sphc,el_name,flags_el);
model.electrodes     = electrodes;
model.flags.fl_elec  = f_el;

% [re]-set the sphere model
%--------------------------------------------------------------------------
try, model           = rmfield(model,'spheres'); end
[spheres,d,L,q,f_Rs] = spm_eeg_inv_Rsph(model,[]);
model.spheres        = spheres;
model.flags.fl_RealS = f_Rs;

% set in D
%--------------------------------------------------------------------------
D = handles.D;
D.inv{D.val}.forward     = model;
D.inv{D.val}.mesh.sMRI   = fullfile(defaults.SWD,'EEGcanonical','smri.img');
D.inv{D.val}.mesh.invdef = fullfile(defaults.SWD,'EEGcanonical','smri_vbm_inv_sn.mat');

handles.D = D;
Forward_Callback(hObject, eventdata, handles);


% --- Executes on button press in Invert.
function Inverse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton CreateMeshes (see GCBO)

if strcmp(handles.D.inv{handles.D.val}.method,'Imaging')
    handles.D = spm_eeg_inv_inverse_ui(handles.D);
else
    handles.D = spm_eeg_inv_ecd_ui(handles.D);
end
Reset(hObject, eventdata, handles);


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

        if strcmp(D.inv{D.val}.method,'ECD')
            sMRI   = D.inv{D.val}.mesh.sMRI;
            resdip = D.inv{D.val}.inverse.resdip;
            spm_eeg_inv_ecd_DrawDip('Init',resdip,sMRI)
            assignin('base','resdip',resdip);
        else
            try
                load(D.inv{D.val}.inverse.resfile);
                AvJev = abs(mean(J,2));
                load(D.inv{D.val}.mesh.tess_ctx);
                
                % Temporary Visualization
                [Finter,Fgraph] = spm('FnUIsetup','Stats: Results');
                figure(Fgraph);
                subplot(3,1,2)

                patch('Vertices',vert,'Faces',face,'FaceVertexCData',AvJev,'FaceColor','flat');
                view(-90,0);
                shading interp
                colormap('pink')
                axis image
                colorbar
                title('Averaged activity over the whole time window');
            catch
                warndlg('please invert model')
            end
            
            return
            warndlg('under construction')
            spm_eeg_inv_visu3D(D);
        end
end


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
    Reset(hObject, eventdata, handles);
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
Reset(hObject, eventdata, handles);

% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spm_eeg_inv_imag_api_OutputFcn(hObject, eventdata, handles);
close(handles.figure1);


% --- Executes on button press in new.
function new_Callback(hObject, eventdata, handles)
% hObject    handle to new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

D  = handles.D;

if ~isfield(D,'inv')
    val   = 1;
    D.val = 1;
else
    val   = length(D.inv) + 1;
end

% if this is the first analysis, or method has changed, get defaults
%--------------------------------------------------------------------------
load('defaults_eeg_inv.mat')
if val == 1
    D.inv      = {invdefts.imag};
else
    D.inv{val} = D.inv{D.val};
    D.inv{val}.inverse = invdefts.imag.inverse;
end

% set D in handles and update analysis specific buttons
%--------------------------------------------------------------------------
D.val     = val;
D         = set_CommentDate(D);
handles.D = D;
Reset(hObject, eventdata, handles);


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.D.val < length(handles.D.inv)
    handles.D.val = handles.D.val + 1;
end
Reset(hObject, eventdata, handles);


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.D.val > 1
    handles.D.val = handles.D.val - 1;
end
Reset(hObject, eventdata, handles);

% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(handles.D.inv)
    str = handles.D.inv{handles.D.val}.comment;
    warndlg({'you are about to delete:',str{1}});
    uiwait
    handles.D.inv(handles.D.val) = [];
    handles.D.val                = handles.D.val - 1;
end
Reset(hObject, eventdata, handles);



% Auxillary functions
%==========================================================================
function Reset(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check to see if a new analysis is required
%--------------------------------------------------------------------------
if ~isfield(handles.D,'inv')
    new_Callback(hObject, eventdata, handles)
    return
end
if ~length(handles.D.inv)
    new_Callback(hObject, eventdata, handles)
    return
end
try
    val = handles.D.val;
    handles.D.inv{val};
catch
    handles.D.val = 1;
    val           = 1;
end

% analysis specification buttons
%--------------------------------------------------------------------------
Q   = handles.D.inv{val};
set(handles.new,      'enable','on')
set(handles.new,      'value',0)
set(handles.next,     'value',0)
set(handles.previous, 'value',0)
set(handles.delete,   'value',0)
if val < length(handles.D.inv)
    set(handles.next,    'enable','on')
end
if val > 1
    set(handles.previous,'enable','on')
end
if val == 1
    set(handles.previous,'enable','off')
end
if val == length(handles.D.inv)
    set(handles.next,    'enable','off')
end
try
    str = sprintf('%i: %s',val,Q.comment{1});
catch
    str = sprintf('%i',val);
end
set(handles.val, 'Value',val,'string',str);



% check anaylsis buttons
%--------------------------------------------------------------------------
set(handles.DataReg,'enable','off')
set(handles.Forward,'enable','off')
set(handles.Inverse,'enable','off')

if length(Q.mesh.sMRI)
    set(handles.DataReg,'enable','on')
    if length(Q.datareg.eeg2mri)
        set(handles.Forward,'enable','on')
           if length(spm_vec(Q.forward))
               set(handles.Inverse,'enable','on')
           end
    end
end

set(handles.fig,'Pointer','arrow')
assignin('base','D',handles.D)
guidata(hObject,handles);


% Set Comment and Date for new inverse analysis
%--------------------------------------------------------------------------
function S = set_CommentDate(D)

clck = fix(clock);
if clck(5) < 10
    clck = [num2str(clck(4)) ':0' num2str(clck(5))];
else
    clck = [num2str(clck(4)) ':' num2str(clck(5))];
end
D.inv{D.val}.date    = strvcat(date,clck);
D.inv{D.val}.comment = inputdlg('Comment/Label for this analysis:');
S = D;

% FRAMES
%--------------------------------------------------------------------------
function bigframe1_Callback(hObject, eventdata, handles)
function bigframe2_Callback(hObject, eventdata, handles)
function bigframe3_Callback(hObject, eventdata, handles)




