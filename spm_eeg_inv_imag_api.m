function varargout = spm_eeg_inv_imag_api(varargin)
% API for EEG/MEG source reconstruction interface
% FORMAT:
%    FIG = SPM_EEG_INV_IMAG_API launch spm_eeg_inv_imag_api GUI.
%    SPM_EEG_INV_IMAG_API('callback_name', ...) invoke the named callback.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_imag_api.m 4260 2011-03-23 13:42:21Z vladimir $

spm('Clear');

% Launch API
%==========================================================================
if nargin < 2 
    
    % open figure
    %----------------------------------------------------------------------
    fig     = openfig(mfilename,'reuse');
    set(fig,'name',[spm('ver') ': ' get(fig,'name')]);
    Rect    = spm('WinSize','Menu');
    S0      = spm('WinSize','0',1);
    set(fig,'units','pixels');
    Fdim    = get(fig,'position');
    set(fig,'position',[S0(1)+Rect(1) S0(2)+Rect(2) Fdim(3) Fdim(4)]);
    handles = guihandles(fig);

    % Use system color scheme for figure:
    %----------------------------------------------------------------------
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    handles.fig = fig;
    guidata(fig,handles);
    
    % intialise with D
    %----------------------------------------------------------------------
    try
        D = spm_eeg_inv_check(varargin{1});
        set(handles.DataFile,'String',D.fname);
        set(handles.Exit,'enable','on')
        cd(D.path);
        handles.D = D;
        Reset(fig, [], handles);
        guidata(fig,handles);
    end
    
% INVOKE NAMED SUBFUNCTION OR CALLBACK
%--------------------------------------------------------------------------
elseif ischar(varargin{1}) 
    if nargout
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    else
        feval(varargin{:}); % FEVAL switchyard
    end

else
    error('Wrong input format.');
end

% MAIN FUNCTIONS FOR MODEL SEPCIFICATION AND INVERSION
%==========================================================================

% --- Executes on button press in CreateMeshes.
%--------------------------------------------------------------------------
function CreateMeshes_Callback(hObject, eventdata, handles)
handles.D = spm_eeg_inv_mesh_ui(handles.D, handles.D.val, 0);
Reset(hObject, eventdata, handles);


% --- Executes on button press in Reg2tem.
%--------------------------------------------------------------------------
function Reg2tem_Callback(hObject, eventdata, handles)
handles.D = spm_eeg_inv_mesh_ui(handles.D, handles.D.val, 1);
Reset(hObject, eventdata, handles);


% --- Executes on button press in Data Reg.
%--------------------------------------------------------------------------
function DataReg_Callback(hObject, eventdata, handles)
handles.D = spm_eeg_inv_datareg_ui(handles.D);
Reset(hObject, eventdata, handles);


% --- Executes on button press in Forward Model.
%--------------------------------------------------------------------------
function Forward_Callback(hObject, eventdata, handles)
handles.D = spm_eeg_inv_forward_ui(handles.D);
Reset(hObject, eventdata, handles);


% --- Executes on button press in Invert.
%--------------------------------------------------------------------------
function Inverse_Callback(hObject, eventdata, handles)
handles.D = spm_eeg_invert_ui(handles.D);
Reset(hObject, eventdata, handles);


% --- Executes on button press in contrast.
%--------------------------------------------------------------------------
function contrast_Callback(hObject, eventdata, handles)
handles.D = spm_eeg_inv_results_ui(handles.D);
Reset(hObject, eventdata, handles);


% --- Executes on button press in Image.
%--------------------------------------------------------------------------
function Image_Callback(hObject, eventdata,handles)
handles.D.inv{handles.D.val}.contrast.display = 1;
handles.D = spm_eeg_inv_Mesh2Voxels(handles.D);
Reset(hObject, eventdata, handles);



% LOAD AND EXIT
%==========================================================================

% --- Executes on button press in Load.
%--------------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
[S, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
if ~sts, return; end
D = spm_eeg_load(S);

[ok, D] = check(D, 'sensfid');

if ~ok
    if check(D, 'basic')
        warndlg(['The requested file is not ready for source reconstruction.'...
            'See Matlab window for details.']);
    else
        warndlg('The meeg file is corrupt or incomplete');
    end
    return
end

set(handles.DataFile,'String',D.fname);
set(handles.Exit,'enable','on');
cd(D.path);
handles.D     = D;
Reset(hObject, eventdata, handles);


% --- Executes on button press in Exit.
%--------------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
D             = handles.D;
D.save;
varargout{1} = handles.D;
assignin('base','D',handles.D)


% FUCNTIONS FOR MANAGING DIFFERENT MODELS
%==========================================================================

% --- Executes on button press in new.
%--------------------------------------------------------------------------
function new_Callback(hObject, eventdata, handles)
D  = handles.D;
if ~isfield(D,'inv')
    val   = 1;
elseif isempty(D.inv)
    val   = 1;
else
    val        = length(D.inv) + 1;
    D.inv{val} = D.inv{D.val};
end

% set D in handles and update analysis specific buttons
%--------------------------------------------------------------------------
D.val     = val;
D         = set_CommentDate(D);
handles.D = D;
set(handles.CreateMeshes,'enable','on')
set(handles.Reg2tem,'enable','on')
Reset(hObject, eventdata, handles);

% --- Executes on button press in next.
%--------------------------------------------------------------------------
function next_Callback(hObject, eventdata, handles)
if handles.D.val < length(handles.D.inv)
    handles.D.val = handles.D.val + 1;
end
Reset(hObject, eventdata, handles);


% --- Executes on button press in previous.
%--------------------------------------------------------------------------
function previous_Callback(hObject, eventdata, handles)
if handles.D.val > 1
    handles.D.val = handles.D.val - 1;
end
Reset(hObject, eventdata, handles);


% --- Executes on button press in clear.
%--------------------------------------------------------------------------
function clear_Callback(hObject, eventdata, handles)
try
    inv.comment = handles.D.inv{handles.D.val}.comment;
    inv.date    = handles.D.inv{handles.D.val}.date;
    handles.D.inv{handles.D.val} = inv;
end
Reset(hObject, eventdata, handles);

% --- Executes on button press in delete.
%--------------------------------------------------------------------------
function delete_Callback(hObject, eventdata, handles)
if ~isempty(handles.D.inv)
    try
        str = handles.D.inv{handles.D.val}.comment;
        warndlg({'you are about to delete:',str{1}});
        uiwait
    end
    handles.D.inv(handles.D.val) = [];
    handles.D.val                = handles.D.val - 1;
end
Reset(hObject, eventdata, handles);



% Auxillary functions
%==========================================================================
function Reset(hObject, eventdata, handles)

% Check to see if a new analysis is required
%--------------------------------------------------------------------------
try
    set(handles.DataFile,'String',handles.D.fname);
end
if ~isfield(handles.D,'inv')
    new_Callback(hObject, eventdata, handles)
    return
end
if isempty(handles.D.inv)
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
Q  = handles.D.inv{val};

% === This is for backward compatibility with SPM8b. Can be removed after
% some time
if isfield(Q, 'mesh') &&...
        isfield(Q.mesh, 'tess_ctx') && ~isa(Q.mesh.tess_ctx, 'char')

    warning(['This is an old version of SPM8b inversion. ',...
        'You can only review and export solutions. ',...
        'Clear and invert again to update']);
    Q = rmfield(Q, {'mesh', 'datareg', 'forward'});
end
% =========================================================================

set(handles.new,      'enable','on','value',0)
set(handles.clear,    'enable','on','value',0)
set(handles.delete,   'enable','on','value',0)
set(handles.next,     'value',0)
set(handles.previous, 'value',0)

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
    try
        str = sprintf('%i: %s',val,Q.comment);
    catch
        str = sprintf('%i',val);
    end
end
set(handles.val, 'Value',val,'string',str);

% condition specification
%--------------------------------------------------------------------------
try
    handles.D.con = max(handles.D.con,1);
    if handles.D.con > length(handles.D.inv{val}.inverse.J);
       handles.D.con = 1;
    end
catch
    try 
       handles.D.con = length(handles.D.inv{val}.inverse.J);
    catch
       handles.D.con = 0;
    end
end
if handles.D.con
    str = sprintf('condition %d',handles.D.con);
    set(handles.con,'String',str,'Enable','on','Value',0)
else
    set(handles.con,'Enable','off','Value',0)
end


% check anaylsis buttons
%--------------------------------------------------------------------------
set(handles.DataReg, 'enable','off')
set(handles.Forward, 'enable','off')
set(handles.Inverse, 'enable','off')
set(handles.contrast,'enable','off')
set(handles.Image,   'enable','off')

set(handles.CheckReg,     'enable','off','Value',0)
set(handles.CheckMesh,    'enable','off','Value',0)
set(handles.CheckForward, 'enable','off','Value',0)
set(handles.CheckInverse, 'enable','off','Value',0)
set(handles.CheckContrast,'enable','off','Value',0)
set(handles.CheckImage,   'enable','off','Value',0)
set(handles.Movie,        'enable','off','Value',0)
set(handles.Vis3D,        'enable','off','Value',0)
set(handles.Image,        'enable','off','Value',0)

set(handles.CreateMeshes,'enable','on')
set(handles.Reg2tem,'enable','on')
if isfield(Q, 'mesh')
    set(handles.DataReg,  'enable','on')
    set(handles.CheckMesh,'enable','on')
    if isfield(Q,'datareg') && isfield(Q.datareg, 'sensors')
        set(handles.Forward, 'enable','on')
        set(handles.CheckReg,'enable','on')
        if isfield(Q,'forward') && isfield(Q.forward, 'vol')
            set(handles.Inverse,     'enable','on')
            set(handles.CheckForward,'enable','on')
        end
    end
end
if isfield(Q,'inverse') && isfield(Q, 'method')
    set(handles.CheckInverse,'enable','on')
    if isfield(Q.inverse,'J')
        set(handles.contrast,    'enable','on')
        set(handles.Movie,       'enable','on')
        set(handles.Vis3D,       'enable','on')
        if isfield(Q,'contrast')
            set(handles.CheckContrast,'enable','on')
            set(handles.Image,        'enable','on')
            if isfield(Q.contrast,'fname')
                set(handles.CheckImage,'enable','on')
            end
        end
    end
end

try
    if strcmp(handles.D.inv{handles.D.val}.method,'Imaging')
        set(handles.CheckInverse,'String','mip');
        set(handles.PST,'Enable','on');
    else
        set(handles.CheckInverse,'String','dip');
        set(handles.PST,'Enable','off');
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


% CHECKS AND DISPLAYS
%==========================================================================

% --- Executes on button press in CheckMesh.
%--------------------------------------------------------------------------
function CheckMesh_Callback(hObject, eventdata, handles)
spm_eeg_inv_checkmeshes(handles.D);
Reset(hObject, eventdata, handles);

% --- Executes on button press in CheckReg.
%--------------------------------------------------------------------------
function CheckReg_Callback(hObject, eventdata, handles)
% check and display registration
%--------------------------------------------------------------------------
spm_eeg_inv_checkdatareg(handles.D);

Reset(hObject, eventdata, handles);

% --- Executes on button press in CheckForward.
%--------------------------------------------------------------------------
function CheckForward_Callback(hObject, eventdata, handles)
spm_eeg_inv_checkforward(handles.D);
Reset(hObject, eventdata, handles);

% --- Executes on button press in CheckInverse.
%--------------------------------------------------------------------------
function CheckInverse_Callback(hObject, eventdata, handles)
if strcmp(handles.D.inv{handles.D.val}.method,'Imaging')
    PST    = str2num(get(handles.PST,'String'));
    spm_eeg_invert_display(handles.D,PST);
    if length(PST) == 3 && get(handles.extract, 'Value')
        handles.D = spm_eeg_inv_extract_ui(handles.D, handles.D.val, PST);        
    end
elseif strcmp(handles.D.inv{handles.D.val}.method, 'vbecd')
    spm_eeg_inv_vbecd_disp('init',handles.D);
end
Reset(hObject, eventdata, handles);

% --- Executes on button press in Movie.
%--------------------------------------------------------------------------
function Movie_Callback(hObject, eventdata, handles)
figure(spm_figure('GetWin','Graphics'));
PST(1) = str2num(get(handles.Start,'String'));
PST(2) = str2num(get(handles.Stop ,'String'));
spm_eeg_invert_display(handles.D,PST);
Reset(hObject, eventdata, handles);

% --- Executes on button press in CheckContrast.
%--------------------------------------------------------------------------
function CheckContrast_Callback(hObject, eventdata, handles)
spm_eeg_inv_results_display(handles.D);
Reset(hObject, eventdata, handles);

% --- Executes on button press in Vis3D.
%--------------------------------------------------------------------------
function Vis3D_Callback(hObject, eventdata, handles)
Exit_Callback(hObject, eventdata, handles)
try
    spm_eeg_inv_visu3D_api(handles.D);
catch
    spm_eeg_review(handles.D,6,handles.D.val);
end
Reset(hObject, eventdata, handles);

% --- Executes on button press in CheckImage.
%--------------------------------------------------------------------------
function CheckImage_Callback(hObject, eventdata, handles)
spm_eeg_inv_image_display(handles.D)
Reset(hObject, eventdata, handles);

% --- Executes on button press in con.
%--------------------------------------------------------------------------
function con_Callback(hObject, eventdata, handles)
try
    handles.D.con = handles.D.con + 1;
    if handles.D.con > length(handles.D.inverse.J);
        handles.D.con = 1;
    end
end
Reset(hObject, eventdata, handles);

% --- Executes on button press in help.
%--------------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
edit spm_eeg_inv_help


% --- Executes on button press in group.
%--------------------------------------------------------------------------
function group_Callback(hObject, eventdata, handles)
spm_eeg_inv_group;
