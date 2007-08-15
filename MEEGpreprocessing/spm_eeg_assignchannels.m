function varargout = spm_eeg_assignchannels(varargin)
% SPM_EEG_ASSIGNCHANNELS M-file for spm_eeg_assignchannels.fig
%      SPM_EEG_ASSIGNCHANNELS, by itself, creates a add SPM_EEG_ASSIGNCHANNELS or raises the existing
%      singleton*.
%
%      H = SPM_EEG_ASSIGNCHANNELS returns the handle to a add SPM_EEG_ASSIGNCHANNELS or the handle to
%      the existing singleton*.
%
%      SPM_EEG_ASSIGNCHANNELS('CALLBACK',hObject,eventData,handles,...)
%      calls the local
%      function named CALLBACK in SPM_EEG_ASSIGNCHANNELS.M with the given input arguments.
%
%      SPM_EEG_ASSIGNCHANNELS('Property','Value',...) creates a add
%      SPM_EEG_ASSIGNCHANNELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spm_eeg_assignchannels_OpeningFunction
%      gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to spm_eeg_assignchannels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spm_eeg_assignchannels

% Last Modified by GUIDE v2.5 13-Aug-2007 13:32:30


if str2num(version('-release'))<14
    errordlg('This tool only work on Matlab 7 and above');
    return;
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @spm_eeg_assignchannels_OpeningFcn, ...
    'gui_OutputFcn',  @spm_eeg_assignchannels_OutputFcn, ...
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

%------------------------------------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = spm_eeg_assignchannels_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
% close(handles.figure1);

%============================ Callback functions ===================================

% ------------------------ Channel type assignment ---------------------------------

% --- Executes on button press in loadDataFile.
function loadDataFile_Callback(hObject, eventdata, handles)
% hObject    handle to loadDataFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


P = spm_select(1, '\.mat$', 'Select EEG mat file');

if isempty(P)
    return;
end

D = spm_eeg_ldata(P);

handles.D = D;
handles.names = D.channels.name;
handles.Nnames = length(handles.names);


% populate list
if isfield(D.channels, 'types')
    handles.types=D.channels.types;
elseif isfield(D.channels, 'eeg')
    if strcmp(handles.D.modality, 'MEG')
        handles = set_types(handles, 'MEG', num2str(D.channels.eeg(:)'));
    else
        handles = set_types(handles, 'EEG', num2str(D.channels.eeg(:)'));
    end
else
    handles=set_types_to_default(handles);
end

handles=update_list(handles);

guidata(hObject, handles);

%--------------------------------------------------------------------------

function MEGlist_Callback(hObject, eventdata, handles)
% hObject    handle to MEGlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MEGlist as text
%        str2double(get(hObject,'String')) returns contents of MEGlist as a double

% retrieve input
str = get(handles.MEGlist, 'String');

handles.D.modality = 'MEG';
handles = set_types(handles, 'MEG', str);
handles=update_list(handles);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

function HEOGlist_Callback(hObject, eventdata, handles)
% hObject    handle to HEOGlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HEOGlist as text
%        str2double(get(hObject,'String')) returns contents of HEOGlist as a double

% retrieve input
str = get(handles.HEOGlist, 'String');

handles = set_types(handles, 'HEOG', str);
handles=update_list(handles);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

function VEOGlist_Callback(hObject, eventdata, handles)
% hObject    handle to VEOGlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VEOGlist as text
%        str2double(get(hObject,'String')) returns contents of VEOGlist as a double

% retrieve input
str = get(handles.VEOGlist, 'String');

handles = set_types(handles, 'VEOG', str);
handles=update_list(handles);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

function Reflist_Callback(hObject, eventdata, handles)
% hObject    handle to Reflist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Reflist as text
%        str2double(get(hObject,'String')) returns contents of Reflist as a double

% retrieve input
str = get(handles.Reflist, 'String');

handles = set_types(handles, 'Ref', str);
handles=update_list(handles);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

function otherlist_Callback(hObject, eventdata, handles)
% hObject    handle to otherlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of otherlist as text
%        str2double(get(hObject,'String')) returns contents of otherlist as a double

% retrieve input
str = get(handles.otherlist, 'String');

handles = set_types(handles, 'other', str);
handles=update_list(handles);


% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

function newlist_Callback(hObject, eventdata, handles)
% hObject    handle to newlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of newlist as text
%        str2double(get(hObject,'String')) returns contents of newlist as a double

str = get(handles.newlist, 'String');

handles = set_types(handles, get(handles.newname, 'String'), str);
handles=update_list(handles);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'types')

    handles = set_types(handles, '', int2str([1:handles.Nnames]));

    handles = update_list(handles);

    % Update handles structure
    guidata(hObject, handles);

end

%--------------------------------------------------------------------------

% --- Executes on button press in default.
function default_Callback(hObject, eventdata, handles)
% hObject    handle to default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'types')
    handles=set_types_to_default(handles);
    handles=update_list(handles);
    guidata(hObject, handles);
end
%--------------------------------------------------------------------------

% --- Executes on button press in typeAssignmentNext.
function typeAssignmentNext_Callback(hObject, eventdata, handles)
% hObject    handle to typeAssignmentNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% check that all channels have been assigned to a type

% if ~alltypes(handles.types)
%     errordlg('Please assign a type to each channel');
% else
%     set_mode(handles, 'coordinates_assignment');
% end

set_mode(handles, 'coordinates_assignment');

% Update handles structure
guidata(hObject, handles);

% ------------------------ Assign 3D coordinates ---------------------------------

% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.textbox, 'String', '');
set(handles.textbox, 'Style', 'edit');
set(handles.textbox, 'HorizontalAlignment', 'left');
set(handles.textbox, 'Max', 2);

sensorFileTypes=get(handles.type3D, 'String');

switch sensorFileTypes{get(handles.type3D, 'Value')}

    case 'SPM EEG'

        % get sensor and fiducial locations
        %----------------------------------------------------------------------
        senspos = load(spm_select(1,'.mat','Select EEG sensors file'));
        name    = fieldnames(senspos);
        senspos = getfield(senspos,name{1});

        fiducials = load(spm_select(1,'.mat','Select EEG fiducials file'));
        name    = fieldnames(fiducials);
        fiducials = getfield(fiducials,name{1});

        handles=matchCoords(handles, fiducials, senspos);

    case 'SPM MEG'
        % get sensor and fiducial locations
        %----------------------------------------------------------------------
        senspos = load(spm_select(1,'.mat','Select MEG sensor locations file'));
        name    = fieldnames(senspos);
        senspos = getfield(senspos,name{1});

        sensorient = load(spm_select(1,'.mat','Select MEG sensor orientations file'));
        name    = fieldnames(sensorient);
        sensorient = getfield(sensorient,name{1});


        fiducials = load(spm_select(1,'.mat','Select MEG fiducials file'));
        name    = fieldnames(fiducials);
        fiducials = getfield(fiducials,name{1});


        handles=matchCoords(handles, fiducials, senspos, sensorient);
    case 'CTF folder'
        [senspos, sensorient, fiducials, labels] = spm_eeg_readCTFsensors;

        if isfield(handles, 'names')
            selection = match_str(lower(labels), lower(handles.names));
            senspos = senspos(selection, :);
            sensorient = sensorient(selection, :);

        else
            % This is for the case when no data file is loaded to enable
            % just creation of a template or coregistration of a location file
            handles.names=labels;
            handles.types=repmat({'MEG'}, 1, length(labels));
        end
        handles=matchCoords(handles, fiducials, senspos, sensorient);
    otherwise

        [filename ok]=spm_select(1,'.*','Select sensor positions file');

        positions_text=[];
        if ok
            fid=fopen(filename, 'r');
            while ~feof(fid)
                positions_text= [positions_text sprintf('%s\n', fgetl(fid))];
            end
            fclose(fid);
            set(handles.textbox, 'String', positions_text);
        end
end

guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in convert.
function convert_Callback(hObject, eventdata, handles)
% hObject    handle to convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

positions_text=get(handles.textbox, 'String');
handles.positions_text=positions_text;

filename=[tempname '.pol'];

fid=fopen(filename, 'w+');
for l=1:size(positions_text, 1)
    fprintf(fid, '%s\n', positions_text(l, :));
end
fclose(fid);

sensor_file_types=get(handles.type3D, 'String');
selected_type=sensor_file_types{get(handles.type3D, 'Value')};
selected_type = strrep(selected_type, 'EEGLAB ', '');
switch selected_type
    case 'FIL Polhemus'
        [fiducials, headshape] = spm_eeg_inv_ReadPolhemus(filename);
    otherwise
        try
            locstruct = eeglab_readlocs(filename , 'filetype', selected_type);
            locstruct = convertlocs(locstruct, 'auto');
        catch
            disp_conversion_message(handles, 'Conversion failed. See message in command window.');
            disp(getfield(lasterror, 'message'));
            delete(filename);
            return
        end

        fidNz=get(handles.fidNz, 'String');
        if isempty(str2num(fidNz))
            indNz = strmatch(lower(fidNz),lower({locstruct.labels}), 'exact');
            fidNz=[locstruct(indNz).X, locstruct(indNz).Y, locstruct(indNz).Z];
        else
            fidNz= str2num(fidNz);
            indNz=[];
        end

        fidNz=fidNz(:)';

        if ~(isnumeric(fidNz) && length(fidNz)==3)
            disp_conversion_message(handles, 'Failed to find the nasion coordinates.');
            return
        end

        fidL=get(handles.fidL, 'String');
        if isempty(str2num(fidL))
            indL = strmatch(lower(fidL),lower({locstruct.labels}), 'exact');
            fidL=[locstruct(indL).X, locstruct(indL).Y, locstruct(indL).Z];
        else
            fidL= str2num(fidL);
            indL=[];
        end

        fidL=fidL(:)';

        if ~(isnumeric(fidL) && length(fidL)==3)
            disp_conversion_message(handles, 'Failed to find the left fiducial coordinates.');
            return
        end


        fidR=get(handles.fidR, 'String');
        if isempty(str2num(fidR))
            indR = strmatch(lower(fidR),lower({locstruct.labels}), 'exact');
            fidR=[locstruct(indR).X, locstruct(indR).Y, locstruct(indR).Z];
        else
            fidR= str2num(fidR);
            indR=[];
        end

        fidR=fidR(:)';

        if ~(isnumeric(fidR) && length(fidR)==3)
            disp_conversion_message(handles, 'Failed to find the right fiducial coordinates.');
            return
        end

        if isfield(handles, 'names')

            if isfield(handles, 'types')
                channel_ind= strmatch('MEG', handles.types, 'exact');

                if isempty(channel_ind)
                    channel_ind= strmatch('EEG', handles.types, 'exact');
                end
            else
                channel_ind = 1:length(handles.names)
            end

            [junk selection] = match_str(lower(handles.names(channel_ind)), lower({locstruct.labels}));

            if length(selection)<length(channel_ind)
                disp_conversion_message(handles, 'The number of points converted is less than the number of EEG channels');
                return;
            end
        else
            % This is for the case when no data file is loaded to enable
            % just creation of a template or coregistration of a location file
            selection = setdiff(1:length(locstruct), [indNz indL indR]);
            handles.names={locstruct(selection).labels};
            handles.types=repmat({'EEG'}, 1, length(selection));
        end

        other_points=setdiff(1:length(locstruct), selection);

        fiducials = [fidNz; fidL; fidR];
        headshape=[[locstruct(selection).X]' [locstruct(selection).Y]' [locstruct(selection).Z]'];

        %other points are added to be used as headshape if necessary
        headshape=[headshape; [locstruct(other_points).X]' [locstruct(other_points).Y]' [locstruct(other_points).Z]'];
end
delete(filename);

handles = matchCoords(handles, fiducials, headshape);

guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in saveCoords.
function saveCoords_Callback(hObject, eventdata, handles)
% hObject    handle to saveCoords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'sensors')
    sensors=handles.sensors;
    [filename, pathname]=uiputfile('*.mat', 'Save sensors file as');
    if any(filename)
        save(fullfile(pathname, filename), 'sensors');
    end
end

if isfield(handles, 'fiducials')
    fiducials=handles.fiducials;
    [filename, pathname]=uiputfile('*.mat', 'Save fiducials file as');
    if any(filename)
        save(fullfile(pathname, filename), 'fiducials');
    end
end

if isfield(handles, 'sensorient') && ~isempty(handles.sensorient)
    sensorient=handles.sensorient;
    [filename, pathname]=uiputfile('*.mat', 'Save sensor orientations file as');
    if any(filename)
        save(fullfile(pathname, filename), 'sensorient');
    end
end

%--------------------------------------------------------------------------

% --- Executes on button press in coordAssignmentPrev.
function coordAssignmentPrev_Callback(hObject, eventdata, handles)
% hObject    handle to coordAssignmentPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_mode(handles, 'channel_type_assignment');

%--------------------------------------------------------------------------

% --- Executes on button press in coordAssignmentNext.
function coordAssignmentNext_Callback(hObject, eventdata, handles)
% hObject    handle to coordAssignmentNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_mode(handles, 'coregistration');

%=========================== Coregistration ===============================

% --- Executes on button press in coregister.
function coregister_Callback(hObject, eventdata, handles)
% hObject    handle to coregister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isfield(handles, 'D')
    D = handles.D;
else
    D=[];
end

if isfield(D, 'inv')
    invbackup=D.inv;
end


D.inv = cell(1);
val=1;

D.inv{val}.mesh.Msize = 1;
[D] = spm_eeg_inv_template(D, val);

allpoints = [handles.sensors; handles.fiducials; handles.headshape];

% The coregistration function doesn't like coordinates smaller than 1
% So the units are normalized to have values in the order of magnitude of 10.
norm_const = 0.1*mean(std(allpoints));

D.inv{val}.datareg.sensors = handles.sensors./norm_const;
D.inv{val}.datareg.fid_eeg = handles.fiducials./norm_const;
D.inv{val}.datareg.headshape = handles.headshape./norm_const;

if isfield(handles, 'sensorient')
    D.inv{val}.datareg.megorient=handles.sensorient;
else
    D.inv{val}.datareg.megorient=sparse(0,3);
end

[RT,sensors_reg,fid_reg,headshape_reg,orient_reg] = spm_eeg_inv_datareg(...
    D.inv{val}.datareg.sensors,...
    D.inv{val}.datareg.fid_eeg,...
    D.inv{val}.datareg.fid_mri,...
    D.inv{val}.datareg.headshape, ...
    D.inv{val}.datareg.scalpvert,...
    D.inv{val}.datareg.megorient,...
    D.inv{val}.mesh.template);


% D.channels.sensors = D.inv{val}.datareg.sensors; % original coordinates
% D.channels.sens_coreg = sensors_reg;
% D.channels.fid_coreg = fid_reg;

D.inv{val}.datareg.sens_coreg=sensors_reg;
D.inv{val}.datareg.fid_coreg=fid_reg;
D.inv{val}.datareg.hsp_coreg=headshape_reg;

handles.sensors_reg = sensors_reg;
handles.fid_reg = fid_reg;
handles.orient_reg=orient_reg;

if isfield(D, 'channels')
    spm_eeg_inv_checkdatareg(D);
end

% remove effects of inv_datareg
if exist('invbackup')==1
    D.inv=invbackup;
else
    D=rmfield(D, 'inv');
end

handles.D=D;
guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in saveRegCoord.
function saveRegCoord_Callback(hObject, eventdata, handles)
% hObject    handle to saveRegCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'sensors_reg')
    sensors=handles.sensors_reg;
    [filename, pathname]=uiputfile('*.mat', 'Save sensors file as');
    if any(filename)
        save(fullfile(pathname, filename), 'sensors');
    end
end

if isfield(handles, 'fid_reg')
    fiducials=handles.fid_reg;
    [filename, pathname]=uiputfile('*.mat', 'Save fiducials file as');
    if any(filename)
        save(fullfile(pathname, filename), 'fiducials');
    end
end

if isfield(handles, 'orient_reg') && ~isempty(handles.orient_reg)
    sensorient=handles.orient_reg;
    [filename, pathname]=uiputfile('*.mat', 'Save sensor orientations file as');
    if any(filename)
        save(fullfile(pathname, filename), 'sensorient');
    end
end

%--------------------------------------------------------------------------

% --- Executes on button press in coregistrationPrev.
function coregistrationPrev_Callback(hObject, eventdata, handles)
% hObject    handle to coregistrationPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_mode(handles, 'coordinates_assignment');

%--------------------------------------------------------------------------

% --- Executes on button press in coregistrationNext.
function coregistrationNext_Callback(hObject, eventdata, handles)
% hObject    handle to coregistrationNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_mode(handles, 'channel_template_preparation');

%----------------------- Assign 2D coordinates -----------------------------

% --- Executes on button press in load_template.
function load_template_Callback(hObject, eventdata, handles)
% hObject    handle to load_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

P = spm_select(1, '\.mat$', 'Select sensor template file', [], fullfile(spm('dir'), 'EEGtemplates'));
ctf=load(P); % must contain Cpos, Cnames

if isfield(handles, 'types')
    % At the moment the code does not support files that have both EEG and MEG
    % sensors thus only one type of sensors is handles and MEG has precedence.
    channel_ind= strmatch('MEG', handles.types, 'exact');

    if isempty(channel_ind)
        channel_ind= strmatch('EEG', handles.types, 'exact');
    end
    if isempty(channel_ind)
        return;
    end
    sel=match_str(ctf.Cnames, handles.names(channel_ind));
else
    sel=1:length(ctf.Cnames);
end

handles.Cpos = ctf.Cpos(:, sel);
handles.Cnames = ctf.Cnames(sel);
handles.Nchannels = length(sel);
handles.Rxy = ctf.Rxy;

set(handles.RxyBox, 'String', num2str(ctf.Rxy));

plot_sensors2D(handles)

setAxisCTF(handles);

guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in project3D.
function project3D_Callback(hObject, eventdata, handles)
% hObject    handle to project3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% At the moment the code does not support files that have both EEG and MEG
% sensors thus only one type of sensors is handles and MEG has precedence.

channel_ind= strmatch('MEG', handles.types, 'exact');

if isempty(channel_ind)
    channel_ind= strmatch('EEG', handles.types, 'exact');
end

if isempty(channel_ind)
    return;
end

cfg.elec.label=handles.names(channel_ind);

if get(handles.useCoregisteredCoord, 'Value')
    cfg.elec.pnt=handles.sensors_reg;
else
    cfg.elec.pnt=handles.sensors;
end

projection_types=get(handles.projectionType, 'String');

cfg.projection=projection_types{get(handles.projectionType, 'Value')};

lay = ft_prepare_layout(cfg);

[Cpos, Cnames, Nchannels, Rxy] =ftlay2ctf(lay);

handles.Cpos = Cpos;
handles.Cnames = Cnames;
handles.Nchannels = Nchannels;
handles.Rxy = Rxy;

plot_sensors2D(handles)

guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in clearTemplate.
function clearTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to clearTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.plotwin);
% Display

cla reset
hold on
handles.Rxy=str2num(get(handles.RxyBox, 'String'));
handles.Cnames={};
handles.Cpos=[];
handles.Nchannels=0;

handles.templateFrame=plot([0 0 1 1 0], [0 1 1 0 0], 'k-');

setAxisCTF(handles);

axis normal off

guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in undo_last.
function undo_last_Callback(hObject, eventdata, handles)
% hObject    handle to undo_last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'lastMoved')
    set(handles.lastMoved(end).point, 'XData', handles.lastMoved(end).coords(1));
    set(handles.lastMoved(end).point, 'YData', handles.lastMoved(end).coords(2));
    set(handles.lastMoved(end).label, 'Position', handles.lastMoved(end).coords);
    if length(handles.lastMoved)>1
        handles.lastMoved=handles.lastMoved(1:(end-1));
    else
        handles=rmfield(handles, 'lastMoved');
    end
    guidata(hObject, handles);
end

%--------------------------------------------------------------------------

% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
% hObject    handle to add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt={'Label:','X-position (0..1):', 'Y-position ([0...1]):'};
name='Define new label';
numlines=1;
defaultanswer={'', '0.5', '0.5'};

answer=inputdlg(prompt, name, numlines, defaultanswer);

if ~isempty(answer)
    handles.Cnames = [handles.Cnames answer(1)];
    handles.Cpos = [handles.Cpos [str2num(answer{2}); str2num(answer{3})]];

    handles.Nchannels = handles.Nchannels+1;
end

plot_sensors2D(handles)
setAxisCTF(handles);
guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'labelSelected') && ~isempty(handles.labelSelected)
    set(gcf, 'WindowButtonDownFcn', '');

    elec_label=get(handles.labelSelected, 'String');
    elec_ind=strmatch(elec_label, handles.Cnames, 'exact');

    delete([handles.labelSelected handles.pointSelected]);

    handles.Cnames(elec_ind)=[];
    handles.Cpos(:, elec_ind) = [];
    handles.Nchannels = handles.Nchannels-1;

    plot_sensors2D(handles)

    setAxisCTF(handles);

    guidata(hObject, handles);

end

guidata(hObject, handles);

%--------------------------------------------------------------------------

function RxyBox_Callback(hObject, eventdata, handles)
% hObject    handle to RxyBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RxyBox as text
%        str2double(get(hObject,'String')) returns contents of RxyBox as a double

try
    handles.Rxy=str2num(get(hObject, 'String'));
catch
    if isfield(handles, 'Rxy')
        set(hObject, 'String', num2str(handles.Rxy));
    else
        set(hObject, 'String', '1.5');
    end
end

setAxisCTF(handles);

guidata(hObject, handles);

%--------------------------------------------------------------------------

% --- Executes on button press in saveTemplate.
function saveTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to saveTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if all(isfield(handles, {'Cnames', 'Cpos'}))
    Cnames=handles.Cnames;
    Cpos=handles.Cpos;
    Nchannels=length(Cnames);
    if isfield(handles, 'Rxy')
        Rxy=handles.Rxy;
    else
        try
            Rxy=str2num(get(handles.RxyBox, 'String'));
        catch
            Rxy=1.5;
        end
    end
    [filename, pathname]=uiputfile('*.mat', 'Save template file as');
    if any(filename)
        save(fullfile(pathname, filename), 'Cnames', 'Cpos', 'Rxy', 'Nchannels');
    end
end

%--------------------------------------------------------------------------

% --- Executes on button press in channelTemplatePrev.
function channelTemplatePrev_Callback(hObject, eventdata, handles)
% hObject    handle to channelTemplatePrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set_mode(handles, 'coregistration');

%--------------------------------------------------------------------------

function labelClick(hObject, eventdata)

handles=guidata(hObject);

if isfield(handles, 'labelSelected') && ~isempty(handles.labelSelected)
    if handles.labelSelected==hObject
        set(handles.labelSelected, 'Color', 'r');
        set(handles.pointSelected,'MarkerFaceColor','r','MarkerSize',2,'MarkerEdgeColor','k');
        set(gcf, 'WindowButtonDownFcn', '');
    else
        handles.pointSelected=[];
        handles.labelSelected=[];
    end
else
    elec_label=get(hObject, 'String');
    set(gcf, 'WindowButtonDownFcn', @labelMove);

    coords = get(hObject, 'Position');

    handles.labelSelected=hObject;
    handles.pointSelected=findobj(gca, 'Type', 'line', 'XData', coords(1), 'YData', coords(2));

    set(handles.labelSelected, 'Color', 'g');
    set(handles.pointSelected,'MarkerFaceColor','g','MarkerSize',2,'MarkerEdgeColor','k');
end

guidata(hObject, handles);

%--------------------------------------------------------------------------

function labelMove(hObject, eventdata, handles)

handles=guidata(hObject);

coords=mean(get(handles.plotwin, 'CurrentPoint'));

set(handles.pointSelected, 'XData', coords(1));
set(handles.pointSelected, 'YData', coords(2));


set(handles.labelSelected, 'Position', coords);
set(handles.labelSelected, 'Color', 'r');
set(handles.pointSelected,'MarkerFaceColor','r','MarkerSize',2,'MarkerEdgeColor','k');
set(gcf, 'WindowButtonDownFcn', '');
set(gcf, 'WindowButtonMotionFcn', @cancelMove);


labelind=strmatch(get(handles.labelSelected, 'String'), handles.Cnames);

if isfield(handles, 'lastMoved')
    handles.lastMoved(end+1).point=handles.pointSelected;
    handles.lastMoved(end).label=handles.labelSelected;
    handles.lastMoved(end).coords=handles.Cpos(:, labelind);
else
    handles.lastMoved.point=handles.pointSelected;
    handles.lastMoved.label=handles.labelSelected;
    handles.lastMoved.coords=handles.Cpos(:, labelind);
end

handles.Cpos(:, labelind)=coords(1:2)';

guidata(hObject, handles);

%------------------------------------------------------------------------

function labelRelease(hObject, eventdata, handles)

handles=guidata(hObject);
set(handles.plotwin, 'ButtonDownFcn',  @labelMove);

%------------------------------------------------------------------------

function cancelMove(hObject, eventdata)
handles=guidata(hObject);

handles.pointSelected=[];
handles.labelSelected=[];
set(gcf, 'WindowButtonMotionFcn', '');
guidata(hObject, handles);

%============================ Utility functions ==========================

function set_mode(handles, mode)

switch mode
    case 'channel_type_assignment'
        set(findall(handles.typeAssignmentPanel, 'Enable', 'off'), 'Enable', 'on');
        set(findall(handles.coordAssignmentPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.coregistrationPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.channelTemplatePanel, 'Enable', 'on'), 'Enable', 'off');
    case 'coordinates_assignment'
        set(findall(handles.typeAssignmentPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.coordAssignmentPanel, 'Enable', 'off'), 'Enable', 'on');
        set(findall(handles.coregistrationPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.channelTemplatePanel, 'Enable', 'on'), 'Enable', 'off');
    case 'coregistration'
        set(findall(handles.typeAssignmentPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.coordAssignmentPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.coregistrationPanel, 'Enable', 'off'), 'Enable', 'on');
        set(findall(handles.channelTemplatePanel, 'Enable', 'on'), 'Enable', 'off');
        %plot_sensors(handles);
    case 'channel_template_preparation'
        set(findall(handles.typeAssignmentPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.coordAssignmentPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.coregistrationPanel, 'Enable', 'on'), 'Enable', 'off');
        set(findall(handles.channelTemplatePanel, 'Enable', 'off'), 'Enable', 'on');
        clearTemplate_Callback(handles.clearTemplate, [], handles);
end

%------------------------------------------------------------------------

function handles = set_types(handles, type, str);



if nargin>1
    ind = eval(['[' str ']']);
    for i = ind
        handles.types{i} = type;
    end
end

handles.types(find(cellfun('isempty', handles.types)))={''};

if length(handles.types)<length(handles.names)
    handles.types((end+1):length(handles.names)) = {''};
end

if strcmp(handles.D.modality, 'MEG')
    handles.D.channels.eeg = strmatch('MEG', handles.types, 'exact');
else
    handles.D.channels.eeg = strmatch('EEG', handles.types, 'exact');
end



%------------------------------------------------------------------------

function handles=update_list(handles);

Nspaces = handles.Nspaces;
% update list
for i = 1:length(handles.names)
    tmp1 = sprintf('%d', i); tmp1 = [tmp1 char(95*ones(1,Nspaces(1) - length(tmp1)))]; % running number
    tmp2 = [handles.names{i} char(95*ones(1,Nspaces(2) - length(handles.names{i})))]; % channel name
    tmp3 = [handles.types{i} char(95*ones(1,Nspaces(3) - length(handles.types{i})))]; % channel type

    textlb{i} = [tmp1 tmp2 tmp3];
end
handles.textlb = textlb;
set(handles.textbox, 'String', handles.textlb);

%------------------------------------------------------------------------

function handles=set_types_to_default(handles)

eeg_types = {'EEG', 'EEG1020', 'EEG1010', 'EEG1005', 'EEGCHWILLA', 'EEGBHAM', 'EEGREF'};
other_types = {'MEG', 'EMG', 'EOG'};

handles.types=repmat({''}, length(handles.names), 1);

for i=1:length(eeg_types)
    handles.types(match_str(handles.names, ft_channelselection(eeg_types(i), handles.names)))={'EEG'};
end

for i=1:length(other_types)
    handles.types(match_str(handles.names, ft_channelselection(other_types(i), handles.names)))=other_types(i);
end

%------------------------------------------------------------------------

function ok = alltypes(types)
ok = ~any(strcmp('', deblank(types)));

%------------------------------------------------------------------------

function handles=matchCoords(handles, fiducials, senspos, sensorient)

inputIsMEG=(exist('sensorient')==1);
dataFileLoaded=isfield(handles, 'types');

if inputIsMEG
    if dataFileLoaded
        channel_ind= strmatch('MEG', handles.types,'exact');

        if size(senspos, 1)~=length(channel_ind)
            disp_conversion_message(handles, 'The number of sensor position should exactly match the number of MEG sensors');
            return
        end
    else
        channel_ind=1:size(senspos, 1);
    end
else
    if dataFileLoaded
        channel_ind= strmatch('EEG', handles.types, 'exact');

        if size(senspos, 1)<length(channel_ind)
            disp_conversion_message(handles, 'The number of points converted is less than the number of EEG channels');
            return
        end
    else
        channel_ind=1:size(senspos, 1);
    end
end

disp_conversion_message(handles, 'ok');

handles.fiducials=fiducials;

if inputIsMEG
    handles.sensors=senspos;
    handles.fiducials=fiducials;
    handles.sensorient=sensorient;
    if ~isfield(handles, 'headshape')
        handles.headshape=sparse(0,3);
    end
else
    output_option=find([get(handles.sensors_and_headshape, 'Value'),...
        get(handles.sensors_only, 'Value'),...
        get(handles.headshape_only, 'Value')]);

    if output_option~=3 && ~isfield(handles, 'names')
        disp_conversion_message(handles, 'Channel labels are necessary to proceed.');
        return;
    end

    switch output_option
        case 1
            handles.sensors=senspos(1:length(channel_ind), :);
            handles.headshape=senspos;
        case 2
            handles.sensors=senspos(1:length(channel_ind), :);
            if ~isfield(handles, 'headshape')
                handles.headshape=sparse(0,3);
            end
        case 3
            handles.headshape=senspos;
            if ~isfield(handles, 'sensors')
                handles.sensors=[];
            end
    end
end

plot_sensors(handles);

%------------------------------------------------------------------------

function disp_conversion_message(handles, msg)

if nargin==1
    set(handles.convertmsg, 'String', '');
    set(handles.convertmsg, 'BackgroundColor', get(handles.coordAssignmentPanel,'BackgroundColor'));
elseif strcmp(msg, 'ok')
    set(handles.convertmsg, 'String', 'OK');
    set(handles.convertmsg, 'BackgroundColor', [0 1 0]);
else
    set(handles.convertmsg, 'BackgroundColor', [1 0 0]);
    set(handles.convertmsg, 'String', msg);
end

%------------------------------------------------------------------------

function plot_sensors(handles)


if isfield(handles, 'types')
    % At the moment the code does not support files that have both EEG and MEG
    % sensors thus only one type of sensors is handles and MEG has precedence.
    channel_ind= strmatch('MEG', handles.types, 'exact');

    if isempty(channel_ind)
        channel_ind= strmatch('EEG', handles.types, 'exact');
    end

    if isempty(channel_ind)
        return;
    end
end

axes(handles.plotwin);
% Display
%--------------------------------------------------------------------------
cla reset


if isfield(handles, 'headshape') && ~isempty(handles.headshape)
    h_hs = plot3(handles.headshape(:,1),handles.headshape(:,2),handles.headshape(:,3),'or');
    set(h_hs,'MarkerFaceColor','r','MarkerSize',4,'MarkerEdgeColor','k');

    axis([min(handles.headshape(:,1)), max(handles.headshape(:,1)), min(handles.headshape(:,2)),...
        max(handles.headshape(:,2)), min(handles.headshape(:,3)), max(handles.headshape(:,3))]);
end


if isfield(handles, 'sensors') && ~isempty(handles.sensors)
    text(handles.sensors(:, 1),handles.sensors(:,2),handles.sensors(:,3),strvcat(handles.names{channel_ind}),...
        'FontSize',8,...
        'Color','r',...
        'FontWeight','bold')

    axis([min(handles.sensors(:,1)), max(handles.sensors(:,1)), min(handles.sensors(:,2)),...
        max(handles.sensors(:,2)), min(handles.sensors(:,3)), max(handles.sensors(:,3))]);
end

hold on

text(handles.fiducials(:, 1), handles.fiducials(:,2),handles.fiducials(:,3), strvcat('Nz', 'L-Ear', 'R-Ear'),...
    'FontSize',8,...
    'Color','m',...
    'FontWeight','bold')


axis vis3d off
%view(-140,70)

drawnow

%------------------------------------------------------------------------

function plot_sensors2D(handles)

axes(handles.plotwin);
% Display

cla

handles.h_lbl=text(handles.Cpos(1, :),handles.Cpos(2, :),strvcat(handles.Cnames{:}),...
    'FontSize',8,...
    'Color','r',...
    'FontWeight','bold');

set(handles.h_lbl, 'ButtonDownFcn', @labelClick);

hold on

handles.h_el =[];
for i=1:size(handles.Cpos, 2)
    handles.h_el(i) = plot(handles.Cpos(1,i),handles.Cpos(2,i),'or');
end

set(handles.h_el,'MarkerFaceColor','r','MarkerSize',2,'MarkerEdgeColor','k');

handles.TemplateFrame=plot([0 0 1 1 0], [0 1 1 0 0], 'k-');

drawnow

guidata(gcf, handles);

%------------------------------------------------------------------------

function [Cpos, Cnames, Nchannels, Rxy] =ftlay2ctf(lay)

Cnames=lay.label(:)';

Cpos=[lay.pos(:,2), -lay.pos(:,1)];

Nchannels=length(Cnames);

Rxy=1.5;

Cpos=(Cpos-repmat(min(Cpos), Nchannels, 1));
Cpos=Cpos./repmat(max(Cpos), Nchannels, 1);
Cpos=Cpos*0.9+0.05;
Cpos=Cpos';

%------------------------------------------------------------------------

function setAxisCTF(handles)

axes(handles.plotwin);
if handles.Rxy>1
    axis([0.5+handles.Rxy.*[-0.5 0.5] 0 1]);
else
    axis([0 1 0.5+[-0.5 0.5]./handles.Rxy]);
end


%================== Initialization Functions ========================================

% --- Executes just before spm_eeg_assignchannels is made visible.
function spm_eeg_assignchannels_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spm_eeg_assignchannels (see VARARGIN)

% Choose default command line output for spm_eeg_assignchannels
handles.output = hObject;

try
    D = varargin{1};
    handles.D = D;
    handles.names = D.channels.name;
    handles.Nnames = length(handles.names);
    % populate list
    if isfield(D.channels, 'types')
        handles.types=D.channels.types;
    else
        handles=set_types_to_default(handles);
    end

    handles=update_list(handles);
catch
    set(handles.loadDataFile, 'Visible', 'on');
end

handles.Nspaces = [6 12 12]; % nr of characters to display 1) number, 2) channel name, 3) channel type

set_mode(handles, 'channel_type_assignment');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spm_eeg_assignchannels wait for user response (see UIRESUME)
%uiwait(handles.figure1);

%--------------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function plotwin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotwin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate plotwin

axis off

%--------------------------------------------------------------------------

function type3D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


set(hObject, 'String', {'SPM EEG', 'SPM MEG', 'CTF folder', 'FIL Polhemus', 'loc', 'sph', 'sfp', 'xyz', ...
    'asc', 'EEGLAB polhemusx', 'EEGLAB polhemusy', 'besa', 'EEGLAB chanedit'});

%--------------------------------------------------------------------------

function EEGlist_Callback(hObject, eventdata, handles)
% hObject    handle to EEGlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EEGlist as text
%        str2double(get(hObject,'String')) returns contents of EEGlist as a double

% retrieve input
str = get(handles.EEGlist, 'String');

if ~strcmp(handles.D.modality, 'MEG')
    handles.D.modality = 'EEG';
end

handles = set_types(handles, 'EEG', str);
update_list(handles);

% Update handles structure
guidata(hObject, handles);

%--------------------------------------------------------------------------


% --- Executes on button press in rotate.
function rotate_Callback(hObject, eventdata, handles)
% hObject    handle to rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rotate3d(handles.plotwin);



% --- Executes on button press in clearAssign.
function clearAssign_Callback(hObject, eventdata, handles)
% hObject    handle to clearAssign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try handles = rmfield(handles, 'sensors' ); end
try handles = rmfield(handles, 'headshape' ); end
try handles = rmfield(handles, 'fiducials' ); end

cla reset
axis vis3d off

disp_conversion_message(handles);

guidata(gcf, handles);


% --- Executes on button press in saveTypes.
function saveTypes_Callback(hObject, eventdata, handles)
% hObject    handle to saveTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'D')
    D=handles.D;
    [filename, pathname]=uiputfile('*.mat', 'Save SPM header file as');
    if any(filename)
        save(fullfile(pathname, filename), 'D');
    end
end

