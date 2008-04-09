function spm_eeg_prep_ui(callback)
% User interface for spm_eeg_prep function performing several tasks
% for preparation of converted MEEG data for further analysis
% FORMAT spm_eeg_prep_ui()
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_prep_ui.m 1326 2008-04-09 12:15:35Z vladimir $

if nargin == 0

    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','MEEG preparation',0);

    delete(findobj(get(Finter,'Children'),'Tag','EEGprepUI'))

    %-Draw top level menu
    % ====== File ===================================
    FileMenu = uimenu(Finter,'Label','File',...
        'Tag','EEGprepUI',...
        'HandleVisibility','on');

    FileOpenMenu = uimenu(FileMenu, ...
        'Label','Open',...
        'Separator','off',...
        'Tag','EEGprepUI',...
        'HandleVisibility', 'on',...
        'Callback', 'spm_eeg_prep_ui(''FileOpenCB'')');

    FileSaveMenu = uimenu(FileMenu, ...
        'Label','Save',...
        'Separator','off',...
        'Tag','EEGprepUI',...
        'HandleVisibility', 'on',...
        'Callback', 'spm_eeg_prep_ui(''FileSaveCB'')');

    % ====== Channel types ===============================
    
    ChanTypeMenu = uimenu(Finter,'Label','Channel types',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on');

    chanTypes = {'EEG', 'VEOG', 'HEOG', 'LFP', 'Other'};

    for i = 1:length(chanTypes)
        CTypesMenu(i) = uimenu(ChanTypeMenu, 'Label', chanTypes{i},...
            'Tag','EEGprepUI',...
            'Enable', 'on', ...
            'HandleVisibility','on',...
            'Callback', 'spm_eeg_prep_ui(''ChanTypeCB'')');
    end

    CTypesReviewMenu = uimenu(ChanTypeMenu, 'Label', 'Review',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on',...
        'Callback', 'spm_eeg_prep_ui(''ChanTypeCB'')');

    % ====== Sensors ===================================
    
    Coor3DMenu = uimenu(Finter,'Label','Sensors',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on');

    LoadEEGSensMenu = uimenu(Coor3DMenu, 'Label', 'Load (EEG)',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on');

    LoadEEGSensMatMenu = uimenu(LoadEEGSensMenu, 'Label', 'From *.mat files',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''LoadEEGSensCB'')');

    LoadEEGSensPolhemusMenu = uimenu(LoadEEGSensMenu, 'Label', 'From FIL polhemus file',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''LoadEEGSensCB'')');

    HeadshapeMenu = uimenu(Coor3DMenu, 'Label', 'Headshape',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on');

    LoadHeadshapeMatMenu = uimenu(HeadshapeMenu, 'Label', 'Load from *.mat file',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on');

    LoadHeadshapePolhemusMenu = uimenu(HeadshapeMenu, 'Label', 'Load from FIL polhemus file',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''HeadshapeCB'')');

    CoregisterMenu = uimenu(Coor3DMenu, 'Label', 'Coregister',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on', ...
        'Callback', 'spm_eeg_prep_ui(''CoregisterCB'')');

    % ====== 2D projection ===================================
    
    Coor2DMenu = uimenu(Finter, 'Label','2D projection',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on');

    EditMEGMenu = uimenu(Coor2DMenu, 'Label', 'Edit existing MEG',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''EditExistingCoor2DCB'')');

    EditEEGMenu = uimenu(Coor2DMenu, 'Label', 'Edit existing EEG',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''EditExistingCoor2DCB'')');

    LoadTemplateMenu = uimenu(Coor2DMenu, 'Label', 'Load template',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on', ...
        'Callback', 'spm_eeg_prep_ui(''LoadTemplateCB'')');

    SaveTemplateMenu = uimenu(Coor2DMenu, 'Label', 'Save template',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''SaveTemplateCB'')');

    Project3DEEGMenu = uimenu(Coor2DMenu, 'Label', 'Project 3D (EEG)',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on', ...
        'Callback', 'spm_eeg_prep_ui(''Project3DCB'')');

    Project3DMEGMenu = uimenu(Coor2DMenu, 'Label', 'Project 3D (MEG)',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''Project3DCB'')');

    AddCoor2DMenu = uimenu(Coor2DMenu, 'Label', 'Add sensor',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on', ...
        'Callback', 'spm_eeg_prep_ui(''AddCoor2DCB'')');

    DeleteCoor2DMenu = uimenu(Coor2DMenu, 'Label', 'Delete sensor',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''DeleteCoor2DCB'')');

    UndoMoveCoor2DMenu = uimenu(Coor2DMenu, 'Label', 'Undo move',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''UndoMoveCoor2DCB'')');

    ApplyCoor2DMenu = uimenu(Coor2DMenu, 'Label', 'Apply',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on', ...
        'Callback', 'spm_eeg_prep_ui(''ApplyCoor2DCB'')');

    Clear2DMenu = uimenu(Coor2DMenu, 'Label', 'Clear',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''Clear2DCB'')');

else
    eval(callback);
end

%-----------------------------------------------------------------------

function FileOpenCB()

try
    D = spm_eeg_load(spm_select(1, 'mat', 'Select M/EEG mat file'));
    setD(D);
end
update_menu;

%-----------------------------------------------------------------------

function FileSaveCB()
D = getD;
if ~isempty(D)
    D.save;
end

update_menu;

%-----------------------------------------------------------------------

function ChanTypeCB

type = get(gcbo, 'Label');

D = getD;

if ~isempty(D)
    chanlist ={};
    for i = 1:D.nchannels
        if strncmp(D.chantype(i), 'MEG', 3)
            chanlist{i} = [num2str(i) '    Label:    ' D.chanlabels(i) '    Type:    ' D.chantype(i) , ' (nonmodifiable)'];
        else
            chanlist{i} = [num2str(i) '    Label:    ' D.chanlabels(i) '    Type:    ' D.chantype(i)];
        end
    end

    if strcmpi(type, 'review')
        listdlg('ListString', chanlist, 'SelectionMode', 'single', 'Name', 'Review channels', 'ListSize', [400 300]);
        return
    else

        [selection ok]= listdlg('ListString', chanlist, 'SelectionMode', 'multiple',...
            'InitialValue', strmatch(type, D.chantype) ,'Name', ['Set type to ' type], 'ListSize', [400 300]);

        % Changing the type of MEG channels in GUI is not allowed.
        selection(strmatch('MEG', chantype(D, selection))) = [];
        if ok && ~isempty(selection)
            S.task = 'settype';
            S.D = D;
            S.ind = selection;
            S.type = type;
            D = spm_eeg_prep(S);
            setD(D);
        end
    end
end

update_menu;

%-----------------------------------------------------------------------

function LoadEEGSensCB

S = [];
S.D = getD;
S.task = 'loadeegsens';

switch get(gcbo, 'Label')
    case 'From *.mat files'
        S.sensfile{1} = spm_select(1,'.mat$','Select EEG fiducials file');
        S.sensfile{2} = spm_select(1,'.mat$','Select EEG sensors file');
        S.source = 'mat';
    case 'From FIL polhemus file'
        S.sensfile = spm_select(1, '\.pol$', 'Select FIL polhemus file');
        S.source = 'filpolhemus'
end

D = spm_eeg_prep(S);

setD(D);

update_menu;

%-----------------------------------------------------------------------

function HeadshapeCB

S = [];
S.D = getD;
S.task = 'headshape';

switch get(gcbo, 'Label')
    case 'Load from *.mat file'
        S.headshapefile = spm_select(1,'.mat$','Select matlab head');
        S.source = 'mat';
    case 'Load from FIL polhemus file'
        S.headshapefile = spm_select(1, '\.pol$', 'Select FIL polhemus file');
        S.source = 'filpolhemus';
end

D = spm_eeg_prep(S);

setD(D);

update_menu;

%-----------------------------------------------------------------------

function CoregisterCB

S = [];
S.D = getD;
S.task = 'coregister';

D = spm_eeg_prep(S);

setD(D);

update_menu;

%-----------------------------------------------------------------------

function EditExistingCoor2DCB

D = getD;

switch get(gcbo, 'Label')
    case 'Edit existing MEG'
        xy = D.coor2D('MEG');
        label = D.chanlabels(strmatch('MEG', D.chantype, 'exact'));
    case 'Edit existing EEG'
        xy = D.coor2D('EEG');
        label = D.chanlabels(strmatch('EEG', D.chantype, 'exact'));
end

plot_sensors2D(xy, label);

update_menu;

%-----------------------------------------------------------------------

function LoadTemplateCB

template = [];

template = load(spm_select(1, '\.mat$', 'Select sensor template file', ...
    [], fullfile(spm('dir'), 'EEGtemplates')));


if isfield(template, 'Cnames') && isfield(template, 'Cpos')
    plot_sensors2D(template.Cpos, template.Cnames);
end

update_menu;

%-----------------------------------------------------------------------

function SaveTemplateCB

handles=getHandles;

Cnames = handles.label;
Cpos = handles.xy;
Rxy = 1.5;
Nchannels = length(Cnames);

[filename, pathname] = uiputfile('*.mat', 'Save channel template as');

save(fullfile(pathname, filename), 'Cnames', 'Cpos', 'Rxy', 'Nchannels');

%-----------------------------------------------------------------------

function Project3DCB

D = getD;

% Ugly, but this is a very specific case that does not justify
% making a separate method at this stage
sD = struct(D);
sensors = sD.sensors;

switch get(gcbo, 'Label')
    case 'Project 3D (EEG)'
        ind = strmatch('EEG', sensors.type, 'exact');
    case 'Project 3D (MEG)'
        ind = strmatch('MEG', sensors.type, 'exact');
end

cfg = [];
cfg.elec = [];
cfg.elec.label = sensors.label(ind);
cfg.elec.pnt = sensors.pnt(ind, :);

lay = ft_prepare_layout(cfg);

[sel1, sel2] = spm_match_str(cfg.elec.label, lay.label);

label =lay.label(sel2)';
xy = [lay.pos(sel2, 2), -lay.pos(sel2, 1)];

nchan = size(xy, 1);

xy =(xy-repmat(min(xy), nchan, 1));
xy = xy./repmat(max(xy), nchan, 1);
xy = xy*0.9+0.05;
xy = xy';

plot_sensors2D(xy, label);

update_menu;

%-----------------------------------------------------------------------

function AddCoor2DCB

newlabel = spm_input('Label?', '+1', 's');

if isempty(newlabel)
    return;
end

coord = spm_input('Coordinates [x y]', '+1', 'r', '0.5 0.5', 2);

handles = getHandles;

if ~isfield(handles, 'xy')
    handles.xy = [];
end

if ~isfield(handles, 'xy')
    handles.xy = [];
end

if ~isfield(handles, 'label')
    handles.label = {};
end

plot_sensors2D([handles.xy coord(:)], ...
    [handles.label newlabel]);

update_menu;

%-----------------------------------------------------------------------

function ApplyCoor2DCB

handles = getHandles;
D = getD;

S = [];
S.task = 'setcoor2d';
S.D = D;
S.xy = handles.xy;
S.label = handles.label;

D = spm_eeg_prep(S);

setD(D);

update_menu;

%-----------------------------------------------------------------------

function update_menu

Finter = spm_figure('GetWin','Interactive');
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'File'), 'Enable', 'on');

IsEEG = 'off';
IsMEG = 'off';
IsSensors = 'off';
IsSensorsEEG = 'off';
IsSensorsMEG = 'off';
if isa(get(Finter, 'UserData'), 'meeg')
    Dloaded = 'on';

    D = getD;

    if ~isempty(strmatch('EEG', D.chantype, 'exact'))
        IsEEG = 'on';
    end

    if ~isempty(strmatch('MEG', D.chantype, 'exact'));
        IsMEG = 'on';
    end

    if D.sensorsdefined
        IsSensors = 'on';
    end

    if D.sensorsdefined('EEG')
        IsSensorsEEG = 'on';
    end

    if  D.sensorsdefined('MEG')
        IsSensorsMEG = 'on';
    end
else
    Dloaded = 'off';
end

handles = getHandles;

IsTemplate = 'off';
IsSelected = 'off';
IsMoved = 'off';
if ~isempty(handles)
    if isfield(handles, 'xy') && size(handles.xy, 1)>0
        IsTemplate = 'on';
    end

    if isfield(handles, 'labelSelected') && ~isempty(handles.labelSelected)
        IsSelected = 'on';
    end

    if isfield(handles, 'lastMoved')
        isMoved = 'on';
    end
end

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Channel types'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Sensors'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', '2D projection'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load (EEG)'), 'Enable', IsEEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Headshape'), 'Enable', IsSensors);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Coregister'), 'Enable', IsSensors);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing EEG'), 'Enable', IsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing MEG'), 'Enable', IsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (EEG)'), 'Enable', IsSensorsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (MEG)'), 'Enable', IsSensorsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Delete sensor'), 'Enable', IsSelected);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Undo move'), 'Enable', IsMoved);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Apply'), 'Enable', IsTemplate);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Clear'), 'Enable', IsTemplate);

delete(setdiff(findobj(Finter), [Finter; findobj(Finter,'Tag','EEGprepUI')]));

figure(Finter);


%-----------------------------------------------------------------------

function D = getD()
Finter = spm_figure('GetWin','Interactive');
D = get(Finter, 'UserData');
if ~isa(D, 'meeg')
    D = [];
end

%-----------------------------------------------------------------------

function setD(D)
Finter = spm_figure('GetWin','Interactive');
set(Finter, 'UserData', D);

%-----------------------------------------------------------------------

function handles = getHandles()
Fgraph = spm_figure('GetWin','Graphics');
handles = get(Fgraph, 'UserData');

%-----------------------------------------------------------------------

function setHandles(handles)
Fgraph = spm_figure('GetWin','Graphics');
set(Fgraph, 'UserData', handles);

%-----------------------------------------------------------------------

function plot_sensors2D(xy, label)

Fgraph = spm_figure('GetWin','Graphics');

figure(Fgraph);
clf

handles = [];

if ~isempty(xy)
    if size(xy, 1) ~= 2
        xy = xy';
    end
    
    
    xy(xy < 0.05) = 0.05;
    xy(xy > 0.95) = 0.95;


    handles.h_lbl=text(xy(1,:), xy(2, :),strvcat(label),...
        'FontSize', 9,...
        'Color','r',...
        'FontWeight','bold');

    set(handles.h_lbl, 'ButtonDownFcn', 'spm_eeg_prep_ui(''LabelClickCB'')');

    hold on

    handles.h_el =[];
    for i=1:size(xy, 2)
        handles.h_el(i) = plot(xy(1,i), xy(2,i), 'or');
    end

    set(handles.h_el,'MarkerFaceColor','r','MarkerSize', 2,'MarkerEdgeColor','k');

end

handles.TemplateFrame = ...
    plot([0.05 0.05 0.95 0.95 0.05], [0.05 0.95 0.95 0.05 0.05], 'k-');

axis off;

handles.xy = xy;
handles.label = label(:)';

setHandles(handles);

update_menu;

%-----------------------------------------------------------------------

function DeleteCoor2DCB

handles = getHandles;
graph = spm_figure('GetWin','Graphics');

if isfield(handles, 'labelSelected') && ~isempty(handles.labelSelected)
    set(graph, 'WindowButtonDownFcn', '');

    label=get(handles.labelSelected, 'String');
    ind=strmatch(label, handles.label, 'exact');

    delete([handles.labelSelected handles.pointSelected]);

    handles.xy(:, ind)=[];
    handles.label(ind) = [];

    plot_sensors2D(handles.xy, handles.label)
end

%-----------------------------------------------------------------------

function UndoMoveCoor2DCB

handles = getHandles;

if isfield(handles, 'lastMoved')
    label = get(handles.lastMoved(end).label, 'String');
    ind = strmatch(label, handles.label, 'exact');
    handles.xy(:, ind) = handles.lastMoved(end).coords(:);

    set(handles.lastMoved(end).point, 'XData', handles.lastMoved(end).coords(1));
    set(handles.lastMoved(end).point, 'YData', handles.lastMoved(end).coords(2));
    set(handles.lastMoved(end).label, 'Position', handles.lastMoved(end).coords);

    if length(handles.lastMoved)>1
        handles.lastMoved = handles.lastMoved(1:(end-1));
    else
        handles = rmfield(handles, 'lastMoved');
    end

    setHandles(handles);
    update_menu;
end

%-----------------------------------------------------------------------

function LabelClickCB

handles=getHandles;
Fgraph = spm_figure('GetWin','Graphics');

if isfield(handles, 'labelSelected') && ~isempty(handles.labelSelected)
    if handles.labelSelected == gcbo
        set(handles.labelSelected, 'Color', 'r');
        set(handles.pointSelected,'MarkerFaceColor', 'r');
        set(Fgraph, 'WindowButtonDownFcn', '');
    else
        handles.pointSelected=[];
        handles.labelSelected=[];
    end
else
    set(Fgraph, 'WindowButtonDownFcn', 'spm_eeg_prep_ui(''LabelMoveCB'')');

    coords = get(gcbo, 'Position');

    handles.labelSelected=gcbo;
    handles.pointSelected=findobj(gca, 'Type', 'line',...
        'XData', coords(1), 'YData', coords(2));

    set(handles.labelSelected, 'Color', 'g');
    set(handles.pointSelected,'MarkerFaceColor', 'g');
end

setHandles(handles);

update_menu;

%-----------------------------------------------------------------------

function LabelMoveCB

handles = getHandles;
Fgraph = spm_figure('GetWin','Graphics');

coords=mean(get(gca, 'CurrentPoint'));

coords(coords < 0.05) = 0.05;
coords(coords > 0.95) = 0.95;

set(handles.pointSelected, 'XData', coords(1));
set(handles.pointSelected, 'YData', coords(2));


set(handles.labelSelected, 'Position', coords);
set(handles.labelSelected, 'Color', 'r');
set(handles.pointSelected,'MarkerFaceColor','r','MarkerSize',2,'MarkerEdgeColor','k');
set(Fgraph, 'WindowButtonDownFcn', '');
set(Fgraph, 'WindowButtonMotionFcn', 'spm_eeg_prep_ui(''CancelMoveCB'')');


labelind=strmatch(get(handles.labelSelected, 'String'), handles.label);

if isfield(handles, 'lastMoved')
    handles.lastMoved(end+1).point = handles.pointSelected;
    handles.lastMoved(end).label = handles.labelSelected;
    handles.lastMoved(end).coords = handles.xy(:, labelind);
else
    handles.lastMoved.point = handles.pointSelected;
    handles.lastMoved.label = handles.labelSelected;
    handles.lastMoved.coords = handles.xy(:, labelind);
end

handles.xy(:, labelind) = coords(1:2)';

setHandles(handles);

update_menu;
%------------------------------------------------------------------------

function CancelMoveCB

Fgraph = spm_figure('GetWin','Graphics');
handles = getHandles;

handles.pointSelected=[];
handles.labelSelected=[];
set(Fgraph, 'WindowButtonMotionFcn', '');

setHandles(handles);

update_menu;
%------------------------------------------------------------------------

function Clear2DCB

plot_sensors2D([], {});

update_menu;

