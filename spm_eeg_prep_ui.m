function spm_eeg_prep_ui(callback)
% User interface for spm_eeg_prep function performing several tasks
% for preparation of converted MEEG data for further analysis
% FORMAT spm_eeg_prep_ui()
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_prep_ui.m 2327 2008-10-10 15:24:53Z jean $

spm('pointer','watch')

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

    CTypesRef2MEGMenu = uimenu(ChanTypeMenu, 'Label', 'MEGREF=>MEG',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on',...
        'Separator', 'on',...
        'Callback', 'spm_eeg_prep_ui(''MEGChanTypeCB'')');

    CTypesPlanar2MEGMenu = uimenu(ChanTypeMenu, 'Label', 'Planar=>MEG',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''MEGChanTypeCB'')');

    CTypesMagnetometer2MEGMenu = uimenu(ChanTypeMenu, 'Label', 'Magnetometer=>MEG',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''MEGChanTypeCB'')');

    CTypesDefaultMenu = uimenu(ChanTypeMenu, 'Label', 'Default',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on',...
        'Callback', 'spm_eeg_prep_ui(''ChanTypeDefaultCB'')');

    CTypesReviewMenu = uimenu(ChanTypeMenu, 'Label', 'Review',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''ChanTypeCB'')');
    % ====== Sensors ===================================

    Coor3DMenu = uimenu(Finter,'Label','Sensors',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on');

    LoadEEGSensMenu = uimenu(Coor3DMenu, 'Label', 'Load EEG sensors',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on');

    LoadEEGSensTemplateMenu = uimenu(LoadEEGSensMenu, 'Label', 'Assign default',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''LoadEEGSensTemplateCB'')');

    LoadEEGSensMatMenu = uimenu(LoadEEGSensMenu, 'Label', 'From *.mat file',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''LoadEEGSensCB'')');

    LoadEEGSensOtherMenu = uimenu(LoadEEGSensMenu, 'Label', 'Convert locations file',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''LoadEEGSensCB'')');

    HeadshapeMenu = uimenu(Coor3DMenu, 'Label', 'Load MEG Fiducials/Headshape',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''HeadshapeCB'')');

    CoregisterEEGMenu = uimenu(Coor3DMenu, 'Label', 'Coregister (EEG)',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Separator', 'on', ...
        'Callback', 'spm_eeg_prep_ui(''CoregisterCB'')');

    CoregisterMEGMenu = uimenu(Coor3DMenu, 'Label', 'Coregister (MEG)',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
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

    % ====== History ===================================

    HistoryMenu = uimenu(Finter, 'Label','History',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on');

    HistoryReviewMenu = uimenu(HistoryMenu, 'Label','Review history',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''HistoryCB'')');

    History2ScriptMenu = uimenu(HistoryMenu, 'Label','Save as script',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''HistoryCB'')');
else
    eval(callback);
end


spm('pointer','arrow')

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

        chanlist{i} = [chanlist{i}{:}];
    end

    if strcmpi(type, 'review')
        listdlg('ListString', chanlist, 'SelectionMode', 'single', 'Name', 'Review channels', 'ListSize', [400 300]);
        return
    else

        [selection ok]= listdlg('ListString', chanlist, 'SelectionMode', 'multiple',...
            'InitialValue', strmatch(type, D.chantype) ,'Name', ['Set type to ' type], 'ListSize', [400 300]);

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

function MEGChanTypeCB

S = [];
S.D = getD;
S.task = 'settype';

switch get(gcbo, 'Label')
    case 'MEGREF=>MEG'
        S.ind = strmatch('MEGREF', S.D.chantype);
        S.type = 'MEG';
        D = spm_eeg_prep(S);
    case 'Planar=>MEG'
        S.type = 'MEG';
        S.ind = strmatch('planar', S.D.origchantypes, 'exact');
        S.D = spm_eeg_prep(S);
        S.type = 'Other';
        S.ind = strmatch('magnetometer' , S.D.origchantypes, 'exact');
        D = spm_eeg_prep(S);
    case 'Magnetometer=>MEG'
        S.type = 'MEG';
        S.ind = strmatch('magnetometer', S.D.origchantypes, 'exact');
        S.D = spm_eeg_prep(S);
        S.type = 'Other';
        S.ind = strmatch('planar', S.D.origchantypes, 'exact');
        D = spm_eeg_prep(S);
end

setD(D);
update_menu;

%-----------------------------------------------------------------------

function ChanTypeDefaultCB

S = [];
S.D = getD;
S.task = 'settype';
S.type = [];
S.ind = 1:S.D.nchannels;
D = spm_eeg_prep(S);

setD(D);
update_menu;

%-----------------------------------------------------------------------

function LoadEEGSensTemplateCB

S = [];
S.D = getD;
S.task = 'defaulteegsens';

D = spm_eeg_prep(S);
setD(D);
update_menu;

%-----------------------------------------------------------------------

function LoadEEGSensCB

S = [];
S.D = getD;
S.task = 'loadeegsens';

switch get(gcbo, 'Label')
    case 'From *.mat file'
        S.sensfile = spm_select(1,'.mat$','Select EEG sensors file');
        S.source = 'mat';
    case 'Convert locations file'
        S.sensfile = spm_select(1, '\.*', 'Select locations file');
        S.source = 'locfile';
end

D = spm_eeg_prep(S);

% ====== This is for the future ==================================
% sens = D.sensors('EEG');
% label = D.chanlabels(strmatch('EEG',D.chantype));
%
% [sel1, sel2] = spm_match_str(label, sens.label);
%
% montage = [];
% montage.labelorg = sens.label;
% montage.labelnew = label;
% montage.tra = sparse(zeros(numel(label), numel(sens.label)));
% montage.tra(sub2ind(size(montage.tra), sel1, sel2)) = 1;
%
% montage = spm_eeg_montage_ui(montage);
%
% S = [];
% S.D = D;
% S.task = 'sens2chan';
% S.montage = montage;
%
% D = spm_eeg_prep(S);

% ============= Assign the fiducials using a mat file or the sensors file

S.task = 'headshape';
S.D = D;
S.regfid = {};
if strcmp(S.source, 'mat')
    S.headshapefile = spm_select(1,'.mat$','Select EEG fiducials file');
    S.fidlabel = spm_input('Fiducial labels:', '+1', 's', 'nas lpa rpa');
else
    S.headshapefile = S.sensfile;
    S.source = 'convert';
end

D = spm_eeg_prep(S);

setD(D);

update_menu;

%-----------------------------------------------------------------------

function HeadshapeCB

S = [];
S.D = getD;
S.task = 'headshape';

S.headshapefile = spm_select(1, '\.*', 'Select fiducials/headshape file');
S.source = 'convert';

shape = fileio_read_headshape(S.headshapefile);
lblshape = shape.fid.label;

fid = fiducials(S.D);
lblfid = fid.fid.label;

if numel(intersect(upper(lblshape), upper(lblfid))) < 3
    if numel(lblshape)<3 || numel(lblfid)<3
        warndlg('3 fiducials are required to load headshape');
        return;
    else
        S.regfid = {};
        for i = 1:length(lblfid)
            [selection ok]= listdlg('ListString',lblshape, 'SelectionMode', 'single',...
                'InitialValue', strmatch(upper(lblfid{i}), upper(lblshape)), ...
                'Name', ['Select matching fiducial for ' lblfid{i}], 'ListSize', [400 300]);

            if ~ok
                continue
            end

            S.regfid = [S.regfid; [lblfid(i) lblshape(selection)]];
        end

        if size(S.regfid, 1) < 3
            warndlg('3 fiducials are required to load headshape');
            return;
        end
    end
else
    [sel1, sel2] = spm_match_str(upper(lblfid), upper(lblshape));
    lblfid = lblfid(sel1);
    lblshape = lblshape(sel2);
    S.regfid = [lblfid(:) lblshape(:)];
end

D = spm_eeg_prep(S);

setD(D);

update_menu;

%-----------------------------------------------------------------------

function CoregisterCB

S = [];
S.D = getD;
S.task = 'coregister';

switch get(gcbo, 'Label')
    case 'Coregister (EEG)'
        S.modality = 'EEG';
    case 'Coregister (MEG)'
        S.modality = 'MEG';
end


D = spm_eeg_prep(S);

% Bring the menu back
spm_eeg_prep_ui;

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

switch get(gcbo, 'Label')
    case 'Project 3D (EEG)'
        modality = 'EEG';
    case 'Project 3D (MEG)'
        modality = 'MEG';
end

[xy, label] = spm_eeg_project3D(D.sensors(modality), modality);

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
IsNeuromag = 'off';
HasSensors = 'off';
HasSensorsEEG = 'off';
HasSensorsMEG = 'off';
HasChannelsMEGREF = 'off';
HasFiducials = 'off';
HasDefaultLocs = 'off';
HasHistory = 'off';
if isa(get(Finter, 'UserData'), 'meeg')
    Dloaded = 'on';

    D = getD;

    if ~isempty(strmatch('EEG', D.chantype, 'exact'))
        IsEEG = 'on';
    end

    if ~isempty(strmatch('MEG', D.chantype, 'exact'));
        IsMEG = 'on';
    end

    if forwinv_senstype(D.chanlabels, 'neuromag') &&...
            isfield(D, 'origchantypes')
        IsNeuromag = 'on';
    end

    if ~isempty(strmatch('MEGREF', D.chantype, 'exact'));
        HasChannelsMEGREF = 'on';
    end

    if ~isempty(D.sensors('EEG')) || ~isempty(D.sensors('MEG'))
        HasSensors = 'on';
    end

    if ~isempty(D.sensors('EEG'))
        HasSensorsEEG = 'on';
    end

    if  ~isempty(D.sensors('MEG'))
        HasSensorsMEG = 'on';
    end

    if  ~isempty(D.fiducials)
        HasFiducials = 'on';
    end

    template_sfp = dir(fullfile(spm('dir'), 'EEGtemplates', '*.sfp'));
    template_sfp = {template_sfp.name};
    ind = strmatch([forwinv_senstype(D.chanlabels) '.sfp'], template_sfp, 'exact');

    if ~isempty(ind)
        HasDefaultLocs = 'on';
    end

    if  ~isempty(D.history)
        HasHistory = 'on';
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

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'MEGREF=>MEG'), 'Enable', HasChannelsMEGREF);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Planar=>MEG'), 'Enable', IsNeuromag);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Magnetometer=>MEG'), 'Enable', IsNeuromag);


set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Assign default'), 'Enable', HasDefaultLocs);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load EEG sensors'), 'Enable', IsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load MEG Fiducials/Headshape'), 'Enable', HasSensorsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Headshape'), 'Enable', HasSensors);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Coregister (EEG)'), 'Enable', HasSensorsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Coregister (MEG)'), 'Enable', HasSensorsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing EEG'), 'Enable', IsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing MEG'), 'Enable', IsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (EEG)'), 'Enable', HasSensorsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (MEG)'), 'Enable', HasSensorsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Delete sensor'), 'Enable', IsSelected);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Undo move'), 'Enable', IsMoved);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Apply'), 'Enable', IsTemplate);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Clear'), 'Enable', IsTemplate);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'History'), 'Enable', HasHistory);

delete(setdiff(findobj(Finter), [Finter; findobj(Finter,'Tag','EEGprepUI')]));

if strcmp(Dloaded, 'on') && isfield(D,'PSD') && D.PSD == 1
    try
        hc = get(Finter,'children');
        hc = findobj(hc,'flat','type','uimenu');
        hc = findobj(hc,'flat','label','File');
        delete(hc)
    end
    uicontrol(Finter,...
        'style','pushbutton','string','OK',...
        'callback','spm_eeg_review_callbacks(''get'',''prep'')',...
        'tooltipstring','Update data informations in ''SPM Graphics'' window',...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'Tag','EEGprepUI');
end


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

% ====== History ===============================

function HistoryCB

D = getD;

h = D.history;

histlist ={};
for i = 1:numel(h)
    switch h(i).fun
        case 'spm_eeg_convert'
            histlist{i} = 'Convert';
        case 'spm_eeg_epochs'
            histlist{i} = 'Epoch';
        case 'spm_eeg_filter'
            histlist{i} = [upper(h(i).args.filter.band(1)) h(i).args.filter.band(2:end)...
                ' filter ' num2str(h(i).args.filter.PHz(:)', '%g %g') ' Hz'];
        case 'spm_eeg_downsample'
            histlist{i}  = ['Downsample to ' num2str(h(i).args.fsample_new) ' Hz'];
        case 'spm_eeg_montage'
            histlist{i} = 'Change montage';
        case 'spm_eeg_artefact'
            histlist{i} = 'Detect artefacts';
        case 'spm_eeg_average'
            histlist{i} = 'Average';
        case 'spm_eeg_average_TF'
            histlist{i} = 'Average time-frequency';
        case 'spm_eeg_average'
            histlist{i} = 'Average';
        case 'spm_eeg_grandmean'
            histlist{i} = 'Grand mean';
        case 'spm_eeg_merge'
            histlist{i} = 'Merge';
        case 'spm_eeg_tf'
            histlist{i} = 'Compute time-frequency';
        case 'spm_eeg_weight_epochs'
            histlist{i} = 'Compute contrast';
        case 'spm_eeg_prep'
            switch h(i).args.task
                case 'settype'
                    histlist{i} = 'Set channel type';
                case {'loadtemplate', 'setcoor2d', 'project3D'}
                    histlist{i} = 'Set 2D coordinates';
                case 'loadeegsens'
                    histlist{i} = 'Load EEG sensor locations';
                case 'defaulteegsens'
                    histlist{i} = 'Set EEG sensor locations to default';
                case 'sens2chan'
                    histlist{i} = 'Specify initial montage';
                case 'headshape'
                    histlist{i} = 'Load fiducials/headshape';
                case 'coregister'
                    histlist{i} = 'Coregister';
                otherwise
                    histlist{i} = ['Prepare: ' h(i).args.task];
            end
        otherwise
            histlist{i} = h(i).fun;
    end
end

switch get(gcbo, 'Label')
    case 'Review history'
        listdlg('ListString', histlist, 'SelectionMode', 'single', 'Name', 'Review history', 'ListSize', [400 300]);
        return
    case 'Save as script'
        [selection ok]= listdlg('ListString', histlist, 'SelectionMode', 'multiple',...
            'InitialValue', 1:numel(histlist) ,'Name', 'Select history entries', 'ListSize', [400 300]);
        if ok
            S = [];
            S.history = h(selection);

            spm_eeg_hist2script(S);
        end
end


