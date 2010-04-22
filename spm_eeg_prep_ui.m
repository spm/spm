function spm_eeg_prep_ui(callback)
% User interface for spm_eeg_prep function performing several tasks
% for preparation of converted MEEG data for further analysis
% FORMAT spm_eeg_prep_ui(callback)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_prep_ui.m 3833 2010-04-22 14:49:48Z vladimir $


spm('Pointer','Watch');

if ~nargin, callback = 'CreateMenu'; end

eval(callback);

spm('Pointer','Arrow');


%==========================================================================
% function CreateMenu
%==========================================================================
function CreateMenu
    
SVNrev = '$Rev: 3833 $';
spm('FnBanner', 'spm_eeg_prep_ui', SVNrev);
Finter = spm('FnUIsetup', 'M/EEG prepare', 0);

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
    'Accelerator','O',...
    'Callback', 'spm_eeg_prep_ui(''FileOpenCB'')');

FileSaveMenu = uimenu(FileMenu, ...
    'Label','Save',...
    'Separator','off',...
    'Tag','EEGprepUI',...
    'Enable','off',...
    'HandleVisibility', 'on',...
    'Accelerator','S',...
    'Callback', 'spm_eeg_prep_ui(''FileSaveCB'')');

FileExitMenu = uimenu(FileMenu, ...
    'Label','Quit',...
    'Separator','on',...
    'Tag','EEGprepUI',...
    'HandleVisibility', 'on',...
    'Accelerator','Q',...
    'Callback', 'spm_eeg_prep_ui(''FileExitCB'')');

% ====== Channel types ===============================

ChanTypeMenu = uimenu(Finter,'Label','Channel types',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on');

chanTypes = {'EEG', 'EOG', 'ECG', 'EMG', 'LFP', 'Other'};

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


%==========================================================================
% function FileOpenCB
%==========================================================================
function FileOpenCB

D = spm_eeg_load;
setD(D);

update_menu;


%==========================================================================
% function FileSaveCB
%==========================================================================
function FileSaveCB

D = getD;
if ~isempty(D)
    D.save;
end

update_menu;


%==========================================================================
% function FileExitCB
%==========================================================================
function FileExitCB

spm_figure('Clear','Interactive');
spm('FigName','M/EEG prepare: done');


%==========================================================================
% function ChanTypeCB
%==========================================================================
function ChanTypeCB

type = get(gcbo, 'Label');

D = getD;

if ~isempty(D)
    chanlist ={};
    for i = 1:D.nchannels
        if strncmp(D.chantype(i), 'MEG', 3) || strncmp(D.chantype(i), 'REF', 3)
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


%==========================================================================
% function MEGChanTypeCB
%==========================================================================
function MEGChanTypeCB

S = [];
S.D = getD;
S.task = 'settype';

switch get(gcbo, 'Label')
    case 'MEGREF=>MEG'
        dictionary = {
            'REFMAG',   'MEGMAG';
            'REFGRAD',  'MEGGRAD';
            };               
                
        ind = spm_match_str(S.D.chantype, dictionary(:,1));

        grad = S.D.sensors('meg');
        if ~isempty(grad)
            % Under some montages only subset of the reference sensors are
            % in the grad
            [junk, sel] = intersect(S.D.chanlabels(ind), grad.label);
            ind = ind(sel);
        end
            
        S.ind = ind;
        
        [sel1, sel2] = spm_match_str(S.D.chantype(S.ind), dictionary(:, 1));

        S.type = dictionary(sel2, 2);

        D = spm_eeg_prep(S);
end

setD(D);
update_menu;


%==========================================================================
% function ChanTypeDefaultCB
%==========================================================================
function ChanTypeDefaultCB

S.D    = getD;
S.task = 'defaulttype';
D      = spm_eeg_prep(S);

setD(D);
update_menu;


%==========================================================================
% function LoadEEGSensTemplateCB
%==========================================================================
function LoadEEGSensTemplateCB

S.D    = getD;
S.task = 'defaulteegsens';

if strcmp(S.D.modality(1, 0), 'Multimodal')
    fid = fiducials(S.D);
    if ~isempty(fid)
        lblfid = fid.fid.label;

        S.regfid = match_fiducials({'nas'; 'lpa'; 'rpa'}, lblfid);
        S.regfid(:, 2) = {'spmnas'; 'spmlpa'; 'spmrpa'};
    else
        warndlg(strvcat('Could not match EEG fiducials for multimodal dataset.', ...
            '           EEG coregistration might fail.'));
    end
end

D = spm_eeg_prep(S);
setD(D);
update_menu;


%==========================================================================
% function LoadEEGSensCB
%==========================================================================
function LoadEEGSensCB

D = getD;

switch get(gcbo, 'Label')
    case 'From *.mat file'
        [S.sensfile, sts] = spm_select(1,'mat','Select EEG sensors file');
        if ~sts, return, end
        S.source = 'mat';
        [S.headshapefile, sts] = spm_select(1,'mat','Select EEG fiducials file');
        if ~sts, return, end
        S.fidlabel = spm_input('Fiducial labels:', '+1', 's', 'nas lpa rpa');
    case 'Convert locations file'
        [S.sensfile, sts] = spm_select(1, '.*', 'Select locations file');
        if ~sts, return, end
        S.source = 'locfile';
end

if strcmp(D.modality(1, 0), 'Multimodal')   
    if ~isempty(D.fiducials)
        S.regfid = {};
        if strcmp(S.source, 'mat')
            fidlabel = S.fidlabel;
            lblshape = {};
            fidnum = 0;
            while ~all(isspace(fidlabel))
                fidnum = fidnum+1;
                [lblshape{fidnum} fidlabel] = strtok(fidlabel);
            end
            if (fidnum < 3)
                error('At least 3 labeled fiducials are necessary');
            end
        else
            shape = ft_read_headshape(S.sensfile);
            lblshape = shape.fid.label;
        end

        fid = fiducials(D);
        lblfid = fid.fid.label;

        S.regfid = match_fiducials(lblshape, lblfid);
    else
        warndlg(strvcat('Could not match EEG fiducials for multimodal dataset.', ...
            '           EEG coregistration might fail.'));
    end
end

S.D    = D;
S.task = 'loadeegsens';
D      = spm_eeg_prep(S);

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

setD(D);

update_menu;


%==========================================================================
% function HeadshapeCB
%==========================================================================
function HeadshapeCB

S = [];
S.D = getD;
S.task = 'headshape';

[S.headshapefile, sts] = spm_select(1, '.*', 'Select fiducials/headshape file');
if ~sts, return, end
S.source = 'convert';

shape = ft_read_headshape(S.headshapefile);
lblshape = shape.fid.label;

fid = fiducials(S.D);

if ~isempty(fid)
    lblfid = fid.fid.label;
    S.regfid = match_fiducials(lblshape, lblfid);
end

D = spm_eeg_prep(S);

setD(D);

update_menu;


%==========================================================================
% function CoregisterCB
%==========================================================================
function CoregisterCB

S = [];
S.D = getD;
S.task = 'coregister';

D = spm_eeg_prep(S);

% Bring the menu back
spm_eeg_prep_ui;

setD(D);

update_menu;


%==========================================================================
% function EditExistingCoor2DCB
%==========================================================================
function EditExistingCoor2DCB

D = getD;

switch get(gcbo, 'Label')
    case 'Edit existing MEG'
        xy = D.coor2D('MEG');
        label = D.chanlabels(strmatch('MEG', D.chantype));
    case 'Edit existing EEG'
        xy = D.coor2D('EEG');
        label = D.chanlabels(strmatch('EEG', D.chantype, 'exact'));
end

plot_sensors2D(xy, label);

update_menu;


%==========================================================================
% function LoadTemplateCB
%==========================================================================
function LoadTemplateCB

[sensorfile, sts] = spm_select(1, 'mat', 'Select sensor template file', ...
    [], fullfile(spm('dir'), 'EEGtemplates'));
if ~sts, return, end

template = load(sensorfile);

if isfield(template, 'Cnames') && isfield(template, 'Cpos')
    plot_sensors2D(template.Cpos, template.Cnames);
end

update_menu;


%==========================================================================
% function SaveTemplateCB
%==========================================================================
function SaveTemplateCB

handles=getHandles;

Cnames = handles.label;
Cpos = handles.xy;
Rxy = 1.5;
Nchannels = length(Cnames);

[filename, pathname] = uiputfile('*.mat', 'Save channel template as');

save(fullfile(pathname, filename), 'Cnames', 'Cpos', 'Rxy', 'Nchannels');


%==========================================================================
% function Project3DCB
%==========================================================================
function Project3DCB

D = getD;

switch get(gcbo, 'Label')
    case 'Project 3D (EEG)'
        modality = 'EEG';
    case 'Project 3D (MEG)'
        modality = 'MEG';
end

if ~isfield(D, 'val')
    D.val = 1;
end

if isfield(D, 'inv') && isfield(D.inv{D.val}, 'datareg')
    datareg = D.inv{D.val}.datareg;
    ind     = strmatch(modality, {datareg(:).modality}, 'exact');
    sens    = datareg(ind).sensors;
else
    sens    = D.sensors(modality);
end

[xy, label] = spm_eeg_project3D(sens, modality);

plot_sensors2D(xy, label);

update_menu;


%==========================================================================
% function AddCoor2DCB
%==========================================================================
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


%==========================================================================
% function ApplyCoor2DCB
%==========================================================================
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


%==========================================================================
% function update_menu
%==========================================================================
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

    if ~isempty(strmatch('MEG', D.chantype));
        IsMEG = 'on';
    end

    if ft_senstype(D.chanlabels, 'neuromag') &&...
            isfield(D, 'origchantypes')
        IsNeuromag = 'on';
    end

    if ~isempty(strmatch('REF', D.chantype));
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
    ind = strmatch([ft_senstype(D.chanlabels) '.sfp'], template_sfp, 'exact');

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

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Save'), 'Enable', 'on');

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Channel types'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Sensors'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', '2D projection'), 'Enable', Dloaded);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'MEGREF=>MEG'), 'Enable', HasChannelsMEGREF);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Assign default'), 'Enable', HasDefaultLocs);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load EEG sensors'), 'Enable', IsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load MEG Fiducials/Headshape'), 'Enable', HasSensorsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Headshape'), 'Enable', HasSensorsMEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Coregister'), 'Enable', HasSensors);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing EEG'), 'Enable', IsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Edit existing MEG'), 'Enable', IsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (EEG)'), 'Enable', HasSensorsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Project 3D (MEG)'), 'Enable', HasSensorsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Delete sensor'), 'Enable', IsSelected);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Undo move'), 'Enable', IsMoved);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Apply'), 'Enable', IsTemplate);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Clear'), 'Enable', IsTemplate);

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
        'tooltipstring','Send changes to ''SPM Graphics'' window',...
        'BusyAction','cancel',...
        'Interruptible','off',...
        'Tag','EEGprepUI');
end


figure(Finter);


%==========================================================================
% function getD
%==========================================================================
function D = getD

Finter = spm_figure('GetWin','Interactive');
D = get(Finter, 'UserData');
if ~isa(D, 'meeg')
    D = [];
end


%==========================================================================
% function setD
%==========================================================================
function setD(D)
Finter = spm_figure('GetWin','Interactive');
set(Finter, 'UserData', D);


%==========================================================================
% function getHandles
%==========================================================================
function handles = getHandles
Fgraph = spm_figure('GetWin','Graphics');
handles = get(Fgraph, 'UserData');


%==========================================================================
% function setHandles
%==========================================================================
function setHandles(handles)
Fgraph = spm_figure('GetWin','Graphics');
set(Fgraph, 'UserData', handles);


%==========================================================================
% function plot_sensors2D
%==========================================================================
function plot_sensors2D(xy, label)

Fgraph = spm_figure('GetWin','Graphics');

spm_clf(Fgraph);

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

    handles.TemplateFrame = ...
    plot([0.05 0.05 0.95 0.95 0.05], [0.05 0.95 0.95 0.05 0.05], 'k-');
    
    axis off;

end

handles.xy = xy;
handles.label = label(:)';

setHandles(handles);

update_menu;


%==========================================================================
% function DeleteCoor2DCB
%==========================================================================
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


%==========================================================================
% function UndoMoveCoor2DCB
%==========================================================================
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


%==========================================================================
% function LabelClickCB
%==========================================================================
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


%==========================================================================
% function LabelMoveCB
%==========================================================================
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


%==========================================================================
% function CancelMoveCB
%==========================================================================
function CancelMoveCB

Fgraph = spm_figure('GetWin','Graphics');
handles = getHandles;

handles.pointSelected=[];
handles.labelSelected=[];
set(Fgraph, 'WindowButtonMotionFcn', '');

setHandles(handles);

update_menu;


%==========================================================================
% function Clear2DCB
%==========================================================================
function Clear2DCB

plot_sensors2D([], {});

update_menu;


%==========================================================================
% function match_fiducials
%==========================================================================
function regfid = match_fiducials(lblshape, lblfid)

if numel(intersect(upper(lblshape), upper(lblfid))) < 3
    if numel(lblshape)<3 || numel(lblfid)<3
        warndlg('3 fiducials are required');
        return;
    else
        regfid = {};
        for i = 1:length(lblfid)
            [selection ok]= listdlg('ListString',lblshape, 'SelectionMode', 'single',...
                'InitialValue', strmatch(upper(lblfid{i}), upper(lblshape)), ...
                'Name', ['Select matching fiducial for ' lblfid{i}], 'ListSize', [400 300]);

            if ~ok
                continue
            end

            regfid = [regfid; [lblfid(i) lblshape(selection)]];
        end

        if size(regfid, 1) < 3
            warndlg('3 fiducials are required to load headshape');
            return;
        end
    end
else
    [sel1, sel2] = spm_match_str(upper(lblfid), upper(lblshape));
    lblfid = lblfid(sel1);
    lblshape = lblshape(sel2);
    regfid = [lblfid(:) lblshape(:)];
end
