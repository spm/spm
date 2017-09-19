function spm_eeg_prep_ui(callback)
% User interface for spm_eeg_prep function performing several tasks
% for preparation of converted MEEG data for further analysis
% FORMAT spm_eeg_prep_ui(callback)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_prep_ui.m 7169 2017-09-19 10:42:27Z vladimir $


spm('Pointer','Watch');

if ~nargin, callback = 'CreateMenu'; end

eval(callback);

spm('Pointer','Arrow');

end
%==========================================================================
% function CreateMenu
%==========================================================================
function CreateMenu

SVNrev = '$Rev: 7169 $';
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

FileHeaderMenu = uimenu(FileMenu, ...
    'Label','Load header',...
    'Separator','off',...
    'Tag','EEGprepUI',...
    'HandleVisibility', 'on',...
    'Accelerator','H',...
    'Callback', 'spm_eeg_prep_ui(''FileHeaderCB'')');

FileImportMenu = uimenu(FileMenu, ...
    'Label','Import from workspace',...
    'Separator','off',...
    'Tag','EEGprepUI',...
    'HandleVisibility', 'on',...
    'Accelerator','I',...
    'Callback', 'spm_eeg_prep_ui(''FileImportCB'')');

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

% ====== Batch inputs ======================================

BatchInputsMenu = uimenu(Finter,'Label','Batch inputs',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on');

BInputsChannelsMenu = uimenu(BatchInputsMenu, 'Label', 'Channel selection',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''ChannelsCB'')');

BInputsTrialsMenu = uimenu(BatchInputsMenu, 'Label', 'Trial definition',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''TrialsCB'')');

BInputsEventsBIDSMenu = uimenu(BatchInputsMenu, ...
    'Label','Load events from BIDS',...
    'Separator','off',...
    'Enable','off',...
    'Tag','EEGprepUI',...
    'HandleVisibility', 'on',...
    'Callback', 'spm_eeg_prep_ui(''EventsBIDSCB'')');

BInputsEventsMenu = uimenu(BatchInputsMenu, 'Label', 'Event list',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''EventsCB'')');

BInputsMontageMenu = uimenu(BatchInputsMenu, 'Label', 'Montage',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on');

BInputsSortCondMenu = uimenu(BatchInputsMenu, 'Label', 'Sort conditions',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''SortCondCB'')');

MontageCustomMenu = uimenu(BInputsMontageMenu, 'Label', 'Custom Montage',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''MontageCB'')');

MontageRereferenceMenu = uimenu(BInputsMontageMenu, 'Label', 'Re-reference',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''RereferenceCB'')');

MontageROIMenu = uimenu(BInputsMontageMenu, 'Label', 'ROI',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''roiCB'')');

% ====== Channel types ===============================

ChanTypeMenu = uimenu(Finter,'Label','Channel types',...
    'Tag','EEGprepUI',...
    'Enable', 'off', ...
    'HandleVisibility','on');

chanTypes = {'EEG', 'EOG', 'ECG', 'EMG', 'LFP', 'PHYS', 'ILAM', 'SRC', 'Other'};

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

CTypesBIDSMenu = uimenu(ChanTypeMenu, 'Label', 'From BIDS',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...   
    'Callback', 'spm_eeg_prep_ui(''ChanTypeBIDSCB'')');

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

LoadMEGSensMenu = uimenu(Coor3DMenu, 'Label', 'Load MEG sensors',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''LoadMEGSensCB'')');

HeadshapeMenu = uimenu(Coor3DMenu, 'Label', 'Load MEG Fiducials/Headshape',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''HeadshapeCB'')');

LoadEEGSensMenu = uimenu(Coor3DMenu, 'Label', 'Load EEG sensors',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'Separator', 'on', ...
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

ReferenceMenu = uimenu(Coor3DMenu, 'Label', 'Define EEG referencing',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on');

ReferenceSelectMenu = uimenu(ReferenceMenu, 'Label', 'Select reference sensors',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''ReferenceCB'')');

ReferenceMontageMenu = uimenu(ReferenceMenu, 'Label', 'Specify montage matrix',...
    'Tag','EEGprepUI',...
    'Enable', 'on', ...
    'HandleVisibility','on',...
    'Callback', 'spm_eeg_prep_ui(''ReferenceCB'')');

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
end

%==========================================================================
% function FileOpenCB
%==========================================================================
function FileOpenCB

D = spm_eeg_load;
setD(D);

update_menu;

end

%==========================================================================
% function FileHeaderCB
%==========================================================================
function FileHeaderCB

S = [];
S.mode = 'header';
[S.dataset, sts] = spm_select(1, '.*', 'Select M/EEG data file');
if ~sts, return; end
D = spm_eeg_convert(S);
setD(D);

update_menu;

end

%==========================================================================
% function FileImportCB
%==========================================================================
function FileImportCB

var = evalin('base', 'whos');
var = {var(ismember({var.class}, {'single', 'double'})).name};

[selection, ok]= listdlg('ListString', var, 'SelectionMode', 'single' ,'Name', 'Select variable' , 'ListSize', [400 300]);

if ok
    data = full(double(evalin('base', var{selection})));
    dim = size(data);
    ndim = length(dim);
    
    if length(dim)>3
        error('Data of up to 3 dimensions are supported.');
    end
    
    prompt = sprintf('%d|', dim);
    prompt(end) = [];
    
    chandim = spm_input('Number of channels:', 1, prompt, 1:length(dim));
    
    nchan = dim(chandim);
    
    otherdim = setdiff(1:ndim, chandim);
    
    if ndim == 3        
        prompt = sprintf('%d|', dim(otherdim));
        prompt(end) = [];
        trialdim = spm_input('Number of trials:','+1', prompt, otherdim);
        
        ntrials = dim(trialdim);
    else
        trialdim = [];
        ntrials  = 1;
    end
        
    sampledim = setdiff(1:ndim, [trialdim, chandim]);
    nsamples  = dim(sampledim);
    
    data    = permute(data, [chandim, sampledim, trialdim]);
    
    fs      = spm_input('Sampling rate (Hz):', '+1','r');
    t_onset = spm_input('Time axis onset (ms):', '+1','r');
    
    fname   = spm_input('Output file name:', '+1','s');
    
    time    = (0:(nsamples-1))/fs - 1e-3*t_onset;
    
    for i = 1:ntrials
        ftdata.trial{i} = data(:, :, i);
        ftdata.time{i}  = time;
    end
    
    for i = 1:nchan
        ftdata.label{i, 1} = sprintf('Ch%d', i); 
    end
    
    D = spm_eeg_ft2spm(ftdata, fname);
    
    spm_eeg_review(D);
    spm_eeg_review_callbacks('edit', 'prep');    
end
end
%==========================================================================
% function FileSaveCB
%==========================================================================
function FileSaveCB

D = getD;
if ~isempty(D)
    D.save;
end

update_menu;

end

%==========================================================================
% function FileExitCB
%==========================================================================
function FileExitCB

spm_figure('Clear','Interactive');
spm('FigName','M/EEG prepare: done');

end

%==========================================================================
% function ChannelsCB
%==========================================================================
function ChannelsCB

D = getD;

[selection, ok]= listdlg('ListString', D.chanlabels, 'SelectionMode', 'multiple' ,'Name', 'Select channels' , 'ListSize', [400 300]);
if ~ok, return; end

[chanfilename, chanpathname] = uiputfile( ...
    {'*.mat', 'MATLAB File (*.mat)'}, 'Save channel selection as');

label = D.chanlabels(selection);

save(fullfile(chanpathname, chanfilename), 'label', spm_get_defaults('mat.format'));

end

%==========================================================================
% function TrialsCB
%==========================================================================
function TrialsCB

D = getD;

S = [];
S.D = D;
S.save = 1;
spm_eeg_definetrial(S);

CreateMenu;

setD(D);

update_menu;

end

%==========================================================================
% function EventsBIDSCB
%==========================================================================
function EventsBIDSCB

D = getD;

S = [];
S.task = 'loadbidsevents';
S.D = D;
[S.filename, sts] = spm_select(1, '.*_events.tsv$', 'Select BIDS tsv file');

if sts
    S.replace = spm_input('Replace existing','1','replace|add',[1 0], 1);
end

D = spm_eeg_prep(S);

setD(D);

update_menu;

end

%==========================================================================
% function EventsCB
%==========================================================================
function EventsCB

D = getD;

ev = D.events;

if iscell(ev)
    sev = ev{1};
    for i = 2:numel(ev)
        sev = spm_cat_struct(sev, ev{i});
    end
else
    sev = ev;
end

selected = spm_eeg_select_event_ui(sev);

events = [];
for i = 1:size(selected, 1)
    events(i).type  = selected{i, 1};
    events(i).value = selected{i, 2};
end

if ~isempty(events)
    [filename, pathname] = uiputfile('*.mat', 'Save event list as');
    
    if ~isequal(filename, 0)
        save(fullfile(pathname, filename), 'events', spm_get_defaults('mat.format'));
    end
end
end

%==========================================================================
% function MontageCB
%==========================================================================
function MontageCB

D = getD;

label = D.chanlabels(D.indchantype({'EEG', 'EOG', 'ECG', 'EMG'}));

montage = [];
montage.labelorg = label;
montage.labelnew = label;
montage.tra = eye(numel(label));

spm_eeg_montage_ui(montage);

end

%==========================================================================
% function RereferenceCB
%==========================================================================
function RereferenceCB

D = getD;

% Get indices for just EEG channels and remove any bad channels
%--------------------------------------------------------------------------
eegchan  = D.indchantype('EEG');
goodind  = D.indchantype('EEG', 'GOOD');

if isempty(goodind)
    error('No good EEG channels.')
end

%-Get reference channel indices
%--------------------------------------------------------------------------
[selection, ok]= listdlg('ListString', D.chanlabels(goodind), 'SelectionMode', 'multiple' ,...
    'Name', 'Select reference channels' , 'ListSize', [400 300]);
if ~ok
    return;
end

refind  = goodind(selection);
badind  = D.indchantype('EEG', 'BAD');

goodind = find(ismember(eegchan, goodind));
badind  = find(ismember(eegchan,  badind));
refind  = find(ismember(eegchan,  refind));

tra                 = eye(length(eegchan));
tra(goodind,refind) = tra(goodind,refind) - 1/length(refind);
tra(badind, refind) = tra(badind,refind)  - 1/length(refind);

montage          = [];
montage.labelorg = D.chanlabels(eegchan);
montage.labelnew = D.chanlabels(eegchan);
montage.tra      = tra;

[filename, pathname] = uiputfile('*.mat', 'Save montage as');

if ~isequal(filename, 0)
    save(fullfile(pathname, filename), 'montage', spm_get_defaults('mat.format'));
end

end

%==========================================================================
% function roiCB
%==========================================================================
function roiCB

D = getD;

[modality, chanind]  = spm_eeg_modality_ui(D, 0, 1);

chanind = setdiff(chanind, D.badchannels);

if isempty(chanind)
    error(['No good ' modality ' channels.']);
end

montage = [];
montage.labelorg = D.chanlabels(chanind)';
montage.labelnew = {};
montage.tra      = [];
while 1
    roilabel = spm_input('ROI Label:', '+1', 's');
    
    [selection, ok]= listdlg('ListString', D.chanlabels(chanind), 'SelectionMode', 'multiple' ,...
        'Name', 'Select reference channels' , 'ListSize', [400 300]);
    if ok
       montage.labelnew{end+1, 1} = roilabel;
       tra = zeros(1, length(chanind));
       tra(selection) = 1;
       tra = tra./sum(tra);
       montage.tra = [montage.tra; tra];        
    end
    
    if spm_input('Add another?','+1','yes|no',[0, 1], 0);
        break;
    end
end
    
[filename, pathname] = uiputfile('*.mat', 'Save montage as');

if ~isequal(filename, 0)
    save(fullfile(pathname, filename), 'montage', spm_get_defaults('mat.format'));
end

update_menu;

end

%==========================================================================
% function SortCondCB
%==========================================================================
function SortCondCB

D = getD;

S   = [];
S.task = 'sortconditions';
S.D = D;
oldcondlist = D.condlist;
S.condlist = cell(size(oldcondlist));
for i = 1:D.nconditions
    str = sprintf('%s|', oldcondlist{:});
    str = str(1:(end-1));
    
    ind = spm_input(['Select condition ' num2str(i)], 1, 'm', str, 1:numel(oldcondlist));
    S.condlist(i) = oldcondlist(ind);
    oldcondlist(ind) = [];
end

[filename, pathname] = uiputfile('*.mat', 'Save conditions list as');

if ~isequal(filename, 0)
    condlist = S.condlist;
    
    save(fullfile(pathname, filename), 'condlist', spm_get_defaults('mat.format'));
end

D = spm_eeg_prep(S);

setD(D);
update_menu;

end

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
        
        [selection,ok]= listdlg('ListString', chanlist, 'SelectionMode', 'multiple',...
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

end
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

end
%==========================================================================
% function ChanTypeDefaultCB
%==========================================================================
function ChanTypeDefaultCB

S.D    = getD;
S.task = 'defaulttype';
D      = spm_eeg_prep(S);

setD(D);
update_menu;

end

%==========================================================================
% function ChanTypeBIDSCB
%==========================================================================
function ChanTypeBIDSCB

S = [];
S.D    = getD;
S.task = 'bidschantype';
[S.filename, sts] = spm_select(1, '.*_channels.tsv$', 'Select BIDS tsv file');

if sts
    D      = spm_eeg_prep(S);
    
    setD(D);
    update_menu;
end

end

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

end
%==========================================================================
% function LoadEEGSensCB
%==========================================================================
function LoadMEGSensCB

D = getD;

S = [];
S.task = 'loadmegsens';
S.D = D;
[S.source,sts] = spm_select(1, '\.*', 'Select a raw MEG data file');
if ~sts, return; end

D = spm_eeg_prep(S);

setD(D);

update_menu;

end
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
                [lblshape{fidnum},fidlabel] = strtok(fidlabel);
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

setD(D);

update_menu;

end
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

end
%==========================================================================
% function CoregisterCB
%==========================================================================
function ReferenceCB

D = getD;

S = [];
S.D = D;
S.task = 'sens2chan';
switch get(gcbo, 'Label')
    case 'Select reference sensors'
        elec = D.sensors('EEG');
        
        [refind, ok]= listdlg('ListString', elec.label, 'SelectionMode', 'multiple' ,'Name', 'Select reference electrodes' , 'ListSize', [400 300]);
        if ~ok
            return;
        end
        
        S.refelec = elec.label(refind);
    case 'Specify montage matrix'
        
        sens  = D.sensors('EEG');
        label = D.chanlabels(D.indchantype('EEG'));
        
        [sel1, sel2] = spm_match_str(label, sens.label);
        
        montage = [];
        montage.labelorg = sens.label;
        montage.labelnew = label;
        montage.tra = zeros(numel(label), numel(sens.label));
        montage.tra(sub2ind(size(montage.tra), sel1, sel2)) = 1;
        
        S.montage = spm_eeg_montage_ui(montage);
end

D = spm_eeg_prep(S);

setD(D);

update_menu;

end
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

end
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

end
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

end
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

save(fullfile(pathname, filename), 'Cnames', 'Cpos', 'Rxy', 'Nchannels', spm_get_defaults('mat.format'));

end
%==========================================================================
% function Project3DCB
%==========================================================================
function Project3DCB

D = getD;

switch get(gcbo, 'Label')
    case 'Project 3D (EEG)'
        modality = 'EEG';
        chanind  = D.indchantype('EEG');
    case 'Project 3D (MEG)'
        modality = 'MEG';
        chanind  = D.indchantype('MEGANY');
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

[sel1, sel2] = spm_match_str(D.chanlabels(chanind), label);

plot_sensors2D(xy(:, sel2), label(sel2));

update_menu;

end
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

end
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

end
%==========================================================================
% function update_menu
%==========================================================================
function update_menu

Finter = spm_figure('GetWin','Interactive');
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'File'), 'Enable', 'on');

IsEEG = 'off';
IsMEG = 'off';
IsEpochable = 'off';
HasEvents   = 'off';
IsEpoched   = 'off';
HasPlanar   = 'off';
HasSensors  = 'off';
HasSensorsEEG = 'off';
ReferenceSelectable = 'off';
HasSensorsMEG = 'off';
HasChannelsMEGREF = 'off';
HasFiducials = 'off';
HasDefaultLocs = 'off';
HasHistory = 'off';
if isa(get(Finter, 'UserData'), 'meeg')
    Dloaded = 'on';
    
    D = getD;
    
    if ~isempty(D.indchantype('EEG'))
        IsEEG = 'on';
    end
    
    if ~isempty(D.indchantype('MEG'));
        IsMEG = 'on';
    end
    
    if ~isempty(D.events)
        HasEvents = 'on';
    end
    
    if isequal(D.type, 'continuous') && ~isempty(D.events);
        IsEpochable = 'on';
    end
    
    if ~isequal(D.type, 'continuous')
        IsEpoched = 'on';
    end
    
    
    [res, list] = modality(D, 1, 1);
    if ismember('MEGPLANAR', list);
        HasPlanar = 'on';
    end
    
    if ~isempty(D.indchantype('REF'));
        HasChannelsMEGREF = 'on';
    end
    
    if ~isempty(D.sensors('EEG')) || ~isempty(D.sensors('MEG'))
        HasSensors = 'on';
    end
    
    if ~isempty(D.sensors('EEG'))
        HasSensorsEEG = 'on';
        
        elec = D.sensors('EEG');
        if ~isfield(elec, 'tra') &&  (size(elec.elecpos, 1) == numel(elec.label))
            ReferenceSelectable = 'on';
        end
    end
    
    if  ~isempty(D.sensors('MEG'))
        HasSensorsMEG = 'on';
    end
    
    if  ~isempty(D.fiducials)
        HasFiducials = 'on';
    end
    
    template_sfp = dir(fullfile(spm('dir'), 'EEGtemplates', '*.sfp'));
    template_sfp = {template_sfp.name};
    ind = strmatch([ft_senstype(D.chanlabels(D.indchantype('EEG'))) '.sfp'], template_sfp, 'exact');
    
    if ~isempty(ind) || ft_senstype(D.chanlabels(D.indchantype('EEG')), 'ext1020')
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
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load events from BIDS'), 'Enable', Dloaded);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Batch inputs'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Trial definition'), 'Enable', IsEpochable);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Event list'), 'Enable', HasEvents);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Montage'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Re-reference'), 'Enable', IsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'ROI'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Sort conditions'), 'Enable', IsEpoched);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Channel types'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Sensors'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', '2D projection'), 'Enable', Dloaded);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'MEGREF=>MEG'), 'Enable', HasChannelsMEGREF);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Assign default'), 'Enable', HasDefaultLocs);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load EEG sensors'), 'Enable', IsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load MEG sensors'), 'Enable', IsMEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load MEG Fiducials/Headshape'), 'Enable', HasSensorsMEG);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Headshape'), 'Enable', HasSensorsMEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Define EEG referencing'), 'Enable', HasSensorsEEG);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Select reference sensors'), 'Enable', ReferenceSelectable);
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

end
%==========================================================================
% function getD
%==========================================================================
function D = getD

Finter = spm_figure('GetWin','Interactive');
D = get(Finter, 'UserData');
if ~isa(D, 'meeg')
    D = [];
end

end
%==========================================================================
% function setD
%==========================================================================
function setD(D)
Finter = spm_figure('GetWin','Interactive');
set(Finter, 'UserData', D);

end
%==========================================================================
% function getHandles
%==========================================================================
function handles = getHandles
Fgraph = spm_figure('GetWin','Graphics');
handles = get(Fgraph, 'UserData');

end
%==========================================================================
% function setHandles
%==========================================================================
function setHandles(handles)
Fgraph = spm_figure('GetWin','Graphics');
set(Fgraph, 'UserData', handles);

end
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
    
    handles.h_el = [];
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

end
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

end
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

end
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

end
%==========================================================================
% function Clear2DCB
%==========================================================================
function Clear2DCB

plot_sensors2D([], {});

update_menu;

end
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
            [selection,ok]= listdlg('ListString',lblshape, 'SelectionMode', 'single',...
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
end
