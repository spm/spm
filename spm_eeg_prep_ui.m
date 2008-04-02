function spm_eeg_prep_ui(callback)

if nargin == 0

    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','MEEG preparation',0);

    delete(findobj(get(Finter,'Children'),'Tag','EEGprepUI'))

    %-Draw top level menu
    %-----------------------------------------------------------------------
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

    Coor2DMenu = uimenu(Finter,'Label','2D projection',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on');

    LoadTemplateMenu = uimenu(Coor2DMenu, 'Label', 'Load template',...
        'Tag','EEGprepUI',...
        'Enable', 'on', ...
        'HandleVisibility','on',...
        'Callback', 'spm_eeg_prep_ui(''LoadTemplateCB'')');

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

function LoadTemplateCB

S=[];
S.task = 'loadtemplate';
S.P = spm_select(1, '\.mat$', 'Select sensor template file', ...
    [], fullfile(spm('dir'), 'EEGtemplates'));
S.D = getD;
D = spm_eeg_prep(S);
setD(D);
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
function update_menu

Finter = spm_figure('GetWin','Interactive');
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'File'), 'Enable', 'on');

if isa(get(Finter, 'UserData'), 'meeg')
    Dloaded = 'on';

    D = getD;

    IsEEG = 'off';
    if ~isempty(strmatch('EEG', D.chantype, 'exact'))
        IsEEG = 'on';
    end

    IsMEG = 'off';
    if ~isempty(strmatch('MEG', D.chantype, 'exact'));
        IsMEG = 'on';
    end
else
    Dloaded = 'off';
end

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Channel types'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Sensors'), 'Enable', Dloaded);
set(findobj(Finter,'Tag','EEGprepUI', 'Label', '2D projection'), 'Enable', Dloaded);

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Load EEG sensor coordinates'), 'Enable', IsEEG);


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
