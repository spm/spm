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

    Coor2DMenu = uimenu(Finter,'Label','2D projection',...
        'Tag','EEGprepUI',...
        'Enable', 'off', ...
        'HandleVisibility','on');

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
        if strcmp(D.chantype(i), 'MEG')
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
        selection(strmatch('MEG', chantype(D, selection), 'exact')) = [];
        if ok && ~isempty(selection)
            D = chantype(D, selection, type);
            setD(D);
        end
    end
end

update_menu;

%-----------------------------------------------------------------------

function update_menu
Finter = spm_figure('GetWin','Interactive');
set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'File'), 'Enable', 'on');

if isa(get(Finter, 'UserData'), 'meeg')
    Dloaded = 'on';
else
    Dloaded = 'off';
end

set(findobj(Finter,'Tag','EEGprepUI', 'Label', 'Channel types'), 'Enable', Dloaded);

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
